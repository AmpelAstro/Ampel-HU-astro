import logging
import re
import subprocess
import time
import warnings
from argparse import ArgumentParser
from collections.abc import Generator
from pathlib import Path

import numpy as np
import requests
from astropy.table import Table
from astropy.utils.exceptions import AstropyWarning
from pymongo import MongoClient
from tqdm import tqdm

logger = logging.getLogger(__name__ if __name__ != "__main__" else "lsdr10_mongodump")
BASE_URL = "https://portal.nersc.gov/cfs/cosmo/data/legacysurvey/dr{DR}/south/sweep/"


def get_data_dir(data_dir: Path, dr: int):
    return data_dir / f"legacy_survey_dr{dr}"


def get_filenames(data_dir: Path, dr: int, sv: int) -> list[tuple[str, str]]:
    data_dir = get_data_dir(data_dir, dr)
    cache_file = data_dir / "lc_filenames.txt"
    version = f"{dr}.{sv}"
    if not cache_file.exists():
        url = BASE_URL.format(DR=dr) + version + "/"
        logger.debug(f"fetching filenames from {url}")
        response = requests.get(url)
        response.raise_for_status()

        filenames = sorted(
            list(
                set(
                    re.findall(
                        r"sweep-\d+[mp]\d+-\d+[mp]\d+\.fits",
                        response.content.decode(),
                    )
                )
            )
        )

        logger.debug(f"caching filenames to {cache_file}")
        cache_file.parent.mkdir(parents=True, exist_ok=True)
        with open(cache_file, "w") as f:
            for filename in filenames:
                f.write(f"{filename}\n")

    logger.debug(f"loading lightcurve filenames from {cache_file}")
    with open(cache_file) as f:
        filenames = [line.strip() for line in f.readlines()]

    logger.debug(f"found {len(filenames)} sweep filenames")

    return [
        (
            version + "/" + filename,
            version + "-photo-z/" + filename.replace(".fits", "-pz.fits"),
        )
        for filename in filenames
    ]


def get_local_path(data_dir: Path, filename: str, dr: int) -> Path:
    local_filename = get_data_dir(data_dir, dr) / filename
    local_filename.parent.mkdir(parents=True, exist_ok=True)
    return local_filename


def iter_filenames(
    data_dir: Path, i: int | list[int], dr: int, sv: int
) -> Generator[tuple[tuple[str, str], tuple[Path, ...]], None, None]:
    filenames = get_filenames(data_dir, dr, sv)
    if isinstance(i, int):
        i = [i]
    for index in i:
        if index < 0 or index >= len(filenames):
            raise IndexError("Index out of range")
        i_filenames = filenames[index]
        sub_paths: tuple[Path, ...] = tuple(
            [get_local_path(data_dir, fn, dr) for fn in i_filenames]
        )  # type: ignore
        yield i_filenames, sub_paths


def download_file_by_index(
    data_dir: Path, i: int | list[int], dr: int, sv: int
) -> tuple[Path, ...] | list[tuple[Path, Path]]:
    paths: list[tuple[Path, ...]] = []
    for i_filenames, sub_paths in iter_filenames(data_dir, i, dr, sv):
        for sub_filename, local_path in zip(i_filenames, sub_paths, strict=True):
            if not local_path.exists():
                url = BASE_URL.format(DR=dr) + sub_filename
                logger.info(f"Downloading {url} to {local_path}...")
                response = requests.get(url, stream=True)
                response.raise_for_status()
                with open(local_path, "wb") as f:
                    for chunk in tqdm(
                        response.iter_content(chunk_size=int(2**20)),
                        desc=f"Downloading {sub_filename}",
                        unit="MB",
                        unit_scale=True,
                    ):
                        f.write(chunk)

        paths.append(tuple(sub_paths))

    return paths[0] if len(paths) == 1 else paths


def get_aggregation_pipeline(target_col: str, merged_col: str):
    return [
        {
            "$lookup": {
                "from": target_col,
                "let": {
                    "release": "$RELEASE",
                    "brickid": "$BRICKID",
                    "objid": "$OBJID",
                },
                "pipeline": [
                    {
                        "$match": {
                            "$expr": {
                                "$and": [
                                    {"$eq": ["$RELEASE", "$$release"]},
                                    {"$eq": ["$BRICKID", "$$brickid"]},
                                    {"$eq": ["$OBJID", "$$objid"]},
                                ]
                            }
                        }
                    }
                ],
                "as": "b_docs",
            }
        },
        {"$unwind": {"path": "$b_docs", "preserveNullAndEmptyArrays": True}},
        {"$replaceRoot": {"newRoot": {"$mergeObjects": ["$$ROOT", "$b_docs"]}}},
        {"$project": {"b_docs": 0}},
        {
            "$merge": {
                "into": merged_col,
                "whenMatched": "replace",
                "whenNotMatched": "insert",
            }
        },
    ]


INDEX = [("RELEASE", 1), ("BRICKID", 1), ("OBJID", 1)]

USE_COLUMNS = [
    # index columns
    "RELEASE",
    "BRICKID",
    "OBJID",
    # position
    "RA",
    "DEC",
    # morphology
    "SHAPE_R",
    "SHAPE_E1",
    "SHAPE_E2",
    "SERSIC",
    "TYPE",
    # Tractor fluxes
    "FLUX_G",
    "FLUX_R",
    "FLUX_I",
    "FLUX_Z",
    # WISE fluxes
    "FLUX_W1",
    "FLUX_W2",
    "FLUX_W3",
    "FLUX_W4",
    # redshift
    # "Z_SPEC",
    # "SURVEY",
    # "Z_PHOT_MEDIAN",
    "Z_PHOT_MEAN_I",
    # "Z_PHOT_L68",
    # "Z_PHOT_U68",
    "Z_PHOT_STD",
]


def mongodumps_by_index(
    data_dir: Path, i: int | list[int], dr: int, sv: int, mogno_uri: str, db_name: str
):
    client = MongoClient(mogno_uri)
    db = client[db_name]
    out_dir = get_data_dir(data_dir, dr) / f"{dr}.{sv}-merged"
    for _, local_paths in tqdm(
        iter_filenames(data_dir, i, dr, sv), desc="making mongo dumps"
    ):
        if len(local_paths) != 2:
            raise NotImplementedError("Only merging of two files is implemented!")

        merged_collection_name = local_paths[0].stem + "-merged"
        collection_names = []
        export_fn = out_dir / local_paths[0].with_suffix(".bson.gz").name

        if not export_fn.exists():
            logger.info(f"Making mongo dumps for {merged_collection_name}...")
            for local_path in local_paths:
                temp_csv = local_path.with_suffix(".csv.temp")
                logger.debug(f"Saving {local_path} to {temp_csv}")
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore", category=AstropyWarning)
                    t = Table.read(local_path)

                # remove unused columns
                remove_cols = [c for c in t.colnames if c not in USE_COLUMNS]
                t.remove_columns(remove_cols)

                if not temp_csv.exists():  # TODO: remove! only testing
                    # write to temporary csv for mongoimport
                    t.write(temp_csv, format="ascii.csv")

                collection_name = local_path.stem
                if collection_name not in db.list_collection_names():  # TODO: remove!
                    import_cmd = [
                        "mongoimport",
                        f'--uri="{mogno_uri}"',
                        f"--db={db_name}",
                        f"--collection={collection_name}",
                        "--type=csv",
                        "--headerline",
                        f"--file={temp_csv!s}",
                    ]
                    logger.debug("running mongoimport")
                    subprocess.run(import_cmd, check=True)
                    logger.debug(f"done, removing {temp_csv}")
                    temp_csv.unlink()
                collection_names.append(collection_name)

            sorted_collection_names = sorted(collection_names, key=lambda x: len(x))

            logger.debug("creating indices")
            for n in sorted_collection_names:
                db[n].create_index(INDEX)

            logger.debug("merging collections")
            db[sorted_collection_names[0]].aggregate(
                get_aggregation_pipeline(
                    sorted_collection_names[1], merged_collection_name
                )
            )

            logger.debug(f"Dumping {merged_collection_name} to {export_fn}")
            export_cmd = [
                "mongodump",
                f'--uri="{mogno_uri}"',
                f"--db={db_name}",
                f"--collection={merged_collection_name}",
                f"-o={export_fn}",
                "--gzip",
            ]
            subprocess.run(export_cmd, check=True)
            logger.debug("done")

            for c in [*collection_names, merged_collection_name]:
                logger.debug(f"dropping collection {c}")
                db.drop_collection(c)


if __name__ == "__main__":
    parser = ArgumentParser(
        description="Download Legacy Survey sweep FITS files by index."
    )
    parser.add_argument(
        "dir",
        type=Path,
        help="Directory to download files to",
    )
    parser.add_argument(
        "indices",
        metavar="N",
        type=int,
        nargs="+",
        help="Indices of the lightcurve files to download.",
    )
    parser.add_argument("dr", type=int, help="data release number")
    parser.add_argument("sv", type=int, help="subversion of data release")
    parser.add_argument(
        "--skip-download", action="store_true", help="skip downloading files"
    )
    parser.add_argument(
        "--mongo-uri",
        type=str,
        help="mongo uri, if not given no mongo action will be executed",
    )
    parser.add_argument(
        "--db-name", type=str, help="mongo db name, required for mongo actions"
    )
    parser.add_argument("--log-level", default="INFO", help="Set the logging level")
    args = parser.parse_args()
    logging.basicConfig(level=args.log_level.upper())
    logging.getLogger("pymongo").setLevel("INFO")

    start = time.time()
    downloaded_files = download_file_by_index(args.dir, args.indices, args.dr, args.sv)
    for file in np.atleast_1d(downloaded_files):
        logger.info(f"Downloaded: {file}")

    if args.mongo_uri:
        mongodumps_by_index(
            args.dir, args.indices, args.dr, args.sv, args.mongo_uri, args.db_name
        )

    end = time.time()
    logger.info(f"Took {end - start} seconds")
