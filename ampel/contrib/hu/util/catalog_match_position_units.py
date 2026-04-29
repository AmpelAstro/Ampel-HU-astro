import json
from pathlib import Path


def get_catalog_position_unit_map() -> dict:
    """
    The json file can be produced with

    curl -s "https://ampel.zeuthen.desy.de/api/catalogmatch/catalogs/" \
      -H "accept: application/json" \
    | jq '
    [
      .[]
      | {
          name,
          ra: [.columns[] | select(
            (.name | ascii_downcase == "ra") or
            (.name | ascii_downcase == "ramean") or
            (.name | ascii_downcase == "raj2000") or
            (.name | ascii_downcase == "radeg") or
            (.name | ascii_downcase == "ra_central")
          )],
          dec: [.columns[] | select(
            (.name | ascii_downcase == "dec") or
            (.name | ascii_downcase == "decmean") or
            (.name | ascii_downcase == "dej2000") or
            (.name | ascii_downcase == "dedeg") or
            (.name | ascii_downcase == "dec_central")
          )]
        }
    ]' > catalog_position_columns.json
    """
    filename = Path(__file__).parent / "catalog_position_columns.json"
    with filename.open() as f:
        mapping = json.load(f)

    # exclude catalogs with multiple or no columns
    mapping = {
        v["name"]: {"ra": v["ra"][0], "dec": v["dec"][0]}
        for v in mapping
        if len(v["ra"]) == 1 and len(v["dec"]) == 1
    }

    # check sanity
    assert all([v["ra"]["unit"] == v["dec"]["unit"] for v in mapping.values()])

    # add NEDz
    mapping["NEDz"] = mapping["NEDz_extcats"]

    return mapping
