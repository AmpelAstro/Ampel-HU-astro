import logging
from contextlib import suppress

import numpy
from astropy.time import Time
from pymongo import GEOSPHERE, MongoClient, UpdateOne
from pymongo.errors import BulkWriteError, OperationFailure

from ampel.contrib.hu.t3.tns.TNSName import TNSName

log = logging.getLogger(__name__)


class TNSMirrorDB:
    def __init__(self, mongo: MongoClient, db_name="TNS", logger=None) -> None:
        db = mongo.get_database(db_name)
        self.collection = db.get_collection("srcs")
        meta = db.get_collection("meta")
        entry = {
            "_id": "pos",
            "key": "pos",
            "is_indexed": "true",
            "pos_format": "geoJSON",
            "type": "sphere2d",
        }
        try:
            meta.update_one({}, {"$set": entry}, upsert=True)
            self.collection.create_index([("pos", GEOSPHERE)])
        except OperationFailure:
            ...
        if logger:
            self.logger = logger
        else:
            self.logger = log

    def add_sources(self, docs):
        """
        :param docs: TNS get/object results
        """
        # NB: bulk_write() an empty list raises InvalidOperation
        if not docs:
            return

        ops = [
            UpdateOne({"_id": doc["_id"]}, {"$set": doc}, upsert=True)
            for doc in map(self.apply_schema, docs)
        ]

        try:
            result = self.collection.bulk_write(ops, ordered=False)
            self.logger.info(
                "update TNS mirror",
                extra={
                    "inserted": result.upserted_count,
                    "modified": result.modified_count,
                },
            )
        except BulkWriteError as bwe:
            pass

    def get_names_for_location(self, ra, dec, radius):
        ra = ra if ra < 180.0 else ra - 360.0
        cursor = self.collection.find(
            {
                "pos": {
                    "$nearSphere": [ra, dec],
                    "$maxDistance": numpy.radians(radius / 3600.0),
                }
            },
            {"_id": 1},
        )
        return [TNSName.from_index(doc["_id"]) for doc in cursor]

    @classmethod
    def apply_schema(cls, doc):
        # convert name to integer for more efficient indexing
        doc["_id"] = int(TNSName.from_str(doc["objname"]))
        doc["pos"] = cls.geojson_key(doc.pop("radeg"), doc.pop("decdeg"))
        for k in ["ra", "dec"]:
            del doc[k]
        with suppress(ValueError):
            doc["discoverydate"] = Time(doc["discoverydate"]).datetime
        return doc

    @staticmethod
    def geojson_key(ra, dec):
        # geoJSON needs longitude between -180 and +180
        ra = ra if ra < 180.0 else ra - 360.0
        return {"type": "Point", "coordinates": [ra, dec]}


if __name__ == "__main__":
    import asyncio
    import logging
    from argparse import ArgumentParser

    from ampel.contrib.hu.t3.tns.TNSClient import TNSClient

    logging.basicConfig(level="INFO")

    parser = ArgumentParser(description="Update local TNS mirror")
    parser.add_argument("mongo_uri")
    parser.add_argument("--db", default="TNS")
    parser.add_argument("--api-key")
    parser.add_argument("--since")
    parser.add_argument("--no-update", default=False, action="store_true")
    parser.add_argument("--timeout", type=float, default=60)
    parser.add_argument("--max-requests", type=int, default=8)
    args = parser.parse_args()

    client = TNSClient(
        args.api_key, args.timeout, args.max_requests, logging.getLogger()
    )
    db = TNSMirrorDB(args.mongo_uri, args.db)
    if args.no_update:
        existing = {doc["objname"] for doc in db.collection.find({}, {"objname": 1})}
    else:
        existing = set()

    async def run():
        chunk = []
        async for doc in client.search(exclude=existing, public_timestamp=args.since):
            chunk.append(doc)
            if len(chunk) >= 100:
                db.add_sources(chunk)
                del chunk[:]
        db.add_sources(chunk)
        del chunk[:]

    asyncio.get_event_loop().run_until_complete(run())
