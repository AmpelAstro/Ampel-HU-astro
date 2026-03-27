from collections.abc import Iterator

from pyarrow import parquet as pq

from ampel.abstract.AbsAlertLoader import AbsAlertLoader


class ParquetAlertLoader(AbsAlertLoader[dict]):
    """
    Alert loader for PyArrow tables. This is used for testing and development purposes.
    """

    path: str

    def __init__(self, **kwargs) -> None:
        super().__init__(**kwargs)
        self._batches = pq.ParquetFile(self.path).iter_batches(1)
        self._batch: None | Iterator[dict] = None

    def _next_batch(self) -> Iterator[dict]:
        # may raise StopIteration if no more batches are available
        return iter(next(self._batches).to_pylist())

    def __next__(self):
        if self._batch is None:
            self._batch = self._next_batch()
        try:
            return next(self._batch)
        except StopIteration:
            self._batch = self._next_batch()
            return next(self._batch)
