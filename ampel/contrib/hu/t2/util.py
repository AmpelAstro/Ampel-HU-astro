from collections.abc import Mapping
from typing import Any, overload

from ampel.types import TBson
from ampel.view.T2DocView import T2DocView


@overload
def get_payload(view: T2DocView, *, code: int | None = None) -> Mapping[str, Any]: ...


@overload
def get_payload(
    view: T2DocView, type: type[TBson], *, code: int | None = None
) -> TBson: ...


def get_payload(
    view: T2DocView,
    type: type[TBson] = Mapping[str, Any],  # type: ignore[assignment]
    *,
    code: int | None = None,
) -> TBson:
    """Get payload as given type, raising ValueError if not found"""
    return view.get_payload(type, code=code, raise_exc=True)
