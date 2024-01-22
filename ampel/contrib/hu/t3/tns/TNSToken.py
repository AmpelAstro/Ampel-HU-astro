from ampel.base.AmpelBaseModel import AmpelBaseModel


class TNSToken(AmpelBaseModel):
    id: int
    name: str
    api_key: str
