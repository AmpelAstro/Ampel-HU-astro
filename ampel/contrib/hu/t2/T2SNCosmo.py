from ampel.abstract.AbsLightCurveT2Unit import AbsLightCurveT2Unit
from ampel.view.LightCurve import LightCurve


class T2SNCosmo(AbsLightCurveT2Unit):
    def process(self, arg: LightCurve) -> None:
        raise NotImplementedError()
