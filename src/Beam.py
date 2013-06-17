
class BeamError(Exception):
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)

class Beam:
    def __init__(self, fromalm=False, fromaipy=False):
        if fromaipy:
            self._fromaipy()
        elif fromalm:
            self._fromalm()
        else:
            raise(BeamError("Beam needs to come from aipy or from a healpix alm!"))

    def _fromaipy(self):
        pass

    def _fromalm(self):
        pass
