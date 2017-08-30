from stoflag import StoFlag


class SpectrumFlag(StoFlag):
    @classmethod
    def BEST(StoFlag):
        flag=StoFlag(1,'Best')
        return flag
    @classmethod
    def GOOD(StoFlag):
        flag=StoFlag(2,'Good')
        return flag
    @classmethod
    def UGLY(StoFlag):
        flag=StoFlag(5,'UGLY')
        return flag
