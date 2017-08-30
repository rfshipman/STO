from stoflag import StoFlag


class ChannelFlag(StoFlag):
    @classmethod
    def GLITCH(StoFlag):
        flag=StoFlag(1,"Glitched Data")
        return flag
    @classmethod
    def DONTUSE(StoFlag):
        flag=StoFlag(2,"Do not use this channel")
        return flag
    @classmethod
    def LINE(StoFlag):
        flag=StoFlag(3,'Line detected')
        return flag
