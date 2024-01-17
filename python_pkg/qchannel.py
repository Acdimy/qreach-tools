class Qchannel:
    def __init__(self, name="", pos=[-1,-1], params=[]) -> None:
        self.name = name
        self.pos = pos
        self.params = params

class Noise(Qchannel):
    def __init__(self, name="", pos=[-1, -1], params=[]) -> None:
        super().__init__(name, pos, params)

class Measurement(Qchannel):
    def __init__(self, name="", pos=[-1, -1], params=[]) -> None:
        super().__init__(name, pos, params)
