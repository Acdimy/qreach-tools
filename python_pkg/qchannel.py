class Qchannel:
    def __init__(self, name="", pos=[-1,-1], params=[]) -> None:
        self.name = name
        self.pos = pos
        self.params = params

class Noise(Qchannel):
    def __init__(self, name="", pos=[-1, -1], params=[], prob=1.) -> None:
        # legal types: ad, Bflip, Pflip, measure
        super().__init__(name, pos, params)
        self.prob = prob
    def check_legalty(self) -> int:
        if self.name == "ad":
            pass
        elif self.name == "Bflip":
            pass
        elif self.name == "Pflip":
            pass
        else:
            pass

class Measurement(Qchannel):
    def __init__(self, name="", pos=[-1, -1], params=[], reset=True) -> None:
        super().__init__(name, pos, params)
        self.params = [reset]

class SequentialCircuit:
    def __init__(self) -> None:
        pass
