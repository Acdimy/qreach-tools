

class Tensor:
    def __init__(self, data=[], index=[]) -> None:
        self.data = data
        self.index_set = index
    def data(self):
        return getData(self.data)

def contract(ts1, ts2) -> Tensor:
    pass

def getData(bdd=None):
    return bdd