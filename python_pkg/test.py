import pyqreach

# Initialize the transition system
ts = pyqreach.TransitionSystem()
ts.setInitLocation(0)

# Add locations to the transition system
# Here we create 5 locations with indices from 0 to 4
loc_list = []
for i in range(4):
    loc = pyqreach.Location(8,i)
    ts.addLocation(loc)
    loc_list.append(loc)

# Define the operations
op0 = pyqreach.QOperation(["00000000"])
op1 = pyqreach.QOperation(["10000000"])
oph = pyqreach.QOperation("H", 8, [0], [])
opx = pyqreach.QOperation("X", 8, [0], [])
opy = pyqreach.QOperation("Y", 8, [0], [])
opz = pyqreach.QOperation("Z", 8, [0], [])
opm0 = pyqreach.QOperation("meas0", 8, [0], [])
opm1 = pyqreach.QOperation("meas1", 8, [0], [])
opi = pyqreach.QOperation("I", 8, [0], [])

# Add the operations to the transition system
# The topological order of the transition system: 
# 0 -> 1, 1 -> 2, 1 -> 3, 2 -> 4, 3 -> 4, 4 -> 4
ts.addRelation(0, 1, oph)
ts.addRelation(1, 2, opm0)
ts.addRelation(1, 3, opm1)
ts.addRelation(2, 1, oph)
ts.addRelation(3, 3, opi)

ts.setAnnotation([[3, op1]])

ts.computingFixedPointPre()
ts.printDims(1)
ts.printDims(2)
ts.printDims(3)
ts.printSupp(0)
