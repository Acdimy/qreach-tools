import pyqreach

# Initialize the transition system
ts = pyqreach.TransitionSystem()
ts.setInitLocation(0)

# Add locations to the transition system
# Here we create 5 locations with indices from 0 to 4
loc_list = []
for i in range(5):
    loc = pyqreach.Location(8,i)
    ts.addLocation(loc)
    loc_list.append(loc)

# Define the operations
op = pyqreach.QOperation(["00000000"])
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
ts.addRelation(2, 4, opi)
ts.addRelation(3, 4, opi)
ts.addRelation(4, 4, opx)

ts.setAnnotation([[4, op]])

ts.computingFixedPointPre()
ts.printDims(4)
ts.printSupp(4)
