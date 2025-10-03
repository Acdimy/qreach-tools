import pyqreach

qnum = 8

# Initialize the transition system
ts = pyqreach.TransitionSystem()
ts.setInitLocation(0)

# Add locations to the transition system
# Here we create 5 locations with indices from 0 to 4
loc_list = []
for i in range(5):
    loc = pyqreach.Location(qnum,i)
    ts.addLocation(loc)
    loc_list.append(loc)

# Define the operations
op0 = pyqreach.QOperation(["00000000"])
op1 = pyqreach.QOperation(["10000000"])
oph = pyqreach.QOperation("H", qnum, [0], [])
opx = pyqreach.QOperation("X", qnum, [0], [])
opy = pyqreach.QOperation("Y", qnum, [0], [])
opz = pyqreach.QOperation("Z", qnum, [0], [])
opm0 = pyqreach.QOperation("meas0", qnum, [0], [])
opm1 = pyqreach.QOperation("meas1", qnum, [0], [])
opi = pyqreach.QOperation("I", qnum, [0], [])

# Add the operations to the transition system
# The topological order of the transition system: 
# 0 -> 1, 1 -> 2, 1 -> 3, 2 -> 4, 3 -> 4, 4 -> 4
ts.addRelation(0, 1, oph)
ts.addRelation(1, 2, opm0)
ts.addRelation(1, 3, opm1)
ts.addRelation(2, 4, opi)
ts.addRelation(3, 4, opi)
ts.addRelation(4, 4, opi)

ts.setAnnotation([[4, op0]])

ts.computingFixedPointPre()
ts.printDims(0)
ts.printSupp(0)
ts.printDims(2)
ts.printDims(3)

print("Transition system location 0 satisfies operation op0:", 
      ts.satisfy(0, op0))
