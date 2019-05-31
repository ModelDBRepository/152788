from compartment import Compartment
from experiment import Experiment

import pdb
# Use pdb.set_trace() to enter debugger somewhere


class SlavedGrowthConeSpeed(Compartment):

  def __init__(self, \
               solver, \
               parent, \
               slaveFile=None, \
               endPoint=None, \
               radie=10e-6, \
               tubulinQuantityPerLength=Experiment.tubulinQuantityPerLength, \
               constraint="arclength"):

    self.solver = solver
    self.slaveFile = slaveFile   # This is the file with growth info
    self.parent = parent
    self.children = []
    self.substances = {}
    self.persistent = True
    self.constraint = constraint

    self.parentArcLength = None

    self.id = Compartment.nextID
    Compartment.nextID += 1

    
    # This is a slaved growth cone, growthSubstance is used when 
    # we have concentration driven growth.
    self.growthSubstance = None 

    if(slaveFile):
      importSlaveInfo(slaveFile)

    if(not endPoint):
      endPoint = Point((0,0,0))

    if(parent):
      if(self.solver.verbose):
        print "Created dend."

      self.parent.addChild(self)
      self.isSphere = False
    else:
      # Trying to create growth cone without soma first...
      print "Error: Make the soma first..."
      pdb.set_trace()

    self.endPoint = endPoint
    self.radie = radie
    self.isGrowthCone = True

    self.tubulinQuantityPerLength = tubulinQuantityPerLength

    # Add the compartment to the list of compartments
    solver.addCompartment(self)
    solver.addGrowthCone(self)

    self.curIdx = 0

  def setSlaveInfo(self,time,speed,arclength):
 
    self.slaveTime = []
    self.slavePoint = []
    self.curIdx = 0

    self.slaveTime = time
    self.slaveSpeed = speed
    self.slaveLength = arclength

  def importSlaveInfo(self,slaveFile):

    fp = open(slaveFile,'rb')

    self.slaveTime = []
    self.slaveSpeed = []
    self.slaveLength = []
    self.curIdx = 0

    for line in fp:
      words = line.split()
      self.slaveTime.append(words[0])
      self.slaveSpeed.append(words[1])
      self.slaveLength.append(words[2])

    print "Read ", len(slaveTime), " data points"


  def growthSpeed(self):

    growthSpeed = 0

    if(self.constraint ==  "speed"):
        # Identify the last time point in the vector that has not passed yet
        while(self.curIdx+1 < len(self.slaveTime) \
              and self.slaveTime[self.curIdx+1] < self.solver.clock.curTime()):
          self.curIdx = self.curIdx + 1

        if(self.solver.verbose):
          print "Growing with speed: " + str(self.slaveSpeed[self.curIdx]) \
            + " (" + str(self) + ")"

        growthSpeed = self.slaveSpeed[self.curIdx]
    elif(self.constraint == "arclength"):
      # Identify the next point in time vector to happen
      while(self.curIdx+1 < len(self.slaveTime) \
            and self.slaveTime[self.curIdx] < self.solver.clock.curTime()):
        self.curIdx = self.curIdx + 1

# Something is wrong with the slaving... !!!!!!!!!!!
      if(self.slaveTime[self.curIdx] > self.solver.clock.curTime()):
        dt = self.slaveTime[self.curIdx] - self.solver.clock.curTime()

        # How fast do we need to grow to reach the target length
        growthSpeed = (self.slaveLength[self.curIdx] - self.arcLength())/dt
      else:
        # We are out of data points
        growthSpeed = 0

      # !!!! Can we speed up the arcLength call somehow

    else:
      print "Unknown growth constraint."
      quit()
   
    return growthSpeed


  def grow(self):
    # Overloading the compartment grow function

    distance = self.growthSpeed()*self.solver.clock.dt

    # Prevent shrinking to smaller than minimal compartment length
    #if(distance + self.length() < self.solver.minCompartmentLength):
    #  distance = self.length() - self.solver.minCompartmentLength

    # Prevent growing more than the amount of available tubulin allows
    if(distance > self.substances["tubulin"].quantity \
                  / self.tubulinQuantityPerLength):
      distance = 0.99*self.substances["tubulin"].quantity \
                   / self.tubulinQuantityPerLength

    # We need to consume the tubulin associated with the growth
    self.substances["tubulin"].quantity -= \
      distance * self.tubulinQuantityPerLength

    growthVector = self.direction()*distance
    self.growDir(growthVector)
    self.solver.markModified(self)

    if(distance > 0):
      if(self.parent.isNeurite() \
        and not self.parent.isBranchPoint() \
        and self.length() > self.solver.minCompartmentLength):

        self.parent.growDir(growthVector)
        self.solver.markModified(self.parent)
    else:
      if(self.parent.isNeurite() \
         and not self.parent.isBranchPoint()):

        growthVectorParent = self.parent.direction()*distance
        self.parent.growDir(growthVectorParent)
        self.solver.markModified(self.parent)

    

  def init(self):
    self.curIdx = 0

    



