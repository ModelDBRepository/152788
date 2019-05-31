# This stop condition is always true.


from solver import Solver
from stopcondition import StopCondition

class StopConditionAlwaysTrue(StopCondition):


  def __init__(self,solver,growthCones,\
               stateSaveFile=None, simSaveFile=None,pertub=1.2):
    self.solver        = solver
    self.growthCones   = growthCones
    self.stateSaveFile = stateSaveFile
    self.simSaveFile   = simSaveFile
    self.pertub = pertub

  def stop(self):
    return True

