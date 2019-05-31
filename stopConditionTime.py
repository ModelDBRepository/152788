# This stop condition activates after time T

from solver import Solver
from stopcondition import StopCondition

import pdb

class StopConditionTime(StopCondition):

  def __init__(self,solver,stopTime, stateSaveFile=None, simSaveFile=None, \
               growthCones=None,pertubation=None):

    self.solver = solver
    self.stopTime = stopTime
    self.stateSaveFile = stateSaveFile
    self.simSaveFile = simSaveFile

    self.growthCones = growthCones
    self.pertub = pertubation

  def stop(self):
    #pdb.set_trace()
    return (self.solver.clock.curTime() >= self.stopTime)

  def action(self):
    if(self.stateSaveFile):
      self.solver.saveState(self.stateSaveFile)

