from solver import Solver

class StopCondition():

  def __init__(self,solver,growthCones,arcDistThresh,\
               stateSaveFile=None, simSaveFile=None,pertub=1.2):
    self.solver        = solver
    self.growthCones   = growthCones
    self.arcDistThresh = arcDistThresh
    self.stateSaveFile = stateSaveFile
    self.simSaveFile   = simSaveFile
    self.pertub = pertub

  def stop(self):
    stopFlag = False

    for gc in self.growthCones:
      if(gc.arcDist() > self.arcDistThresh):
        if(self.solver.verbose):
          print "Stop condition true : ", self.arcDistThresh
        stopFlag = True

    return stopFlag

  def action(self):
    if(self.stateSaveFile):
      self.solver.saveState(self.stateSaveFile)

  def load(self):

    if(self.solver.verbose):
      print "Reloading old simulation ", self.stateSaveFile

    self.solver = self.solver.loadState(self.stateSaveFile)

    if(self.simSaveFile):

      if(self.solver.verbose):
        print "Setting save file ", self.simSaveFile

      if(self.solver.output):
        print("Flux status: " + str(self.solver.output.saveFlux))
        self.solver.setSaveFile(self.simSaveFile,self.solver.output.saveFlux)
      else:
        self.solver.setSaveFile(self.simSaveFile)

      Solver.nIterTransientRemoval = 0

      # Do not reset the solver clock
      # self.solver.clock.reset()

      self.solver.init() # Need to open the new save file
      return self.solver

  def pertubation(self):
    gc = self.solver.growthCones[0]

    # Increase polymerisation rate by 20%  
    gc.tubulinPolymerisation *= self.pertub




