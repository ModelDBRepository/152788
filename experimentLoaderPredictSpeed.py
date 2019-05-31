# This script reads in parameters from an external file, then does
# a prediction simulation of the growth cones.
#
#

import math
from point import Point
from solver import Solver
from compartment import Compartment
from substance import Substance
from experiment import Experiment
from clock import Clock
from experimentPredictGrowthSpeed import ExperimentPredictGrowthSpeed
from stepActionSpeed import StepActionSpeed
from stepActionSomaClamp import StepActionSomaClamp

import pdb       # pdb.set_trace()


import os.path

class ExperimentLoaderPredictSpeed(ExperimentPredictGrowthSpeed):

  
  def __init__(self,solver,loaderFile, overwriteFlag = False):

    # Additional parameters
    Solver.nIterTransientRemoval = -1 # Estimate based on L = 2*sqrt(D*t)
    solver.solverType = 'impEuler'

    self.clockDt = 10 
    self.clockStart = 0.0
    self.clockEnd = 1.6e5

    self.solver = solver

    self.useAlwaysTrueStopCond = False

    self.clampSomaConcentration = False

    Experiment.tubulinConcentrationSoma = 10e-3
    Experiment.tubulinConcentrationNeurite = 5e-3
    Experiment.tubulinQuantityGrowthCone = 5e-3

    solver.preventSplittingOfGrowthCones = False

    # Allow the file to override some parameters...
    self.loadInfo(loaderFile)

    solver.setClock(Clock(self.clockDt,self.clockStart,self.clockEnd))

 
    if(not overwriteFlag \
        and (os.path.isfile(self.saveFileName))):
      self.aborted = True
      print("One or more of output files already exist - aborting.")
      print(self.saveFileName + "\n")
      return(None)
    else:
      self.aborted = False

      # Lets create a dummy file to prevent others from writing to it
      fp = open(self.saveFileName,'w')
      fp.write("Job assigned to worker, this file will be overwritten shortly.")
      fp.close()


    # We want to keep track of net influx from parent
    saveFlux = True

    solver.setSaveFile(self.saveFileName, saveFlux)    

    # Set up the morphology
    self.readMorphFile(Experiment.morphFile)
 
    # Growth cone handling
    self.growthConeInfo = []
    ctr = 0
    for gcFile in Experiment.growthConeFiles:
      if(ctr in Experiment.predictGCnumber):
        print "GC " + str(ctr) + " is predicted."
        predictFlag = True
      else:
        print "GC " + str(ctr) + " is slaved."
        predictFlag = False

      # pdb.set_trace()

      self.readGrowthConeFile(gcFile, predictFlag)
      ctr = ctr + 1

    # This function is called by solver every step
    StepActionSpeed(solver,self.growthConeInfo, Experiment.predictGCnumber)

    if(self.clampSomaConcentration):
      # We want to clamp the soma concentration
      if(self.solver.compartments[0].isSoma()):
        print "Clamping soma concentration"
        soma = self.solver.compartments[0]
        self.solver.addStepAction(StepActionSomaClamp(self.solver,soma, \
                                  "tubulin", \
                                   Experiment.tubulinConcentrationSoma))
      else:
        print "Why god why is not the first compartment a soma?"
        pdb.set_trace()
    else:
      print "Not clamping soma concentration"




  def loadInfo(self, loaderFile):

    print "Loading ", loaderFile

    fp = open(loaderFile,'r')

    for line in fp:
      exec(line)

    fp.close()
