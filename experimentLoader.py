# This script reads in variations of experiments from a file
# so that they can be run.
#
# The file contains statements that modify the default settings
# for example:
#
# tubulinDiffusionConstant = 1e-12
# 
# see setupPredictSpeedSimulation.m

import math
from point import Point
from solver import Solver
from compartment import Compartment
from substance import Substance
from experiment import Experiment
from clock import Clock
from stepActionSomaClamp import StepActionSomaClamp


import os.path

from stopcondition import StopCondition
from stopConditionAlwaysTrue import StopConditionAlwaysTrue
from stopConditionTime import StopConditionTime

class ExperimentLoader(Experiment):

  def __init__(self,solver,loaderFile, overwriteFlag = False):

    # Additional parameters
    Solver.nIterTransientRemoval = None # Estimate based on L = 2*sqrt(D*t)
    Solver.removeInitStateTransients = False
    solver.solverType = 'impEuler'

    self.clockDt = 10 # 0.1 #1.0 #10.0 # 50.0
    self.clockStart = 0.0
    self.clockEnd = 2e5 #1e6 #1e5 #1e6 #5e6 #1000.0 #10000.0

    self.solver = solver

    self.useAlwaysTrueStopCond = False
    self.useStopCondTime = False
    self.stopCondTime = 3600
    self.clampSomaConcentration = False

    Experiment.tubulinConcentrationSoma = 10e-3
    Experiment.tubulinConcentrationNeurite = 5e-3
    Experiment.tubulinQuantityGrowthCone = 5e-3

    self.polyRateModifier = 1.2

    # Growth to allow for system to equalize
    self.preGrowth = 5e-6

    # Allow the file to override some parameters...
    self.loadInfo(loaderFile)

    solver.setClock(Clock(self.clockDt,self.clockStart,self.clockEnd))

    # Define at what points in the initialisation simulation we want to
    # save states for the different simulations
    self.stopInfo = []

    arcDist = self.distA + self.distB


    growthConeLength = 9e-6

    self.stateFile = self.saveFileName + ".pickle"
    self.saveFileOriginal = self.saveFileName + ".original.txt"
    self.saveFilePerturbed = self.saveFileName + ".perturbed.txt"
   
    if(not overwriteFlag \
        and (os.path.isfile(self.stateFile) \
             or os.path.isfile(self.saveFileOriginal) \
             or os.path.isfile(self.saveFilePerturbed))):
      self.aborted = True
      print("One or more of output files already exist - aborting.")
      print(self.stateFile + "\n" + self.saveFileOriginal + "\n" \
            + self.saveFilePerturbed)
      return(None)
    else:
      self.aborted = False

      # Lets create a dummy file to prevent others from writing to it
      fp = open(self.stateFile,'w')
      fp.write("Job assigned to worker, this file will be overwritten shortly.")
      fp.close()

    # Give the loader the opportunity to define another stop cond, and
    # then prevent this one from be created.
    if(self.useAlwaysTrueStopCond):
      self.stopCond = StopConditionAlwaysTrue(self.solver, \
                                              self.solver.growthCones, \
                                              self.stateFile, \
                                              self.saveFilePerturbed, \
                                              self.polyRateModifier)  
    elif(self.useStopCondTime):
      self.stopCond = StopConditionTime(self.solver, \
                                        self.stopCondTime, \
                                        self.stateFile, \
                                        self.saveFilePerturbed, \
                                        self.solver.growthCones, \
                                        self.polyRateModifier)
    else:
      self.stopCond = StopCondition(solver, solver.growthCones, arcDist, \
                                      self.stateFile, self.saveFilePerturbed, \
                                      self.polyRateModifier)

    self.stopInfo.append(self.stopCond)

    # We want to keep track of net influx from parent
    saveFlux = True

    solver.setSaveFile(self.saveFileOriginal, saveFlux)    


    # Set up the morphology
    self.soma = Compartment(solver, \
                            Experiment.somaParent, \
                            Experiment.somaCoord, \
                            Experiment.somaRadie)

    Experiment.tubulinQuantitySoma = Experiment.tubulinConcentrationSoma \
                                      * self.soma.volume()
    Experiment.tubulinSomaProductionRate = Experiment.tubulinQuantitySoma \
                                       *Experiment.tubulinDegradationConstant


    Substance( "tubulin",self.soma,self.solver, \
               Experiment.tubulinConcentrationSoma*self.soma.volume(), \
               Experiment.tubulinDiffusionConstant, \
               Experiment.tubulinDegradationConstant, \
               Experiment.tubulinActiveTransportRate, \
               Experiment.tubulinSomaProductionRate )

    neuriteEndA = Experiment.somaCoord \
                  + Point((self.distA + Experiment.somaRadie,0,0))

    if(self.preGrowth + growthConeLength > self.distB):
      print("ERROR: Too large preGrowth")
      self.preGrowth = (self.distB - growthConeLength) * 0.9;
    else:
      print("Using pre-growth: " + str(self.preGrowth))

    neuriteEndB1 = neuriteEndA + Point((1,1,0)) \
                    * (self.distB-self.preGrowth-growthConeLength)/math.sqrt(2.0)
    neuriteEndB2 = neuriteEndA + Point((1,-1,0)) \
                    * (self.distB-self.preGrowth-growthConeLength)/math.sqrt(2.0)

    growthConeEnd1 = neuriteEndB1 + Point((1,1,0))*growthConeLength/math.sqrt(2.0)
    growthConeEnd2 = neuriteEndB2 + Point((1,-1,0))*growthConeLength/math.sqrt(2.0)

    if(self.solver.verbose):
      print "neuriteEndA=", neuriteEndA
      print "neuriteEndB1=", neuriteEndB1
      print "neuriteEndB2=", neuriteEndB2

    neuriteA = Compartment(solver, \
                           self.soma, \
                           neuriteEndA, \
                           Experiment.neuriteRadie)

    neuriteB1 = Compartment(solver, \
                            neuriteA, \
                            neuriteEndB1, \
                            Experiment.neuriteRadie)

    neuriteB2 = Compartment(solver, \
                            neuriteA, \
                            neuriteEndB2, \
                            Experiment.neuriteRadie)

    growthCone1 = Compartment(solver, \
                              neuriteB1, \
                              growthConeEnd1, \
                              Experiment.neuriteRadie)

    growthCone2 = Compartment(solver, \
                              neuriteB2, \
                              growthConeEnd2, \
                              Experiment.neuriteRadie)


    Substance("tubulin", neuriteA, self.solver, \
              Experiment.tubulinConcentrationNeurite*neuriteA.volume())

    Substance("tubulin", neuriteB1, self.solver, \
              Experiment.tubulinConcentrationNeurite*neuriteB1.volume())

    Substance("tubulin", neuriteB2, self.solver, \
              Experiment.tubulinConcentrationNeurite*neuriteB2.volume())

    Substance("tubulin", growthCone1, self.solver, \
              Experiment.tubulinConcentrationNeurite*growthCone1.volume())

    Substance("tubulin", growthCone2, self.solver, \
              Experiment.tubulinConcentrationNeurite*growthCone2.volume())


    growthCone1.makeGrowthCone(Experiment.neuriteGrowthPoly, \
                               Experiment.neuriteGrowthDepoly, \
                               Experiment.tubulinQuantityPerLength)

    growthCone2.makeGrowthCone(Experiment.neuriteGrowthPoly, \
                               Experiment.neuriteGrowthDepoly, \
                               Experiment.tubulinQuantityPerLength)

    self.neuriteA  = neuriteA
    self.neuriteB1 = neuriteB1
    self.neuriteB2 = neuriteB2
    self.growthCone1 = growthCone1
    self.growthCone2 = growthCone2

    if(self.clampSomaConcentration):
      soma = self.solver.compartments[0]
      self.solver.addStepAction(StepActionSomaClamp(self.solver, \
                                soma, "tubulin", \
                                Experiment.tubulinConcentrationSoma))



  def loadInfo(self, loaderFile):

    fp = open(loaderFile,'r')

    for line in fp:
      exec(line)

    fp.close()
