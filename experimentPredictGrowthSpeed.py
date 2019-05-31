import sys
import getopt

import os.path

import math
import numpy
import clock
import substance
import compartment
import saveData

from solver import Solver
# from solverDense import Solver
from point import Point

from experiment import Experiment
from substance import Substance
from compartment import Compartment
from stepActionSpeed import StepActionSpeed
from clock import Clock

class ExperimentPredictGrowthSpeed(Experiment):

  # Should this code really use a morph file...

  def __init__(self,solver,morphFile,growthConeFiles, \
               predictGCnumber, saveFile):

    self.solver = solver

    # Additional parameters
    Solver.nIterTransientRemoval = -1 # Estimate based on L = 2*sqrt(D*t)
    solver.solverType =  'impEuler'
    solver.preventSplittingOfGrowthCones = False

    self.clockDt = 10.0 # 10 # 0.1 #1.0 #10.0 # 50.0
    self.clockStart = 0.0
    self.clockEnd = 1.6e5 # This should be set to the last time entry 
                           # in growth cone file

    solver.setClock(Clock(self.clockDt,self.clockStart,self.clockEnd))

    self.growthConeInfo = []

    self.readMorphFile(morphFile)

    ctr = 0

    for gcFile in growthConeFiles:
      if(ctr in predictGCnumber):
        predictFlag = True
      else:
        predictFlag = False

      self.readGrowthConeFile(gcFile, predictFlag)
      ctr = ctr + 1

    # This function is called by solver every step
    StepActionSpeed(solver,self.growthConeInfo, predictGCnumber)

    # We want to keep track of net influx from parent
    saveFlux = True

    solver.setSaveFile(saveFile, saveFlux)    



  def readMorphFile(self,morphFile):

    fp = open(morphFile,'r')

    coords = []

    for line in fp:
      words = line.split()
      coords.append(Point((float(words[0]),float(words[1]),0)))  

    # Direction of initial part of neurite
    dirVec = coords[1]-coords[0]
    dirVec /= dirVec.norm()

    # Put soma before the neurite
    somaPos = coords[0] - dirVec*Experiment.somaRadie

    self.soma = Compartment(self.solver, \
                            Experiment.somaParent, \
                            somaPos, \
                            Experiment.somaRadie)

    Experiment.tubulinQuantitySoma = Experiment.tubulinConcentrationSoma \
                                     * self.soma.volume()
    Experiment.tubulinSomaProductionRate = Experiment.tubulinQuantitySoma \
                                         * Experiment.tubulinDegradationConstant

    Substance( "tubulin",self.soma,self.solver, \
                Experiment.tubulinConcentrationSoma*self.soma.volume(), \
                Experiment.tubulinDiffusionConstant, \
                Experiment.tubulinDegradationConstant, \
                Experiment.tubulinActiveTransportRate, \
                Experiment.tubulinSomaProductionRate )

    neuriteParent = self.soma

    ctr = 0

    for endPoint in coords:
      if(ctr == 0):
        # Skip first one
        ctr = ctr + 1
        continue

      neurite = Compartment(self.solver, \
                            neuriteParent, \
                            endPoint, \
                            Experiment.neuriteRadie)

      Substance("tubulin",neurite,self.solver, \
                Experiment.tubulinConcentrationNeurite*neurite.volume())

      # The next neurite attaches to this one
      neuriteParent = neurite

      ctr = ctr + 1

  def readGrowthConeFile(self,growthConeFile, predictFlag):

    speed = []
    time = []
    arclength = []
    startDist = None

    fp = open(growthConeFile,'r')

    for line in fp:
      words = line.split()
      time.append(float(words[0]))
      speed.append(float(words[1]))  
      arclength.append(float(words[2]))

      if(not startDist):
        startDist = float(words[2])

    self.growthConeInfo.append([time, startDist, speed, predictFlag, arclength])
