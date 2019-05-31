# This experiment loads a SWC experiment, runs a growth simulation with
# all branches. Then it restarts the simulation and this time runs a second 
# simulation with one growth cone having different polymerisation rate.
#
# I am interested to see how the effect on the other growth cones speed 
# is with distance.
#

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
from stepActionSomaClamp import StepActionSomaClamp
from clock import Clock
from importSWC import ImportSWC

import pdb


class ExperimentDistanceDependence(Experiment):

  def __init__(self,solver,SWCfile,modifyGCnumber,saveFile):

    self.solver = solver

    # Additional parameters
    #Solver.nIterTransientRemoval = -1 # Estimate based on L = 2*sqrt(D*t)


    # Testing, try running without transient removal just to see how it looks
    Solver.nIterTransientRemoval = 0

    solver.solverType = 'impEuler'
    solver.preventSplittingOfGrowthCones = False

    self.clockDt = 0.5 # 2 # 10 # 0.1 #1.0 #10.0 # 50.0
    self.clockStart = 0.0
    self.clockEnd = 3600*10 # 3600*10 #3600*24*1 # 3600*24*5 # 505

    solver.setClock(Clock(self.clockDt,self.clockStart,self.clockEnd))

    SWCread = ImportSWC(SWCfile, solver, self)
    SWCread.makeTipsGrowthCones()

    # We want to keep track of net influx from parent
    saveFlux = True

    solver.setSaveFile(saveFile, saveFlux)    

    # We want to clamp the soma concentration
    if(solver.compartments[0].isSoma()):
      soma = solver.compartments[0]
      solver.addStepAction(StepActionSomaClamp(solver,soma,"tubulin", \
                                        Experiment.tubulinConcentrationSoma))
    else:
      print "Why god why is not the first compartment a soma?"
      pdb.set_trace()

    # Modify polymerisation rate of one growth cone
    if(modifyGCnumber >= 0):
      GCmod = solver.growthCones[modifyGCnumber]
      GCmod.tubulinPolymerisation = 2*GCmod.tubulinPolymerisation
      print("Changing polymerisation for growth cone " + modifyGCnumber.__str__() + "\n")
