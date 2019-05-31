# This script contains all the default parameters to run one experiment

import math
from point import Point
from solver import Solver
from compartment import Compartment
from substance import Substance
from clock import Clock

class Experiment():

  # Define a default soma

  somaParent = None
  somaCoord = Point((0,0,0))
  somaRadie = 10e-6

  neuriteEndPoint = Point((30e-6,0,0))
  neuriteRadie = 0.25e-6 # Axon, diameter around 0.5 micrometer

  # These constants are estimated from neurite growth.
  # Assuming no growth at 5 mM (Walker et al 1998)
  # Assuming 0.22mm/hour at 10 mM for DENDRITES
  #  neuriteGrowthPoly = 1.2e-5*1e-3 # in m/(s mM), SI unit is mM
  #  neuriteGrowthDepoly = 6.1e-11 # in m/s # UNIT ERROR?

  # Dendrite growth speed
  #  neuriteGrowthPoly =  0.22/(5e-3*3600*1e3)# in m/(s mM), SI unit is mM
  #  neuriteGrowthDepoly = 0.22/(3600*1e3) # in m/s

  # Axon growth speed 33e-6/3600 m/s, Myers, ..., Baas 2006
  # No growth at 5e-3 mM tubulin (Walker et al 1988)
  neuriteGrowthPoly =  33e-6/(5e-3*3600)# in m/(s mM), SI unit is mM
  neuriteGrowthDepoly = 33e-6/3600 # in m/s

  tubulinConcentrationNeurite = 5e-3
  tubulinConcentrationSoma = 10e-3

  # Odde 1997: Assuming 1640 subunits/micrometer
  # gives 2.724e-15 mol/m (for ONE microtubuli) -
  # - he seems to be using same assumption as below:
  # Alternatively:
  # Assuming 8nm per tubulin, 13 ring --> 2.7e-15 mol/m
  # (Kinoshita et al 2001, gives 8nm per tubulin)
  # BUT, how many microtubulin in a growth cone?!!
  # 15 microtubuli in Mamalian axon Mudrakola, Zhang, Cui 2009

  tubulinQuantityPerLength = 2.7e-15*15

  # tubulinDiffusionConstant = 8.3e-12 # in m^2/s = 30000 mum^2/h
  # tubulinDiffusionConstant = 4.54e-13 # in m^2/s # 4.54e-9 cm^2/s Okabe,Hirokava 1990
  # tubulinDiffusionConstant = 4.3e-11 # in m^2/s Kruoglova, Vercamnen, Engelborghs 2004 (assuming 1mM Mg2+, Hille says 0.1-2mM free Mg2+ in neurons)
  tubulinDiffusionConstant = 8.591e-12 # in m^2/s Galbraith,..., Gallant 1999), squid giant axon

  tubulinDegradationConstant = 5.67e-7 # s^-1 (Miller & Samuels 1997)
          # about 1 week, Forgue, Dahl 1978 ==> 1.6e-6

  # tubulinActiveTransportRate = 0 # 2.8e-8 m/s = 100e-6 m/h
  tubulinActiveTransportRate = 440e-9*6e-3

  # Max scale factor is 20!! Here we use 6e-3, so below upper limit.

  # This gives us a diffusion length of 2*sqrt(D*t) = /t=50s/ = 9.5e-6

  #
  # !!! How do we handle tubulin consumption when the neurite diameter growths
  #     How many microtubuli should we assume in growth cone
  #

  # Assume 10e-3 mM tubulin
  tubulinQuantitySoma = \
    10e-3*(4*math.pi*(somaRadie**3)/3)

  tubulinQuantityNeurite = \
    5e-3*(20e-6*math.pi*(neuriteRadie**2))
    #10e-3*(20e-6*math.pi*(neuriteRadie**2))

  tubulinQuantityGrowthCone = \
    5e-3*(10e-6*math.pi*(neuriteRadie**2))
    # 10e-3*(10e-6*math.pi*(neuriteRadie**2))

  # Produce enough to maintain steady state at 10mM (approximately)
  tubulinSomaProductionRate = tubulinQuantitySoma*tubulinDegradationConstant

  # There is autoregulation of beta-tubulin, see Laferriere, MacRae, Brown 1997

  # This function sets up the simulation with all the compartments
  # make a subclass to do another experiment
  def __init__(self, solver):

    # clockDt = 10.0
    # clockStart = 0.0
    # clockEnd = 100.0
    # solver.setClock(Clock(clockDt,clockStart,clockEnd))

    self.solver = solver

    self.soma = Compartment(solver, \
                            Experiment.somaParent, \
                            Experiment.somaCoord, \
                            Experiment.somaRadie)

    neuriteParent = self.soma

    self.neurite = Compartment(solver,neuriteParent, \
                               Experiment.neuriteEndPoint, \
                               Experiment.neuriteRadie)

    self.neurite.makeGrowthCone(Experiment.neuriteGrowthPoly, \
                                Experiment.neuriteGrowthDepoly, \
                                Experiment.tubulinQuantityPerLength)

    Substance("tubulin",self.soma,self.solver, \
              Experiment.tubulinConcentrationSoma*self.soma.volume(), \
              Experiment.tubulinDiffusionConstant, \
              Experiment.tubulinDegradationConstant, \
              Experiment.tubulinActiveTransportRate, \
              Experiment.tubulinSomaProductionRate)

    Substance("tubulin", self.neurite, self.solver, \
              Experiment.tubulinConcentrationNeurite*self.neurite.volume())

    self.solver.setSaveFile('output/output.txt')    

