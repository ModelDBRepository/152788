import math
import numpy
import clock
import substance
import compartment
import saveData
from solver import Solver
from point import Point

from experiment import Experiment
from multiplexperiments import MultipleExperiments

import pickle

pickleFile = 'output/two-GC-transient-inital-state-1.pickle.txt'
outputFile = 'output/two-GC-output-1.txt'

solver = Solver()
solver = solver.loadState(pickleFile)

solver.setSaveFile(outputFile)
solver.nIterTransientRemoval = 0
solver.init()

gc = solver.growthCones[0]
gc.tubulinPolymerisation *= 1.2

solver.clearStopConditions()
solver.run()
