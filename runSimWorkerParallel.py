import sys
import getopt

import multiprocessing
import subprocess

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
from experimentLoader import ExperimentLoader

def main():

  # Add error checking
  if(len(sys.argv) < 2):
    print("Usage: " + sys.argv[0] + " summaryFile\n")
    exit()

  summaryFile = sys.argv[1]

  print("Summary file: " + summaryFile + "\n")

  fp = open(summaryFile,'r')

  jobList = []
  results = []

  for line in fp:

    words = line.split()
    jobID = words[0]
    infoFile = words[1]

    jobList.append(infoFile)

  fp.close()

  # print jobList

  pool = multiprocessing.Pool(None)
  r = pool.map_async(runSim, jobList, callback=results.append)
  r.wait()
  print results

  print "All jobs done."


def runSim(infoFile):

  print "Processing ", infoFile

  maxElem = 5 #51 #51 #2010 #1010

  clockDt = 10.0 # 2000.0 #50.0
  outputDt = 500.0
  maxCompLen = 10e-6 # 200e-6 #9e-6 # 2.5e-6
  minCompLen = 1e-6 # 50e-6 #1e-6 #0.5e-6

  solver = Solver(maxElem,clockDt,outputDt,maxCompLen,minCompLen)
  experiment = ExperimentLoader(solver,infoFile)

  if(not experiment.aborted):

    solver.init()

    # Save a set of different starting conditions
    for stopCond in experiment.stopInfo:
      solver.addStopCondition(stopCond.stop)
      solver.run()
      stopCond.action()

    solver.run()
    solver.finish()

    # Reload old starting conditions, but this time change the
    # polymerisation rate for one of the growth cones
    for stopCond in experiment.stopInfo:
      # Reload the old solver
      solver = stopCond.load()
      solver.clearStopConditions()
      stopCond.pertubation()  
      solver.run()
      solver.finish()

      print "Data point done: " + infoFile

    return True

  else:

    return False




if __name__ == "__main__":
    main()
