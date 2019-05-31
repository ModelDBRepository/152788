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
from experimentLoaderPredictSpeed import ExperimentLoaderPredictSpeed

def main():

  # Add error checking
  if(len(sys.argv) < 3):
    print("Usage: " + sys.argv[0] + " summaryFile workerID\n")
    exit()

  summaryFile = sys.argv[1]
  workerID = sys.argv[2]

  print("Summary file: " + summaryFile + "\n" \
        + "Worker ID: " + workerID + "\n")

  maxElem = 5 #51 #51 #2010 #1010

  clockDt = 10.0 # This value gets overwritten...
  outputDt = 500.0
  maxCompLen = 10e-6 # 200e-6 #9e-6 # 2.5e-6
  minCompLen = 1e-6 # 50e-6 #1e-6 #0.5e-6

  fp = open(summaryFile,'r')

  for line in fp:

    words = line.split()
    jobID = words[0]
    infoFile = words[1]

    if(jobID == workerID):
      solver = Solver(maxElem,clockDt,outputDt,maxCompLen,minCompLen)

      experiment = ExperimentLoaderPredictSpeed(solver,infoFile)

      if(not experiment.aborted):
        solver.init(True) # Find steady state before starting
        solver.run()
        solver.finish()

        print "Data point done."

  fp.close()


if __name__ == "__main__":
    main()
