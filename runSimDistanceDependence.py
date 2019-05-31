#
# #  python runSimDistanceDependence.py /home/hjorth/DATA/Martine/martineNeurolucida/SWC/cell4.swc 0
# -1 for reference

# #  python runSimDistanceDependence.py /home/hjorth/DATA/Martine/martineNeurolucida/SWC/cell1.swc 0
# -1 for reference (79 max?)


#  python runSimDistanceDependence.py /home/hjorth/DATA/Neuromorpho.org/allman/CNG-version-full/11o_pyramidal12aFI.CNG.swc 0
# -1 to do reference case (0-17 GC)

# python runSimDistanceDependence.py /home/hjorth/DATA/Neuromorpho.org/allman/CNG-version-full/17o_pyramidal13aFI.CNG.swc 0
# -1, 0 to 18

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
from experimentDistanceDependence import ExperimentDistanceDependence

def main():

  if(len(sys.argv) < 3):
    print("Usage: " + sys.argv[0] + " SWCfile growthConeNumber\n")
    exit()

  SWCfile = sys.argv[1]
  GCmodId = int(sys.argv[2])

  if(GCmodId >= 0):
    saveFile = SWCfile + "-GCmod-" + GCmodId.__str__() + "-out.txt"
  else:
    saveFile = SWCfile + "-NOGCmod-out.txt"

  print("SWC file: " + SWCfile + "\n" \
        + "GCmodID: " + GCmodId.__str__())

  maxElem = 100 # 1500 #100 # Initial size does not affect segfault...

  clockDt = 2.0 # This one gets overwritten by the experiment
  outputDt = 500.0
  maxCompLen = 11e-6 # 10e-6
  minCompLen = 5e-6 # 4e-6  # 1e-6
 
  solver = Solver(maxElem,clockDt,outputDt,maxCompLen,minCompLen, \
                  True, False) # Verbose,Do not use sparse

  solver.verbose = True
  experiment = ExperimentDistanceDependence(solver,SWCfile,GCmodId,saveFile)


  solver.init(False) # Do not remove transients just yet
  print "Setting guessed tubulin gradient"
  solver.setEstimatedGradient("tubulin", \
                              Experiment.tubulinConcentrationSoma, \
                              Experiment.tubulinConcentrationNeurite)
  solver.init(True) # Remove transients
  solver.run()
  solver.finish()


if __name__ == "__main__":
    main()



