#
# This script allows the user to read in a SWC morphology
#

import os
import math

from point import Point
from solver import Solver
import numpy
from substance import Substance
from compartment import Compartment
from experiment import Experiment

import pdb       # pdb.set_trace()



class ImportSWC():

  def __init__(self, SWCfile, solver, experiment,permitMerge=True):

    print "Loading " + SWCfile

    self.fileName = SWCfile
    self.solver = solver
    self.soma = None

    fp = open(self.fileName,"r")

    str = fp.readline()
    IDlookup = dict()

    somaPerifery = []

    while(str):

      if(str[0] == '#'):
        # Skip comment lines
        str = fp.readline()
        continue

      tokens = str.split()

      id = int(tokens[0])
      type = int(tokens[1])
      x = float(tokens[2])*1e-6
      y = float(tokens[3])*1e-6
      z = float(tokens[4])*1e-6
      r = float(tokens[5])*1e-6
      parentid = int(tokens[6])

      # pdb.set_trace()

      if(self.solver.verbose):
        print("id: "+id.__str__()+" type: " + type.__str__() + " (" + x.__str__() + "," + y.__str__() + "," + z.__str__() + ")" \
            + " r: " + r.__str__() + " pid: " + parentid.__str__() + "\n")

      if(type == 1):
        # Soma

        somaPerifery.append(Point((x,y,z)))

        if(not self.soma):
          # Soma not created yet, add it

          self.soma = Compartment(solver, None, \
                                  Point((x, y, z)), r)
          IDlookup[id] = self.soma

          Experiment.tubulinQuantitySoma = \
              Experiment.tubulinConcentrationSoma \
              * self.soma.volume()
          Experiment.tubulinSomaProductionRate = \
              Experiment.tubulinQuantitySoma \
              * Experiment.tubulinDegradationConstant

          Substance( "tubulin",self.soma,self.solver, \
                     Experiment.tubulinConcentrationSoma*self.soma.volume(), \
                     Experiment.tubulinDiffusionConstant, \
                     Experiment.tubulinDegradationConstant, \
                     Experiment.tubulinActiveTransportRate, \
                     Experiment.tubulinSomaProductionRate )

        else:
          # The soma already exists, only increase its volume

          # Point to the old soma
          IDlookup[id] = self.soma

          # Use the somaPerifery points to calculate the new soma
          centerPoint = Point((0,0,0))        

          for p in somaPerifery:
            centerPoint = centerPoint + p / len(somaPerifery)

          maxR = 0

          for p in somaPerifery:
            pDiff = centerPoint - p
            maxR = max(maxR,pDiff.norm())

         
          self.soma.endPoint = centerPoint
          self.soma.radie = maxR

          Experiment.tubulinQuantitySoma = \
              Experiment.tubulinConcentrationSoma \
              * self.soma.volume()

          Experiment.tubulinSomaProductionRate = \
              Experiment.tubulinQuantitySoma \
              * Experiment.tubulinDegradationConstant

          tubId = self.soma.substances["tubulin"].id 

          solver.quantity[tubId,0] = Experiment.tubulinQuantitySoma

          self.soma.substances["tubulin"].productionRate \
            = Experiment.tubulinSomaProductionRate


      else:

        # Not soma, so it is a neurite... 

        parent = IDlookup[parentid]
     
        if(parent == self.soma):
          # We need to verify that all compartments that are attached to the 
          # soma do not end within the soma.

          p = Point((x,y,z)) - self.soma.endPoint
          if(p.norm() - self.soma.radie <= 0):
            print "Found a compartment completely within soma radius, removed."
            # pdb.set_trace()

            IDlookup[id] = self.soma
            str = fp.readline()

            continue

        if(parent.isNeurite() and (not parent.isBranchPoint()) \
           and parent.length() < solver.minCompartmentLength \
           and permitMerge):
          print "importSWC: Merged compartment with parent"

          # The parent is too short, add this compartment to it
          oldLength = parent.length()
          oldRadie = parent.radie
          oldVolume = parent.volume()
          parent.endPoint = Point((x, y, z))

          # Scale radie if needed
          newLength = parent.length()

          if(newLength == 0):
            print "Merged compartment has length 0"
          else:
            parent.radie = (oldRadie*oldLength+r*(newLength-oldLength)) \
                            / newLength

          newVolume = parent.volume()

          # Scale quantity in all compartments
          for (name,subs) in parent.substances.iteritems():
            if(name == "tubulin"):
              # solver.quantity[subs.id,0] *= newVolume/oldVolume
              solver.quantity[subs.id,0] = \
                newVolume * Experiment.tubulinConcentrationNeurite
            else:
              print "Unknown substance " + name
              pdb.set_trace()

          # Bypass the new compartment which was merged
          IDlookup[id] = parent

        else:
          neurite = Compartment(solver, parent, \
                                Point((x, y, z)), r)    

          IDlookup[id] = neurite

          Substance("tubulin",neurite,self.solver, \
                    Experiment.tubulinConcentrationNeurite*neurite.volume())


      str = fp.readline()

    # Inspect to make sure all is ok - OK!
 
    # pdb.set_trace()


  # This function adds growth cones to the tip of the neurites

  def makeTipsGrowthCones(self):

   for comp in self.solver.compartments:
     if(len(comp.children) == 0 and not comp.isSoma()):
       comp.makeGrowthCone(Experiment.neuriteGrowthPoly, \
                           Experiment.neuriteGrowthDepoly, \
                           Experiment.tubulinQuantityPerLength)
