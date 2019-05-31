import math
import numpy

from point import Point
from substance import Substance

# Debugger
import pdb



#import simulate

#
# The compartment keeps track of the following
# * Location of the cylinder end point (or center for spheres)
# * The parent
# * A list of the children
# * Radie of the compartment
#
# The following values are derived from the available data
#
# * Length
# * Volume
# * Cross sectional area
# * Distance between compartment centres
#
#
# !!! The diameter should increase at branches
# d_p^e = d_c1^e + d_c2^e, e=2 cross area matching
# (Graham, van Ooyen)
#
#
# Zhao, ..., Murai 1994 (peripheral nerves!!)
# Number of MT per (unmyelinated) axon: 21 +/- 2
# 65 +/- 6 MT per micrometer^2 = 65e12 m^-2
#
# Mammalian axon 500-1000 nm in diameter have about 15 microtubule
# Mudrakola,Zhang,Cui 2009 (Cell)
#
# Diameter of unmyelinated axon: 0.57 +/- 0.03
#
#


import point

class Compartment():
  '''
  Implements the Compartment class
  '''

  nextID = 0
  
  ##############################################################################

  def __init__(self, solver, parent=None, endPoint=None, radie=10e-6):
    self.solver = solver
    self.parent = parent
    self.children = []
    self.substances = {}

    self.parentArcLength = None

    self.id = Compartment.nextID
    Compartment.nextID += 1
    
    if(not endPoint):
      endPoint = Point((0,0,0))

    if(parent):
      if(self.solver.verbose):
        print "Created dend."

      self.parent.addChild(self)
      self.isSphere = False
    else:
      if(self.solver.verbose):
        print "Created soma"

      # We have a soma
      self.isSphere = True

    self.endPoint = endPoint
    self.radie = radie
    self.isGrowthCone = False

    # Add the compartment to the list of compartments
    solver.addCompartment(self)

    # Only for growth cones
    # self.persistent = False

  ##############################################################################

  def addSubstance(self,substance):
    if(substance.name in self.substances):
      print("Substance: " + substance.name + " already exist in compartment")

    self.substances[substance.name] = substance
    substance.compartment = self

  ##############################################################################

  def getSubstance(self,name):
    if name in self.substances:
      return self.substances[name]
    else:
      return None

  ##############################################################################

  def addChild(self,comp):
    if(not (comp in self.children)):
      self.children.append(comp)

  ##############################################################################

  def makeGrowthCone(self,tubulinPolymerisation, \
                          tubulinDepolymerisation, \
                          tubulinQuantityPerLength, \
                          persistent = False):

    if(self.children or not self.parent):
      print('Error: This compartment should not be a growth cone!')
    else:
      if(self.solver.verbose):
        print 'Making compartment a growth cone'

    self.tubulinPolymerisation = tubulinPolymerisation
    self.tubulinDepolymerisation = tubulinDepolymerisation
    self.tubulinQuantityPerLength = tubulinQuantityPerLength
    self.persistent = persistent # Means it will never be removed

    self.isGrowthCone = True

    # Tell the solver that this is a growth cone
    self.solver.addGrowthCone(self)

    # We need to request a quantity element to know how much tubulin
    # we can use each time step for growth

    # Setting compartment to None to mark that this is a virtual substance

    self.growthSubstance = Substance("tubulin",None,self.solver,0.0)
    
      
  ##############################################################################

## Derived properties, calculated from more basal variables
## 

  def center(self):
    if(self.isSphere):
      return self.endPoint
    else:
      return (self.endPoint + self.parent.endPoint)/2

  ##############################################################################

  def length(self):
    # Calculate length based on end point
    if(self.isSphere):
      # It is a sphere, no length
      return None

    elif(self.parent.isSphere):
      # Parent is a sphere, subtract parent radie
      p = self.endPoint-self.parent.endPoint
      return p.norm() - self.parent.radie

    else:
      # Parent and child are cylinders, just subtract end points
      p = self.endPoint-self.parent.endPoint
      return p.norm()

  ##############################################################################

  def arcLength(self):

    if(self.isSphere):
      return 0.0

    self.parentArcLength = self.parent.arcLength()

    return self.length() + self.parentArcLength

  ##############################################################################

  # Distance to soma
  def arcLengthTEST(self,recalculate=True):

    if(self.isSphere):
      return 0.0

    if(recalculate):
      self.parentArcLength = self.parent.arcLength(recalculate)

    # Use the cached version of parent arc length to avoid excessive
    # recursion. Normally in each time step when saving data to file
    # the parents' values will be updated before the children

    if(not self.parentArcLength):
      # We did not have a value for the parent, fall back to recursion
      self.parentArcLength = self.parent.arcLength(recalculate)

    return self.length() + self.parentArcLength

  ##############################################################################

  # This only calculates mantle surface
  def area(self):
    if(self.isSphere):
      return 4*math.pi*(self.radie**2)
    else:
      return 2*math.pi*self.radie*self.length()

  ##############################################################################

  def volume(self):
    if(self.isSphere):
      return 4*math.pi*(self.radie**3)/3
    else:
      return math.pi*(self.radie**2)*self.length()

  ##############################################################################

  def crossSectionAreaWithParent(self):
    ownCrossSection = math.pi*(self.radie**2)
    parentCrossSection = math.pi*(self.parent.radie**2)

    return min(ownCrossSection,parentCrossSection)

  ##############################################################################

  # Distance from center to parents center, used for diffusion
  def centerDistanceToParent(self):
    p = self.center() - self.parent.center()
    return p.norm()

  ##############################################################################

#  def isGrowthCone(self):
#    # If it does have a parent (then it is a dendrite or axon)
#    # If it does not have children then it is the end point
#    return ((not self.children) and self.parent)

  ##############################################################################

  def isNeurite(self):
    return self.parent

  ##############################################################################

  def isBranchPoint(self):
    # Does the compartment have more than one child?
    return (len(self.children) > 1)

  ##############################################################################

  def isSoma(self):
    return not self.parent

  ##############################################################################

## Helper functions, only used internally
## 

  # Which direction is neurite pointing in.
  def direction(self):
    if(self.parent):
      d = self.endPoint - self.parent.endPoint
      d /= d.norm()
      return d
    else:
      return None

  ##############################################################################

  # Helper function for grow
  def growDir(self, growthVector):
    #print "Growth vector: ", growthVector, " distance ", growthVector.norm()
    #print "old end point", self.endPoint, " parent: ", self.parent.endPoint
    self.endPoint += growthVector
    # for comp in self.solver.growthCones:
    #   print "GC: ", comp.endPoint

  ##############################################################################

  # This function grows the neurite by extending the second compartment
  # (this avoids too large fluctuations in concentration due to small
  # volumes. If there is no parent compartment which is a dendrite
  # then it grows itself.

  def grow(self,distance=None):
    if(self.isGrowthCone):
      if(distance):
       # We are being forced to grow a certain distance

        if(self.solver.noNegativeConcentrations \
           and distance*self.tubulinQuantityPerLength \
            > self.substances["tubulin"].quantity):   

          print "Growth capped by tubulin: " + str(distance)

          distance = min(distance, (self.substances["tubulin"].quantity*0.99) \
                                    / self.tubulinQuantityPerLength)

          print "Reduced forced growth: " + str(distance)

      else:
        # Default, use the quantity in growthSubstance to grow, if negative
        # we will shrink (it was allocated by the solver for this timestep)

        distance = self.growthSubstance.quantity \
                   / self.tubulinQuantityPerLength 

      if(distance > 0):
        # Growing
        growthVector = self.direction()*distance
        self.growDir(growthVector)
        self.solver.markModified(self)

        # Only move the 2nd compartment also if:
        # * 2nd compartment is neurite, but not branch point
        # * growth cone is not too small

        if(self.parent.isNeurite() \
          and not self.parent.isBranchPoint() \
          and self.length() > self.solver.minCompartmentLength):

          self.parent.growDir(growthVector)
          self.solver.markModified(self.parent)

          # Set the quantity to 0, spent it on growing
          self.growthSubstance.quantity = 0

      elif(distance < 0):
        # Shrinking
        if(self.persistent):
          # We are not allowed to remove this growth cone

          if(math.isnan(distance) or math.isnan(self.length())):
            print "NaN detected!"
            pdb.set_trace()

          if(distance+self.length() < self.solver.minCompartmentLength \
             and (self.parent.isSoma() or self.parent.isBranchPoint())):
            requestedDistance = distance

            # The max is to make sure the concentration does not go negative
            distance = min(0.999*self.solver.minCompartmentLength-self.length(),
                           requestedDistance \
                           +0.999*self.substances["tubulin"].quantity \
                            / self.tubulinQuantityPerLength)

            oldQuant = self.substances["tubulin"].quantity

            self.substances["tubulin"].quantity \
              += (requestedDistance-distance)*self.tubulinQuantityPerLength

            print "Requested shrinkage is too large, reduced to " \
                     + str(distance)
            if(self.substances["tubulin"].quantity < 0):
              print "Negative concentration, debug mode!!"
              pdb.set_trace()
            else:
              print "Quantity left: ", self.substances["tubulin"].quantity

        else:
          # Non-persistent compartment, should we remove it?

          if(-distance > self.length()):
            print "Removing growth cone!"
            self.parent.substances["tubulin"].quantity \
              += self.tubulinQuantityPerLength*self.length() \
                 + self.growthSubstance.quantity
            self.solver.markModified(self.parent)
            self.solver.removeGrowthCone(self)
            return

        growthVector = self.direction()*distance
        self.growDir(growthVector)
        self.solver.markModified(self)

        # Set the quantity to 0, accounted for the deficit by shrinking
        self.growthSubstance.quantity = 0

        if(self.parent.isNeurite() \
           and not self.parent.isBranchPoint()):

          growthVectorParent = self.parent.direction()*distance
          self.parent.growDir(growthVectorParent)
          self.solver.markModified(self.parent)

      else:
        # Zero growth requested, only modify if persistent flag is on
        # in case P=0,D=0 due to too low concentration earlier
        self.solver.markModified(self)

    else:
      print("Warning growth called for non-growth cone.")

  ##############################################################################

  def growOLD(self,distance):
    print "Entering OLD grow function"


    if(self.isGrowthCone):

      if(distance > 0):

        # Should we restrict growth to prevent negative concentrations
        if(self.solver.noNegativeConcentrations \
           and distance*self.tubulinQuantityPerLength \
            > self.substances["tubulin"].quantity):   

          print "Growth capped by tubulin: " + str(distance)

          distance = min(distance, (self.substances["tubulin"].quantity*0.99) \
                                    / self.tubulinQuantityPerLength)


          print "Reduced growth: " + str(distance)

      else:
        # Shrinking

        if(self.persistent):
          # If persistent compartment, do not allow it to shrink below
          # a certain size.

          if(distance+self.length() < self.solver.minCompartmentLength):
            distance = 0 #self.solver.minCompartmentLength - self.length()
            print "Requested shrinkage too large, reduced to " + str(distance)

        else:
          # Non-persistent compartment, should we remove it?

          if(-distance > self.length()):
            print "Removing growth cone!"
            self.parent.substances["tubulin"].quantity \
              += self.tubulinQuantityPerLength*self.length()

            self.solver.markModified(self.parent)
            self.solver.removeGrowthCone(self)
            return

      # Move growth cone accordingly
      gVec = self.direction()*distance
      self.growDir(gVec)      

      # Conversion factor for tubulin to microtubulin length
      # Assumes 8nm per tubulin, 13 ring, 1 microtubuli --> 2.7e-15 mol/m
      # Mudrakola, Zheng, Cui 2009: 15 MT in an axon at any given point
      # tubulinQuantityPerLength = 2.7e-15*15

      self.substances["tubulin"].quantity \
          -=  self.tubulinQuantityPerLength*distance


      # Should we move second compartment along also?

      # This call is only made by the growth cone to extend the 2nd compartment
      # * If the parent is a soma or if it is a branch point, then we do not
      #   extend 2nd compartment, instead we have to extend growth cone.
      # * If growth cone is small, then make it grow instead.

      # If parent is soma, then do not try and move it.
      # If parent is a branch point, then do not try and move it
      # If growth cone was supposed to grow, but is smaller than
      # min compartment length, then do not move parent, instead let
      # growth cone increase in size.

      if(self.parent.isNeurite() and not self.parent.isBranchPoint() \
           and not (distance > 0  \
                    and self.length() < self.solver.minCompartmentLength)):
        if(distance > 0):
          # If positive direction, follow growth cone
          self.parent.growDir(gVec)
        else:
          # if negative direction, retract along own axis
          gVecParent = self.parent.direction()*distance
          self.parent.growDir(gVecParent)
        self.solver.markModified(self.parent)
      else:
        self.solver.markModified(self)

    else:
      print("Warning growth called for non-growth cone.")

  ##############################################################################

  # This function decides the growth speed based on the tubulin concentration
  def growthSpeedOLD(self):
    if(not self.isGrowthCone):
      return None
    else:
      speed = self.tubulinPolymerisation*self.substances["tubulin"].conc \
                - self.tubulinDepolymerisation
      # Temporary change, using quantity for Arjen's simple model
      #speed = self.tubulinPolymerisation*self.substances["tubulin"].quantity \
      #          - self.tubulinDepolymerisation
      # print "Growing ", self, " at speed ", speed
      return speed

  ##############################################################################

  # Print a summary of the compartment to stdout
  def printComp(self):
    print "ID: ", self, "coords", self.endPoint, "parent:", self.parent
    print "Length: ", self.length(), "Radie: ", self.radie

    for subs in self.substances.itervalues():
      subs.printSubs()

  ##############################################################################

  def arcDist(self):
    if(self.isSoma()):
      return 0
    else:
      return self.parent.arcDist() + self.length()

  ##############################################################################

  def clearSolver(self):
    self.solver = None

  ##############################################################################

  def setSolver(self,solver):
    self.solver = solver

  ##############################################################################

  # Override this functino to allow a compartment to spawn a new growth cone
  def conditionalSpawnGrowthCone(self):
    # If this was true we should call removeCompartmentTransportReaction
    # to make sure we clear the corresponding parts of the matrix
    return False

  ##############################################################################
