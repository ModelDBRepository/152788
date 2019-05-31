#
# This script takes care of the solving of all the differential equations
#
# Compartments keep track of where they are in space, their size,
# parents and children.
#

#
# !!! Production rates are not properly handled during split and merges !!!
#     or growth for that matter. 
#

import numpy

#from scipy import sparse
#from scipy.sparse.linalg import spsolve
#from scipy.sparse import csr_matrix

import clock
import saveData
import pickle
import math

from compartment import Compartment
from substance import Substance

import sys
sys.setrecursionlimit(10000)


# Python debugger
import pdb
# Use pdb.set_trace() to enter debugger somewhere

class Solver():

  nIterTransientRemoval = -1 #1000

  # newCompartmentLength is no longer used!!

  def __init__(self,maxNumElement=500,clockDt=10.0, \
               outputDt=500.0, \
               maxCompartmentLength = 2.5e-6, \
               minCompartmentLength = 0.1e-6, \
               verbose=False, \
               useSparse=False, \
               removeInitStateTransients=False):

    print "Setting up solver"


    # This function is called every step, can be specified by user
    # to do other things.
    self._stepAction = []

    self.verbose = verbose
    self.useSparse = useSparse

    # The data structure is cyclic, pickle needs deeper recursion
    pickle.sys.setrecursionlimit(10000)

    # Size of solver matrix, affects speed
    self.maxNumElement = maxNumElement 

    self.solverType = 'impEuler' #'RK4' #'euler' #'RK2'

    # These are just default values for the clock, usually overwritten 
    # by the experiment class.
    clockDt = 10.0 #50.0
    clockStart = 0.0
    clockEnd = 1e5
    self.clock = clock.Clock(clockDt,clockStart,clockEnd)
    self.outputDt = outputDt
    
    self.maxCompartmentLength = maxCompartmentLength 
    self.minCompartmentLength = minCompartmentLength 

    nElem = self.maxNumElement

    self.nextId = 0


    self.quantity = numpy.zeros((nElem,1)) 
    self.production = numpy.zeros((nElem,1)) 

    print "useSparse=",useSparse,maxNumElement
    
    if(self.useSparse):
      self.transport = csr_matrix((nElem,nElem))
      self.reaction  = csr_matrix((nElem,nElem))
    else:
      self.transport = numpy.zeros((nElem,nElem))
      self.reaction  = numpy.zeros((nElem,nElem))

    self.removeInitStateTransients = removeInitStateTransients

    self.compartments = []
    self.growthCones = []

    self.output = None

    self.disableGrowth = False
    self.disableOutput = False 

    self.modifiedCompartments = []
    self.disabledPersistentGrowthCones = []

    self.stopCondition = []
    self.preventSplittingOfGrowthCones = True

    # Clear up old substances
    Substance._substanceList = {}

    self.noNegativeConcentrations = True
    self.allowGrowthConeMerge = False


  ##############################################################################

  def setClock(self,clock):
    self.clock = clock

  ##############################################################################

  def clearStepAction(self):
    self._stepAction = []

  ##############################################################################

  def addStepAction(self, stepAction):
    self._stepAction.append(stepAction)

  ##############################################################################

  def deleteStepAction(self, stepAction):
    self._stepAction.remove(stepAction)

  ##############################################################################

  def addCompartment(self,comp):
    if(not (comp in self.compartments)):
      self.compartments.append(comp)
    else:
      if(self.verbose):
        print('Compartment already present in solvers list')

  ##############################################################################

  def markModified(self, comp):
    if(not (comp in self.modifiedCompartments)):
      self.modifiedCompartments.append(comp)

  ##############################################################################

  def clearModified(self):
    self.modifiedCompartments = []

  ##############################################################################

  def addGrowthCone(self,comp):
    if(not (comp in self.compartments)):
      self.compartments.append(comp)

    if(not (comp in self.growthCones)):
      self.growthCones.append(comp)

  ##############################################################################

  def removeGrowthCone(self,comp):

    # We need to remove the reaction and transport before we change the
    # children or parents
    self.removeCompartmentTransportReaction(comp)

    self.growthCones.remove(comp)
    self.compartments.remove(comp)

    for (name,subs) in comp.substances.iteritems():

      parentId = comp.parent.substances[name].id 

      # Merge remaining quantity back into parent
      self.quantity[parentId,0] += self.quantity[subs.id,0]
      self.quantity[subs.id,0] = 0
      self.production[subs.id,0] = 0

    # Remove reference
    comp.parent.children.remove(comp)

    for child in comp.children:
      comp.parent.children.append(child)
      self.markModified(child)

    # Mark for recalculation
    self.markModified(comp.parent)

  ##############################################################################

  # This function should also keep track of freed substanceIDs so that
  # we can recycle them in case there is growth, retraction, growth, retr...
  # This will be useful later if we start modelling filopodia, spines etc
  # which growth the retract repeatedly.

  def getNewSubstanceId(self):
    if(self.nextId >= self.quantity.shape[0]):
      # raise Exception('Sparse matrix out of space','Increase maxNumComp') 
      newMaxNumElement = int(1.5*self.maxNumElement)
      self.increaseMaxNumElement(newMaxNumElement)

    newId = self.nextId
    self.nextId += 1

    return newId

  ##############################################################################

  def setSaveFile(self,name, saveFlux=False):
    #if(self.verbose):
    print "Writing output to ", name
    self.saveFlux = saveFlux

    if(saveFlux):
      print("Saving net influx from parent compartment also")

    self.output = saveData.SaveData(name,self,self.outputDt,saveFlux)

  ##############################################################################

  def init(self,removeTransients=None):

    if(removeTransients == None):
      # No value for removeTransient specified, use the simulation default
      removeTransients = self.removeInitStateTransients
    else:
      # User specified a new default value, use it and store it
      self.removeInitStateTransients = removeTransients

    if(self.verbose):
      print('Initialising the solver')

    # We temporarilly allow growth cones to be merged during the init phase
    self.allowGrowthConeMerge = True

    allDone = False
    updateCounter = 0

    oldGCsplitFlag = self.preventSplittingOfGrowthCones
    self.preventSplittingOfGrowthCones = False

    while(not allDone):
      updateCounter += 1

      print "Solver: Looping update spatial discretisation"
      allDone = True

      for comp in self.compartments:
        if((not comp.isSoma()) and  comp.length() < 0):
          print "Compartment ", comp, " has negative length!"
          print "Length: ", comp.length()
          quit()

        modifiedFlag = self.updateSpatialDiscretisation(comp)
        # print "Modified flag: ", modifiedFlag

        if(modifiedFlag):
          allDone = False

        if(updateCounter > 100):
          print "The spatial discretisation keep iterating, weird."
          print "This might be because the neurite is smaller than the min compartment length"
          pdb.set_trace()


    self.preventSplittingOfGrowthCones = oldGCsplitFlag

    print "Done updating discretisation"

    for comp in self.compartments:
      if(comp.volume() == 0):
        print "Missed a compartment with zero volume"
        self.updateSpatialDiscretisation(comp)

    print "All done... "

    for comp in self.compartments:
      if(comp.volume() == 0):
        print "Missed a compartment with zero volume - AGAIN?!!"
        self.updateSpatialDiscretisation(comp)


    for gc in self.growthCones:
      print "Growth cone length: ", gc.length()*1e6, " micrometers"


    self.allowGrowthConeMerge = False

    self.setProduction()
    self.calculateTransportReactionMatrix()

    # Only remove transients if we are at the start time
    #if(abs(self.clock.curTime() - self.clock.start) < self.clock.dt/2.0):
    if(removeTransients):
      print("Trying to remove transients!")
      #self.removeTransients()
      self.findInitStateAuto()
    else:
      print("Not removing transients!")

    if(self.output):
      self.output.init()

    for sa in self._stepAction:
      sa.init()

    # Save the simulation state to file...
    if(1):
      stateFile = self.output.fileName + ".PRE"
      self.saveState(stateFile)

    print "Max number of substances: ", self.maxNumElement, \
          " number of compartments: ", len(self.compartments)

  ##############################################################################

  def step(self):

    # print "Stepping time: " + str(self.clock.curTime())

    for sa in self._stepAction:
      # Call user supplied function
      sa.update()

    # print "Step action done."

    # Enter debug mode
    # pdb.set_trace()

    # Write the first time step before stepping forward
    if(abs(self.clock._t - self.clock.start) < self.clock.dt/2.0):
      if((not self.disableOutput) and self.output):
        self.output.write()

    # Disable growth if we want to find initial concentrations
    if(not self.disableGrowth):
      # Move growth cones
      for gc in self.growthCones:
        gc.grow()

      # New, allow spawning of new growth cones
      for comp in self.compartments:
        if(comp.conditionalSpawnGrowthCone()):
          # New growth cone spawned, mark compartment as modified
          self.modifiedCompartments.append(comp)

    for comp in self.compartments:
      self.updateSpatialDiscretisation(comp)

    # Make sure we update the growth cones that are temporarilly disabled
    # so that they can start growing again if concentration increases
    for comp in self.disabledPersistentGrowthCones:
      self.markModified(comp)

    for comp in self.modifiedCompartments:
      self.setTransportMatrixTwoSided(comp)
      self._setReactionProductionMatrix(comp)
      # comp.arcLength(False)

    # Mark that we have updated the modified compartments

    self.clearModified()

        
    if(False): #!!!!! Verification: Should be False for faster computation
      tTemp = self.transport
      rTemp = self.reaction  
      # Production matrix now also changed... for the growth cones
      # and their tubulin quantities 

      # Recalculate ENTIRE transport matrix every time
      self.setProduction()
      self.calculateTransportReactionMatrix()

      # Compare the locally modified to the complete recalculation
      # they should be the same.
      if(abs(self.transport-tTemp).sum()/abs(self.transport).sum() > 1e-15):

        print "self.transport"
        print self.transport[:,:].toarray()
        print "old self.transport"
        print tTemp[:,:].toarray()
        print "Difference: ", abs(self.transport-tTemp).sum()
        print (self.transport[:,:] - tTemp[:,:]).toarray()

        d = (self.transport[:,:] - tTemp[:,:])
        pdb.set_trace()


      if(abs(self.reaction-rTemp).sum() > 0):
        print "self.reaction"
        print self.reaction[:,:].toarray()
        print "old self.reaction"
        print rTemp[:,:].toarray()
        print "Difference: ", abs(self.reaction-rTemp).sum()
        print (self.reaction[:,:] - rTemp[:,:]).toarray()

    # Check that we have no NAN
    if(numpy.isnan(self.quantity).__array__().any()):
      print "Found a NaN, pre step"
      pdb.set_trace()


    # Step quantities one step using selected solver
    Solver.method[self.solverType](self)

    # Check that we have no NAN
    if(numpy.isnan(self.quantity).__array__().any()):
      print "Found a NaN, post step"
      pdb.set_trace()

    self.clock.tick()

    if((not self.disableOutput) and self.output):
        self.output.write()

    #self.printSumQuantity()
    #pdb.set_trace()

    #print "Cur time: ", self.clock.curTime()

    # self.printComps()

    # !!! The radie should increase when neurite grows, 
    #     not implemented currently.

  ##############################################################################

  def printSumQuantity(self):
 
    totLen = 0
    for comp in self.compartments:
      cl = comp.length()
      if(cl):
        totLen += comp.length()

    if(self.growthCones):
      QperL = self.growthCones[0].tubulinQuantityPerLength
      print "Quantity sum: ", sum(self.quantity)+totLen*QperL

  ##############################################################################

  # This function allows users to specify a condition at which the
  # simulation stops, making it possible to stop when for instance
  # the growth cones have reached a certain distance from the soma.

  def stop(self):

    for stopCond in self.stopCondition:
      if(stopCond()):
        self.stopCondition.remove(stopCond)
        return True

    return False

  ##############################################################################

  def finish(self):
    if(self.verbose):
      print("Closing open files...")

    self.output.close()

  ##############################################################################

  def addStopCondition(self, stopCondition):
    self.stopCondition.append(stopCondition)

  ##############################################################################

  def clearStopConditions(self):
    self.stopCondition = []

  ##############################################################################

  # This clamps the concentration in the growth cones and the soma
  # and iterates N times to find the correct concentration gradient.
  #
  #
  # Make sure we really want to remove transients.
  # If we are reloading then we do not want to do it.

  def removeTransients(self, nIter=None, disableOutput=True):

    if(nIter == None):
      print "Prematurely exited removeTransients."
      return

    # Check the max arc length of a neurite and figure 

    maxArcLen = 0
    minDiffusion = 1

    for subs in Substance._substanceList.itervalues():
      minDiffusion = min(minDiffusion,subs.diffusionConstant)

    for comp in self.compartments:
      maxArcLen = max(maxArcLen,comp.arcLength())

    # Estimate how long time it takes for diffusion to reach the entire
    # neuron. We also have to consider the number of compartments.
    if(not nIter):
      nIter = int(2*max(self.maxNumElement, \
                        math.ceil(maxArcLen**2/(4.0*minDiffusion) \
                        / self.clock.dt)))
    
      print("Estimated needing " + str(nIter) \
            + " iterations to remove transients")

    self.setDisableGrowth(True)
    self.disableOutput = disableOutput

    # Get the soma and growth cone concentrations
    clampSubstances = []
    clampValues = []

    for comp in self.compartments:
      if(comp.isSoma()):
        for subs in comp.substances.itervalues():
          clampSubstances.append(subs)
          clampValues.append(subs.quantity)

      if(comp.isGrowthCone):
        for subs in comp.substances.itervalues():
          clampSubstances.append(subs)
          clampValues.append(subs.quantity)

    # Loop nIterTransientRemoval steps to remove transients
    for i in range(0,nIter):
      self.step()

      # At every step reset the concentration back to clamped values
      for j in range(0,len(clampSubstances)):
        clampSubstances[j].quantity = clampValues[j]
        

    # Reset clock
    self.clock._t = self.clock._start

    self.setDisableGrowth(False)
    self.disableOutput = False

  ##############################################################################

  def findInitState(self):
 
    nIter = 0

    if(Solver.nIterTransientRemoval > 0):
      nIter = Solver.nIterTransientRemoval

    elif(Solver.nIterTransientRemoval == -1):
      # Try and estimate the number of iterations needed based on
      # diffusion speed of first molecule.
      
      maxArcLen = 0
      minDiffusion = 1

      for subs in Substance._substanceList.itervalues():
        minDiffusion = min(minDiffusion,subs.diffusionConstant)

      for comp in self.compartments:
        maxArcLen = max(maxArcLen,comp.arcLength())

      nIter = math.ceil(maxArcLen**2/(4.0*minDiffusion) / self.clock.dt)
      print("Estimated needing " + str(nIter) + " iterations to remove transients")

    if(nIter > 0):

      if(self.verbose):
        print('--- Removing transients')
        print("Doing " + str(nIter) + " iterations to remove transients.")

      self.setDisableGrowth(True)
      self.disableOutput = True

      # qOld = self.quantity
      # self.step()
      # 
      # qDist = self.vectNorm(qOld-self.quantity)/self.vectNorm(self.quantity) 
      # 
      # while(qDist > 1e-3):
      #  qOld = self.quantity
      #  self.step()
      #
      #  qDist = self.vectNorm(qOld-self.quantity)/self.vectNorm(self.quantity)
      #  print('qDist: ', qDist)

      # Loop nIterTransientRemoval steps to remove transients
      for i in range(0,nIter):
        self.step()

      if(self.verbose):
        print('--- Transients removed (maybe)')

      # Reset clock
      self.clock._t = self.clock._start

      self.setDisableGrowth(False)
      self.disableOutput = False

  ##############################################################################

  # We try and set a concentration gradient which is closer to the true
  # steady state value. This is to make findInitStateAuto converge faster

  def setEstimatedGradient(self,subsName,concSoma,concGC):

    print "Guessing concentration gradient"


    for gc in self.growthCones:
      neuriteLength = gc.arcLength()
      gc.substances[subsName].conc = concGC

      comp = gc.parent

      while(not comp.isSoma()):
        
        concComp = concGC \
              + (concSoma-concGC)*math.exp(-comp.arcLength()/neuriteLength)

        comp.substances[subsName].conc = concComp

        comp = comp.parent

      if(comp.isSoma()):
        comp.substances[subsName].conc = concSoma       

    print "Done with first gradient guestimate"
    # Note we overwrite the conc values proximally of a branch point
    # several times. This can lead to non-monotonic gradients if the
    # find init state does not run enough times

  ##############################################################################

  def setDisableGrowth(self,flag):
    self.disableGrowth = flag

    # The matrixes need to be recalculated
    self.setProduction()
    self.calculateTransportReactionMatrix()
     
  ##############################################################################

  def findInitStateAuto(self,maxDiff=1e-4,minIter=1000,maxIter=100000):

    print "Searching for init-state"

    self.setDisableGrowth(True)
    self.disableOutput = True

    oldProd = self.production
    # self.setProduction()
    self.production = numpy.zeros((self.maxNumElement,1))
    self.calculateTransportReactionMatrix()

    # Get the soma and growth cone concentrations
    clampSubstances = []
    clampValues = []

    for comp in self.compartments:
      if(comp.isSoma()):
        for subs in comp.substances.itervalues():
          clampSubstances.append(subs)
          clampValues.append(subs.quantity)

      if(comp.isGrowthCone):
        for subs in comp.substances.itervalues():
          clampSubstances.append(subs)
          clampValues.append(subs.quantity)


    qDiff = 1
    iter = 0


    # Temporarilly disable step action during init of steady state
    oldStepAction = self._stepAction
    self._stepAction  = []

    while((qDiff > maxDiff or iter < minIter) and iter < maxIter):

      qOld = self.quantity

      self.step()
      iter += 1

      # At every step reset the concentration back to clamped values
      for j in range(0,len(clampSubstances)):
        clampSubstances[j].quantity = clampValues[j]

      qDiff = 0
      for i in range(0,self.quantity.shape[0]):
        if(qOld[i,0] > 0):
          qDiff = max(qDiff, abs((self.quantity[i,0]-qOld[i,0])/qOld[i,0]))
     
      # pdb.set_trace()

       
      if(iter % 1000 < 1e-3):
        print "Iter ", iter, " max relative diff: ", qDiff

      #if(iter > 30000):
      #  print "Many many iterations... why?"
      #  pdb.set_trace()

    print "Iter ", iter, " max relative diff: ", qDiff

    # Restore step action for real simulation
    self._stepAction = oldStepAction

    # Reset clock
    self.clock._t = self.clock._start

    self.setDisableGrowth(False)
    self.disableOutput = False

    self.production = oldProd

  ##############################################################################

  def vectNorm(self,spvect):
    tmp = spvect.transpose()*spvect
    return tmp.sum()

  ##############################################################################

  def run(self):

    print("Current time: " + str(self.clock.curTime()) \
            + " End time: " + str(self.clock.end))
    while((not self.clock.done()) and (not self.stop())):
      self.step()

    self.clearStopConditions()

  ##############################################################################

  def stepEuler(self):

    if(self.clock.curTime() == 0):
      if(self.verbose):
        print('Using Euler forward.')

    self.quantity += ( self.transport * self.quantity \
                       + self.reaction * self.quantity \
                       + self.production ) * self.clock.dt

  ##############################################################################

  def stepRK2(self):

    if(self.clock.curTime() == 0):
      if(self.verbose):
        print('Using RK2.')

    midStep = self.quantity \
               + ( self.transport * self.quantity \
                   + self.reaction * self.quantity \
                   + self.production) * self.clock.dt/2.0

    midRate = (self.transport * midStep \
               + self.reaction * midStep \
               + self.production)

    self.quantity = self.quantity + midRate * self.clock.dt

  ##############################################################################

  def stepRK4(self):

    if(self.clock.curTime() == 0):
      if(self.verbose):
        print('Using RK4.')

    k1 = ( self.transport * self.quantity \
           + self.reaction * self.quantity \
           + self.production ) 

    k2 = ( self.transport * (self.quantity + 0.5*self.clock.dt*k1) \
           + self.reaction * (self.quantity + 0.5*self.clock.dt*k1) \
           + self.production )  

    k3 = ( self.transport * (self.quantity + 0.5*self.clock.dt*k2) \
           + self.reaction * (self.quantity + 0.5*self.clock.dt*k2) \
           + self.production )  

    k4 = ( self.transport * (self.quantity + self.clock.dt*k3) \
           + self.reaction * (self.quantity + self.clock.dt*k3) \
           + self.production )  

    self.quantity = self.quantity \
                   + self.clock.dt*(k1 + 2*k2 + 2*k3 + k4)/6.0


  ##############################################################################

  def implicitEuler(self):

    if(self.clock.curTime() == 0):
      if(self.verbose):
        print('Using backward euler.')

       
    b = self.quantity+self.production*self.clock.dt

    if(self.useSparse):
      A = sparse.identity( self.maxNumElement ) \
           - ( self.transport + self.reaction )*self.clock.dt

      self.quantity = numpy.mat(spsolve(A,b)).transpose()
    else:
      A = numpy.eye( self.maxNumElement ) \
           - ( self.transport + self.reaction )*self.clock.dt

      self.quantity = numpy.linalg.solve(A,b)

  ##############################################################################

  # Dictionary for methods

  method = { 'euler': stepEuler, 
             'RK2': stepRK2,
             'RK4': stepRK4,
             'impEuler': implicitEuler }


  ##############################################################################

  def increaseMaxNumElement(self,maxNumElement):

    if(maxNumElement > self.maxNumElement):

      if(self.verbose):      
        print "Increasing maxNumElement to ", maxNumElement

      oldMax = self.maxNumElement
      nPad = maxNumElement - oldMax
      # Update the max num element counter
      self.maxNumElement = maxNumElement

      oldQuantity = self.quantity
      oldTransport = self.transport
      oldReaction = self.reaction
      oldProduction = self.production


      paddingVector = numpy.zeros((nPad,1))
      self.quantity   = numpy.bmat([[oldQuantity],[paddingVector]])
      self.production = numpy.bmat([[oldProduction],[paddingVector]])

      if(self.useSparse):
        paddingMat = csr_matrix((nPad,nPad))

        self.transport  = sparse.bmat([[oldTransport,None], \
                                      [None,paddingMat]],'csr')
        self.reaction   = sparse.bmat([[oldReaction,None], \
                                      [None,paddingMat]],'csr')
      else:
        padA = numpy.zeros((oldMax,nPad))
        padB = numpy.zeros((nPad,oldMax))
        paddingMat = numpy.zeros((nPad,nPad))

        self.transport = numpy.bmat([[oldTransport,padA],[padB,paddingMat]])
        self.reaction = numpy.bmat([[oldReaction,padA],[padB,paddingMat]])


  ##############################################################################

  # Calculate the movement matrix (diffusion and active transport)

  def calculateTransportReactionMatrix(self):
    #self.transport[:,:] = 0 # Slow!! faster to just reallocate... hmmm..
    if(self.useSparse):
      self.transport = csr_matrix(self.transport.shape)
      self.reaction = csr_matrix(self.reaction.shape)
    else:
      self.transport = numpy.zeros(self.transport.shape)
      self.reaction = numpy.zeros(self.reaction.shape)

    for comp in self.compartments:
      self.setTransportReactionMatrixOneSided(comp)

    self.clearModified()

    # pdb.set_trace()

  ##############################################################################

  def removeCompartmentTransportReaction(self,comp):

    for subs in comp.substances.itervalues():
      if(comp.parent):
        parSubsID = comp.parent.substances[subs.name].id

        # The if checks are here because if we are using sparse
        # matrices then it is very expensive to add a new element, 
        # so would be stupid to add a new element only to make it a 0

        if(self.transport[parSubsID,subs.id]):
          self.transport[parSubsID,subs.id] = 0
  
        if(self.transport[subs.id,parSubsID]):
          self.transport[subs.id,parSubsID] = 0

        if(self.reaction[parSubsID,subs.id]):
          self.reaction[parSubsID,subs.id] = 0

        if(self.reaction[subs.id,parSubsID]):
          self.reaction[subs.id,parSubsID] = 0

      for child in comp.children:
        childSubsID = child.substances[subs.name].id

        if(self.transport[childSubsID,subs.id]):
          self.transport[childSubsID,subs.id] = 0

        if(self.transport[subs.id,childSubsID]):
          self.transport[subs.id,childSubsID] = 0

        if(self.reaction[childSubsID,subs.id]):
          self.reaction[childSubsID,subs.id] = 0

        if(self.reaction[subs.id,childSubsID]):
          self.reaction[subs.id,childSubsID] = 0

    if(comp.isGrowthCone and comp.growthSubstance):
      gID = comp.growthSubstance.id
      tID = comp.substances["tubulin"].id     

      self.reaction[tID,tID] = 0
      self.reaction[gID,gID] = 0
      self.production[tID,0] = 0
      self.production[gID,0] = 0


  ##############################################################################

  # Updated so it only clears the necessary values, modifying sparse is slow

  def clearCompartmentTransport(self,comp):

    for subs in comp.substances.itervalues():
      if(comp.parent):
        parSubsID = comp.parent.substances[subs.name].id

        if(self.transport[parSubsID,subs.id]):
          self.transport[parSubsID,subs.id] = 0

        if(self.transport[subs.id,parSubsID]):
          self.transport[subs.id,parSubsID] = 0

        if(self.transport[subs.id,subs.id]):
          self.transport[subs.id,subs.id] = 0

        if(self.transport[parSubsID,parSubsID]):
          self.transport[parSubsID,parSubsID] = 0

      for child in comp.children:
        childSubsID = child.substances[subs.name].id

        if(self.transport[childSubsID,subs.id]):
          self.transport[childSubsID,subs.id] = 0

        if(self.transport[subs.id,childSubsID]):
          self.transport[subs.id,childSubsID] = 0

        if(self.transport[subs.id,subs.id]):
          self.transport[subs.id,subs.id] = 0

        if(self.transport[childSubsID,childSubsID]):
          self.transport[childSubsID,childSubsID] = 0

  ##############################################################################

  # This function is used to make sure conservation of molecules
  # is correct (used when a neighbouring compartment has been modified)
  def _transportMatrixConservation(self,comp,subsName):
    subs = comp.substances[subsName]

    if(self.transport[subs.id,subs.id]):
      self.transport[subs.id,subs.id] = 0;

    self.transport[subs.id,subs.id] = -self.transport[:,subs.id].sum()


  ##############################################################################

  def _setTransportMatrix(self,compA,compB,subsName):

    if(not compA):
      # No parent, nothing to do
      return

    if(compA.volume() <= 0 or compB.volume() <= 0):
      print "Volume is zero or negative!"
      pdb.set_trace()

    # This function always assumes that compA is the parent
    area = compB.crossSectionAreaWithParent()
    dist = compB.centerDistanceToParent()

    subsA = compA.substances[subsName]
    subsB = compB.substances[subsName]

    diffConst = subsA.diffusionConstant

    f = diffConst*area/dist

    ATrate = subsA.activeTransportRate*area/compA.volume()

    # Flow from A to B
    self.transport[subsB.id,subsA.id] = \
      self.transport[subsB.id,subsA.id] + f/compA.volume() + ATrate
    self.transport[subsA.id,subsA.id] = \
      self.transport[subsA.id,subsA.id] - f/compA.volume() - ATrate

    # Flow from B to A
    self.transport[subsA.id,subsB.id] = \
      self.transport[subsA.id,subsB.id] + f/compB.volume()    
    self.transport[subsB.id,subsB.id] = \
      self.transport[subsB.id,subsB.id] - f/compB.volume()    
    
  ##############################################################################

  def _setReactionProductionMatrix(self,comp):
    for subs in comp.substances.itervalues():
      # Degradation
      self.reaction[subs.id,subs.id] = -subs.degradationConstant

    if(comp in self.disabledPersistentGrowthCones):
      print "Comp growth disabled, in set reaction production matrix"
      # pdb.set_trace()

    if(comp.isGrowthCone and comp.growthSubstance):
      if(not self.disableGrowth): # Global setting, used for init state

        # print "GC setting reaction production matrix"
        # pdb.set_trace()

        # We can have a problem with presistent growth cones when 
        # the compartment gets too small, and the concentration tells it
        # to keep shrinking. Then we are not going to allow shrinkage.
        # 
        if(comp.persistent \
           and comp.substances["tubulin"].conc < (comp.tubulinDepolymerisation \
                                                  /comp.tubulinPolymerisation) \
           and comp.length() <= self.minCompartmentLength):

          print "Warning GC is small, has low conc and is persistent: ", \
                comp.substances["tubulin"].conc

          aID = comp.substances["tubulin"].id
          bID = comp.parent.substances["tubulin"].id
          print "transport = ", self.transport[aID,bID]

          #

          gID = comp.growthSubstance.id
          tID = comp.substances["tubulin"].id

          # Temporarilly make growth cone dormant, when the
          # concentration of tubulin is low, and the length is short
          # and the compartment is persistent.

          # self.reaction[tID,tID] -= 0
          self.reaction[gID,tID] = 0

          self.production[tID,0] = 0
          self.production[gID,0] = 0

          # We have to make sure that these formulas gets overwritten
          # once the concentration starts rising again

          if(not comp in self.disabledPersistentGrowthCones):
            self.disabledPersistentGrowthCones.append(comp)
            print "Shutting down growth cone"

        else:

          # No need to check it every step anymore
          if(comp.persistent and comp in self.disabledPersistentGrowthCones):
             self.disabledPersistentGrowthCones.remove(comp)
             print "Allowing growth cone to grow again"

          # Calculating how much of the tubulin quantity is converted to growth
          gID = comp.growthSubstance.id
          tID = comp.substances["tubulin"].id

          self.reaction[tID,tID] -= comp.tubulinQuantityPerLength \
                                    * comp.tubulinPolymerisation \
                                    / comp.volume()

          self.reaction[gID,tID] = comp.tubulinQuantityPerLength \
                                   * comp.tubulinPolymerisation \
                                   / comp.volume()

          self.production[tID,0] = comp.tubulinQuantityPerLength \
                                 * comp.tubulinDepolymerisation

          self.production[gID,0] = - comp.tubulinQuantityPerLength \
                                   * comp.tubulinDepolymerisation

      else:
        # Growth disabled, make sure the fields are cleared
        gID = comp.growthSubstance.id
        tID = comp.substances["tubulin"].id

        if(self.reaction[tID,tID]):
          self.reaction[tID,tID] = 0

        if(self.reaction[gID,tID]):
          self.reaction[gID,tID] = 0

        if(self.production[tID,0]):
          self.production[tID,0] = 0

        if(self.production[gID,0]):
          self.production[gID,0] = 0


      # Chemical reactions
      # !!! To be added

  ##############################################################################

  def setProduction(self):
    for comp in self.compartments:
      for subs in comp.substances.itervalues():
         self.production[subs.id,0] = subs.productionRate
  
  ##############################################################################

  # This does the reactions in the compartment, but only the 
  # transport to the parent
  def setTransportReactionMatrixOneSided(self,compartment):
    # This function does *NOT* clear the matrix before adding elements

    for name in compartment.substances.iterkeys():
      self._setTransportMatrix(compartment.parent,compartment,name)      

    self._setReactionProductionMatrix(compartment)

  ##############################################################################

  def setTransportMatrixTwoSided(self,compartment):

    self.clearCompartmentTransport(compartment)

    parent = compartment.parent
    children = compartment.children

    for name in compartment.substances.iterkeys():
      # Should we clear row first?
      if(parent):
        self._setTransportMatrix(parent,compartment,name)
        self._transportMatrixConservation(parent,name)

      for child in children:
        self._setTransportMatrix(compartment,child,name)
        self._transportMatrixConservation(child,name)

  ##############################################################################

  # Update the columns and rows that corresponds to the substances 
  # in this compartment
  def setTransportReactionMatrixTwoSided(self,compartment):

    # First we clear the old transports/reactions that might have been 
    self.removeCompartmentTransportReactions(compartment)

    parent = compartment.parent
    children = compartment.children

    for name in compartment.substances.iterkeys():
      # Should we clear row first?
      self._setTransportMatrix(parent,self,name)      
      self._transportMatrixConservation(parent,name)

      for child in children:
        self._setTransportMatrix(self,child,name)
        self._transportMatrixConservation(child,name)

    self._setReactionProductionMatrix(compartment)

  ##############################################################################

  def updateSpatialDiscretisation(self,compartment):

    if(compartment.isSoma()):
      # It is soma, do nothing
      return False

    ## Are the compartment too big
    # Used to be a while statement, that was ok when we split off small
    # parts because it allowed us to discretise everything in one iteration
    # but now that we divide compartment in the middle, there is no point
    # in doing a while statement here.
    if(compartment.length() > self.maxCompartmentLength):

      # We must split compartment
      if(self.verbose):
        print("Split compartment, length: " + str(compartment.length()))

      return self._splitCompartment(compartment) # Returns true if modified

    if(compartment.length() < self.minCompartmentLength):

      # We must merge compartments
      if(self.verbose):
        print("Merge compartment " + str(compartment) \
               + " length: " + str(compartment.length()))

      return self._mergeCompartments(compartment) # Returns true if modified

    # Recalculations of transport matrix is done after, called from step

  ##############################################################################

  def _splitCompartment(self,compartment):

    # This version always splits compartments in the middle

    # Disallow splitting of growth cones!
    if(compartment in self.growthCones):
      if(self.preventSplittingOfGrowthCones):
        print("Prevented splitting of growth cone!")
        return False
      elif(compartment.length() < 1.05*self.maxCompartmentLength):
        print("Prevented splitting of growth cone (but will do it eventually!)")
        print "GC length: " + str(compartment.length()) \
               + " maxCompLength: " + str(self.maxCompartmentLength)
        return False

    # This is necessary if we have split neighbours on both sides
    self.removeCompartmentTransportReaction(compartment)

    # Save for later use
    oldParent = compartment.parent

    newCompartmentLength = compartment.length()/2

    if(compartment.parent.isSoma()):
      newCompEndPoint = compartment.direction() \
          * (newCompartmentLength + compartment.parent.radie) \
          + compartment.parent.endPoint
    else:
      newCompEndPoint = compartment.direction()*newCompartmentLength \
                         + compartment.parent.endPoint

    newComp = Compartment(self,compartment.parent,newCompEndPoint, \
                          compartment.radie)

    # We need to remove the reaction and transport before we change the
    # children or parents
    self.removeCompartmentTransportReaction(compartment)

    # Update the structure so that it is 
    # oldParent->newCompartment->oldCompartment (smaller)

    compartment.parent.children.remove(compartment)
    compartment.parent.addChild(newComp)
    newComp.addChild(compartment)
    compartment.parent = newComp

    # If there is a gradient across the compartment then we want
    # to preserve it when we split the compartment, we do this by
    # assuming the gradient is linear between the parent and child
    # of the compartment.
    # Only do it if there are children, and if it is not the first time step
    # 
   
    for subs in compartment.substances.itervalues():

      if(len(compartment.children) > 0 \
         and not abs(self.clock._t - self.clock.start) < self.clock.dt/2.0):

        parentConc = oldParent.substances[subs.name].conc
        branchGrad = []

        for child in compartment.children:
          childConc = child.substances[subs.name].conc
          distParChild = (oldParent.center()-child.center()).norm()
          branchGrad.append((childConc-parentConc)/distParChild)

          if(self.verbose):
            print "Par conc: " + str(parentConc) \
                + " child conc: " + str(childConc) \
                + " PC center dist " + str(distParChild) \
                + " grad: " + str((childConc-parentConc)/distParChild)

        subsGrad = sum(branchGrad)/len(branchGrad) # Mean!

      else:

        # No child, so do not do gradient interpolation
        subsGrad = 0

      if(self.verbose):
        print "New comp length: " + str(newComp.length())
        print "Parent end point: " + str(newComp.parent.endPoint)
        print "New comp endpoint: " + str(newComp.endPoint)

      # Distance between centres of new and original old compartment
      halfLen = newComp.length()/2 

      if(self.verbose):
        print "subsGrad: " + str(subsGrad) + " halfLen " + str(halfLen)

      # Calculate the gradient, and split the quantities accordingly
      # gradient between the new compartments parent and child
      quantDiff = subsGrad*halfLen*compartment.volume()

      if(self.verbose):
        print "quantDiff " + str(quantDiff)

      # Need the zero check if we allow negative concentrations...
      # which we only do if we are testing validity of parameters.
      if(abs(quantDiff) > 0.5*subs.quantity and subs.quantity > 0):
        # Too large gradient
        print "WARNING: Too steep gradient"
        quantDiff *= 0.499999*subs.quantity/abs(quantDiff)

      if(subs.quantity > 0):
        newSubQuantity = 0.5*subs.quantity + quantDiff
        oldSubQuantity = 0.5*subs.quantity - quantDiff
      else:
        newSubQuantity = 0
        oldSubQuantity = 0
        

      if(self.noNegativeConcentrations \
         and (newSubQuantity < 0 or oldSubQuantity < 0)):
        print("Negative quantities. What are you doing wrong!!")
        pdb.set_trace()


      # Note the production rate is not split, it remains in parent!!
      Substance(subs.name, newComp, self, \
                newSubQuantity, \
                subs.diffusionConstant, \
                subs.degradationConstant, \
                subs.activeTransportRate, \
                subs.productionRate*0)



      # Remove the split fraction from parent
      subs.quantity = oldSubQuantity

      # Update the production matrix, however the new substance
      # does not have any production currently.
      # Perhaps the production should be inherited from parent substance

      newsub = newComp.substances[subs.name]
      self.production[newsub.id,0] = newsub.productionRate

    # Mark that these two compartments needs to be recalculated
    self.markModified(compartment)
    self.markModified(newComp)

    # Make sure that newComp is a good enough size...
    self.updateSpatialDiscretisation(newComp)

    return True

  ##############################################################################

  def _mergeCompartments(self,compartment):

    parent = compartment.parent

    print "Called merge compartments"

    # Disallow merging of growth cones!
    if(compartment in self.growthCones):

      # pdb.set_trace()

      if(self.allowGrowthConeMerge \
         and parent \
         and parent.isNeurite() \
         and (not parent.isBranchPoint())):

        # Special case, we are going to merge a growth cone. We have to
        # merge parent into this compartment. Normally the child is merged
        # into the parent, but here the child is special.

        self.removeCompartmentTransportReaction(parent)
        self.removeCompartmentTransportReaction(compartment)

        # Check that parent exists, that parent is a neurite and that
        # the parent only has the growth cone as child. Done.

        for (name,subs) in compartment.substances.iteritems():

          parentId = parent.substances[name].id 

          self.quantity[subs.id,0] += self.quantity[parentId,0]
          self.quantity[parentId,0] = 0

          # Remove production in deleted compartment
          self.production[parentId,0] = 0

        # Adjust the radie
        compartment.radie = math.sqrt((parent.length()*(parent.radie**2) \
                                        + compartment.length() \
                                           *(compartment.radie**2)) \
                                        / (parent.length() \
                                           + compartment.length()))


        # The compartment is only child of old parent, bypass old parent
        compartment.parent = parent.parent
        parent.parent.children.remove(parent)
        parent.parent.children.append(compartment)

        # Clean up references
        parent.parent = None
        parent.children = []

        # Remove old parent from solver
        self.compartments.remove(parent)

        # Removed last references to the compartment, wait for 
        # garbage collection to send it to compartment heaven.

        self.markModified(compartment)
        print "Merging growth cone, new length is ", compartment.length(), "m"

        # Indicate we updated the discretization
        return True

      else:
        print("Prevented merging of growth cone! Length: " + str(compartment.length()) + "m")

        if(compartment.volume() == 0):
          print "Growth cone volume is 0, Error!"
          pdb.set_trace()

        return False


    # Only merge with if parent exists and if it is not the soma
    # Modified, some SWC files have compartments with 0 volume... ouch!
    if(parent and parent.isNeurite() and (not parent.isBranchPoint())):

      # Clear row and column corresponding to substance in solver matrixes
      self.removeCompartmentTransportReaction(compartment)
      self.removeCompartmentTransportReaction(parent)

      # Adjust the radie
      parent.radie = math.sqrt((parent.length()*(parent.radie**2) \
                                 + compartment.length() \
                                   *(compartment.radie**2)) \
                                / (parent.length() + compartment.length()))

      # Elongate the parent compartment (assuming it is neurite)
      parent.endPoint = compartment.endPoint

      for (name,subs) in compartment.substances.iteritems():

        parentId = parent.substances[name].id 

        self.quantity[parentId,0] += self.quantity[subs.id,0]
        self.quantity[subs.id,0] = 0

        # Remove production in deleted compartment
        self.production[subs.id,0] = 0

        for child in compartment.children:

          # Bypass the compartment itself
          parent.addChild(child)
          child.parent = parent

        # Remove last references to the compartment and wait for 
        # garbage collection to send it to compartment heaven.
        parent.children.remove(compartment)
        self.compartments.remove(compartment)

        # Remove all references to neighbours from the removed compartment
        compartment.parent = None
        compartment.children = []

        self.markModified(parent)

        if(parent.volume() == 0):
          print "How is this possible..?"
          pdb.set_trace()

      return True # Updated discretization
    else:
      return False # Nothing updated

  ##############################################################################

  def printComps(self):
    for comp in self.compartments:
      comp.printComp()

  ##############################################################################

  def saveState(self, fileName):
    if(self.verbose):
      print "Saving state to ", fileName

    # We have to close the files before pickling
    if(self.output):
      self.output.close()    

    fp = open(fileName,'wb')
    pickle.dump(self,fp)
    pickle.dump(Substance,fp)
    fp.close()

    # Reopen the output file so that we can keep writing to it
    if(self.output):
      self.output.reopen()

  ##############################################################################

  def loadState(self, fileName):
    if(self.verbose):
      print "Load state from", fileName

    fp = open(fileName, 'rb')
    oldSolver = pickle.load(fp)
    oldSubstance = pickle.load(fp)

    # Need to restore the class variables in Substance
    Substance._substanceList = oldSubstance._substanceList

    fp.close()
    return oldSolver

  ##############################################################################

