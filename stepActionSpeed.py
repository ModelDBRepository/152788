from solver import Solver
from point import Point
from experiment import Experiment
from substance import Substance
from slavedGrowthConeSpeed import SlavedGrowthConeSpeed
from compartment import Compartment


import pdb


# This function can be added to solver, will then be run every step

class StepActionSpeed():

  def __init__(self, solver, growthConeData, predictGCnumber):

    self.solver = solver
    solver.addStepAction(self)

    self.predictGCflag = [0]*len(growthConeData)

    # There must be a nicer way of doing this...
    for idx in predictGCnumber:
      self.predictGCflag[idx-1] = 1

    self.growthConeInfo = growthConeData

    self.slavedGrowthCones = []

  def update(self):
    # This checks if there is time to add a new growth cone, and if so
    # it finds the closest compartment.
    
    for gcInfo in self.growthConeInfo:
      time = gcInfo[0]
      somaDist = gcInfo[1]
      allSpeeds = gcInfo[2]
      predictFlag = gcInfo[3]
      arclength = gcInfo[4]

      if(self.solver.verbose):
        print "Solver clock: " + str(self.solver.clock.curTime()) \
            + " GC arrive time " + str(time[0]) + "\n"

      if(self.solver.clock.curTime() >= time[0]):
        print "Time to add GC " + str(self.solver.clock.curTime())

        # Find closest compartment to spawn growth cone from
        minDist = float('infinity')
        closestCompartment = None

        for comp in self.solver.compartments:

          if(comp in self.solver.growthCones):
            # Do not attach to another growth cone
            continue

          compDist = abs(somaDist - comp.arcLength())

          if(compDist < minDist):
            minDist = compDist
            closestCompartment = comp
            # Create a 1 micrometer long new growth cone
            coords = comp.endPoint + Point((0.9e-6,0,0))

        print "Asked dist: " + str(somaDist)
        print "Start dist: " + str(closestCompartment.arcLength())
        print "Min distance: " + str(minDist)

        # Add growth cone to compartment
        if(predictFlag):
          growthCone = Compartment(self.solver, \
                                   closestCompartment, \
                                   coords, \
                                   Experiment.neuriteRadie)

          growthCone.makeGrowthCone(Experiment.neuriteGrowthPoly, \
                                    Experiment.neuriteGrowthDepoly, \
                                    Experiment.tubulinQuantityPerLength)

          growthCone.persistent = True

        else:
          print "Adding slave cone"
          # Slave mode
          growthCone = SlavedGrowthConeSpeed(self.solver, \
                                             closestCompartment, \
                                             None, \
                                             coords, \
                                             Experiment.neuriteRadie, \
                                             Experiment.tubulinQuantityPerLength)
          # print "WARNING: Temporarilly setting tubulinQuantityPerLength to zero!!!"

          growthCone.setSlaveInfo(time,allSpeeds,arclength)

          self.slavedGrowthCones.append(growthCone)

        quantityRatio = growthCone.volume() \
            / (growthCone.volume() + closestCompartment.volume())

        for subs in closestCompartment.substances.itervalues():
          # Add substance to new compartment, steal some from parent
          Substance(subs.name, growthCone, self.solver, \
                      quantityRatio*subs.quantity)
          subs.quantity = subs.quantity * (1 - quantityRatio)

        # Remove growth cone info
        self.growthConeInfo.remove(gcInfo)
 
        # Verify that the old GC is removed...
        # pdb.set_trace() 


  def init(self):
    for gc in self.slavedGrowthCones:
      gc.init()
