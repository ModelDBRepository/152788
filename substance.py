# Python debugger
import pdb
# Use pdb.set_trace() to enter debugger somewhere


class Substance(object):
  '''
  Implements a Substance
  '''

  _substanceList = {}
  
  # Production rate must be specified on a compartment to compartment basis
  # Active transport rate is parent to child

  def __init__(self, name, compartment, solver, quantity=0.0, \
               diffusionConstant=None, \
               degradationConstant=None, \
               activeTransportRate=None, \
               productionRate=0):

    self.name = name
    self.compartment = compartment

    if(compartment):
      compartment.addSubstance(self)

    self.solver = solver
    self.id = solver.getNewSubstanceId()
    self.quantity = quantity

    self.productionRate = productionRate

    # Check if the substance already exist, if not register it
    if self.name in Substance._substanceList:

      oldSubs = Substance._substanceList[name]

      if(diffusionConstant \
         and diffusionConstant != oldSubs.diffusionConstant):
        print("Prevented changing of diffusion constant from ", \
              oldSubs.diffusionConstant, " to ", diffusionConstant)

      if(degradationConstant \
         and degradationConstant != oldSubs.degradationConstant):
        print("Prevented changing of degradation constant from ", \
              oldSubs.degradationConstant, " to ", degradationConstant)

      # Use the original diffusion and degradation constants
      self.diffusionConstant = oldSubs.diffusionConstant
      self.degradationConstant = oldSubs.degradationConstant

      if(activeTransportRate):
        self.activeTransportRate = activeTransportRate
      else:
        self.activeTransportRate = oldSubs.activeTransportRate

    else:

      if(not diffusionConstant):
        print("No diffusion constant specified for ", name, " using 0")
        diffusionConstant = 0

      if(not degradationConstant):      
        print("No degradation constant specified for ", name, " using 0")
        degradationConstant = 0

      if(not activeTransportRate):
        activeTransportRate = 0

      self.diffusionConstant = diffusionConstant
      self.degradationConstant = degradationConstant
      self.activeTransportRate = activeTransportRate

      Substance._substanceList[self.name] = self


  # Quantity is stored in the solver, we access the value there when needed

  def get_quantity(self):
    return self.solver.quantity[self.id,0]

  def set_quantityOLD(self, quantity):
    try:
      self.solver.quantity[self.id,0] = quantity
    except IndexError:
      print "What is going on ... substance.py"
      pdb.set_trace()

  def set_quantity(self, quantity):
    self.solver.quantity[self.id,0] = quantity

  quantity = property(fget=get_quantity,fset=set_quantity)


  # We do not store concentration explicity, instead we calculate
  # it from the quantity and volume

  def get_conc(self):
    return self.quantity/self.compartment.volume()

  def set_conc(self,conc):
    self.quantity = conc*self.compartment.volume()

  conc = property(fget=get_conc, fset=set_conc)


  def printSubs(self):
    print(self.name, "conc: ", self.conc, " quant: ", self.quantity)


  def clearSolver(self):
    self.solver = None

  def setSolver(self,solver):
    self.solver = solver
