from solver import Solver

class StepActionSomaClamp():

  def __init__(self, solver, soma, substanceName, concentration):

    self.solver = solver
    self.soma = soma
    self.substanceName = substanceName
    self.concentration = concentration

  def update(self):

    self.soma.substances[self.substanceName].conc = self.concentration

  def init(self):
    print "Initialising Step Action Soma Clamp"
    self.soma.substances[self.substanceName].conc = self.concentration
