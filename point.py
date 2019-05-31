# Helper class that takes care of the 

from math import sqrt

class Point():

  coords = []

  def __init__(self, coords):
    self.coords = coords

  def __add__(self, other):
    return Point((self.coords[0] + other.coords[0], \
                  self.coords[1] + other.coords[1], \
                  self.coords[2] + other.coords[2]))

  def __sub__(self, other):
    return Point((self.coords[0] - other.coords[0], \
                  self.coords[1] - other.coords[1], \
                  self.coords[2] - other.coords[2]))

  # Multiplication with scalar only works if scalar is written after
  def __mul__(self, other):
    return Point((self.coords[0]*float(other), \
                  self.coords[1]*float(other), \
                  self.coords[2]*float(other)))

  def __div__(self, other):
    return Point((self.coords[0]/float(other), \
                  self.coords[1]/float(other), \
                  self.coords[2]/float(other)))

  # !!! Should I use this, or does it screw up picle
  def __repr__(self):
    return "(" + str(self.coords[0]) + "," \
               + str(self.coords[1]) + ","\
               + str(self.coords[2]) + ")"

  def __getitem__(self,key):
    return self.coords[key]

  def norm(self):
    return sqrt(self.coords[0]**2 + self.coords[1]**2 + self.coords[2]**2)

  
