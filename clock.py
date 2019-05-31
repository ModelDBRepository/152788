# Keeps track of the simulation time

from time import time

class Clock():
  '''
  Simulation clock
  '''

  def __init__(self,dt=1.0,start=0.0,end=20.0):
    self._dt = float(dt)
    self._t = float(start)
    self._start = float(start)
    self._end = float(end)

  def reset(self):
    self._t = float(self._start)

  def tick(self):
    self._t += self._dt


  def set_dt(self,dt):
    self._dt = float(dt)

  def get_dt(self):
    return self._dt

  dt = property(fget=get_dt,fset=set_dt)


  def set_start(self,start):
    self._start = float(start)

  def get_start(self):
    return self._start

  start = property(fget=get_start,fset=set_start)


  def set_end(self,end):
    self._end = float(end)

  def get_end(self):
    return self._end

  end = property(fget=get_end,fset=set_end)


  def done(self):
    return (self._t >= self._end)

  def curTime(self):
    return self._t


  
  
