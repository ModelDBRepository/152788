import numpy 
import substance

class SaveData():
  '''
  Saves morphological structure to file.
  '''

  def __init__(self, fileName, solver, outputDt, saveFlux=False):
    self.solver = solver
    self.fileName = fileName
    self.fp = None
    self.outputDt = outputDt

    # Should we store net flux in from parent?
    self.saveFlux = saveFlux


  def __del__(self):
    if(self.fp):
      self.fp.close()


  def init(self):

    header = []
    header.append('time;id;parent id;start coords;end coords;radie;dist')

    for name in substance.Substance._substanceList:
      header.append(';')
      header.append(str(name))

      if(self.saveFlux):
        header.append(";flux")


    header.append('\n')

    self.fp = open(self.fileName,"w")

    self.fp.write(''.join(header))


  def close(self):
    if(self.fp):
      print "Closing ", self.fileName
      self.fp.close()
      self.fp = None
    else:
      print "File ", self.fileName, " already closed."

  def reopen(self):
    self.fp = open(self.fileName,"a+b")

  def write(self):

    if(self.solver.clock.curTime() % self.outputDt \
        >= self.solver.clock.dt / (2 * self.outputDt)):

     # print "Skipping write: ", self.solver.clock.curTime()

     # Not time to write output, skip
     return

    if(self.solver.verbose):
      print "* Write called: ", self.solver.clock.curTime()

    for comp in self.solver.compartments:
      line = []
      line.append(str(self.solver.clock.curTime()))
      line.append(';')
      # line.append(str(id(comp)))
      line.append(str(comp.id))      
      line.append(';')
      if(comp.parent):
        # line.append(str(id(comp.parent)))
        line.append(str(comp.parent.id))        
        line.append(';')
        if(comp.parent.parent):
          # Start point is parents end point
          line.append("(%.5g,%.5g,%.5g)" % \
                      (comp.parent.endPoint[0], \
                       comp.parent.endPoint[1], \
                       comp.parent.endPoint[2]))
        else:
          # If the parent is soma, start point is outside the sphere
          cDir = comp.direction()

          line.append("(%.5g,%.5g,%.5g)" % \
                      (comp.parent.endPoint[0]+comp.parent.radie*cDir[0], \
                       comp.parent.endPoint[1]+comp.parent.radie*cDir[1], \
                       comp.parent.endPoint[2]+comp.parent.radie*cDir[2]))
        line.append(';')
      else:
        # No parent, use ID -1 and use own end point for start point
        line.append('-1')
        line.append(';')
        line.append("(%.5g,%.5g,%.5g)" % \
                  (comp.endPoint[0], \
                   comp.endPoint[1], \
                   comp.endPoint[2]))
        line.append(';')      

      # line.append(str(comp.endPoint))
      line.append("(%.5g,%.5g,%.5g)" % \
                  (comp.endPoint[0], comp.endPoint[1], comp.endPoint[2]))
      line.append(';')
      line.append("%.5g" % (comp.radie))
      line.append(';')
      # False means do not use recursion - faster.
      #line.append("%.5g" % (comp.arcLength(False)))
      line.append("%.5g" % (comp.arcLength()))

      # We use the sim substance to make sure they are in the right order
      for name in substance.Substance._substanceList.iterkeys():
        line.append(";")
        line.append(str(comp.substances[name].conc))

        if(self.saveFlux):
          # We are also saving net influx from parent
          selfSubs = comp.substances[name]

          if(comp.parent):
            parentSubs = comp.parent.substances[name]

            inFlux = self.solver.transport[selfSubs.id,parentSubs.id] \
                      * parentSubs.quantity
            outFlux = self.solver.transport[parentSubs.id,selfSubs.id] \
                      * selfSubs.quantity
          else:
            # No parent, no influx...
            inFlux = 0.0
            outFlux = 0.0

          line.append(";")
          line.append("%.5g" % (inFlux-outFlux))

        #line.append(str(comp.substances[name].quantity))

      line.append('\n')

      self.fp.write(''.join(line))
