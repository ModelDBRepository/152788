# This script writes the simulation parameters to a text file
# to help with book keeping.

class SaveInfo():

  def __init__(self, fileName, experiment):
    self.fileName = fileName
    self.experiment = experiment
    self.fp = None

  def __del__(self):
    if(self.fp):
      self.fp.close()

  # The description is a string prepended
  # to the info, for example "Pre change of polymerisation rate"
  # or "Post change of polymerisation rate"

  def writeInfo(self,description,appendFlag=False):

    if(appendFlag):
      if(not self.fp):
        self.fp = open(self.fileName,"a+")

    else:
      if(self.fp):
        self.fp.close()

      self.fp = open(self.fileName,"w")      

    # Write the actual info to file

    data = []
    data.append(description + '\n')
    data.append('tubulinQuantityPerLength ' \
                    + self.experiment.tubulinQuantityPerLength + '\n')
    data.append('tubulinDiffusionConstant ' \
                    + self.experiment.tubulinDiffusionConstant + '\n')
    data.append('tubulinDegradationConstant ' \
                    + self.experiment.tubulinDegradationConstant + '\n')
    data.append('tubulinActiveTransportRate ' \
                    + self.experiment.tubulinActiveTransportRate + '\n')
    data.append('tubulinSomaProductionRate ' \
                    + self.experiment.tubulinSomaProducationRate + '\n')
    data.append('neuriteGrowthPoly ' \
                    + self.experiment.neuriteGrowthPoly + '\n')
    data.append('neuriteGrowthDepoly ' \
                    + self.experiment.neuriteGrowthDepoly + '\n')


    self.fp.write(''.join(data))
     
