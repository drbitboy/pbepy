import numpy
class KINETXBASE:
  """Base class for KinetX data (delivery and knowledge uncertainties"""

  ### To be overridden in systems with a primary and multiple satellites
  PRIMARY_BARYCENTER = False

  def __init__(self): pass

  def kget(self
          ,KnSig=False      ### Get Knowledge sigmas
          ,KnEig=False      ### Get Knowledge eigenvalues (sigmas^2)
          ,DeSig=False      ### Get Delivery sigmas
          ,DeEig=False      ### Get Delivery eigenvalues (sigmas^2)
          ):

    if KnSig:                                       ### Knowledge sigmas
      try   : return self.sigmas
      except: return numpy.sqrt(self.eigVals)

    if KnEig:                                  ### Knowledge eigenvalues
      try   : return self.eigVals
      except: return self.sigmas * self.sigmas


    if DeSig:                                        ### Delivery sigmas
      try   : return self.delivSigmas
      except:
        try   : return numpy.sqrt(self.delivEigVals)
        except:
          msg = ('### WARNING:  DELIVERY (CONTROL) VALUES NOT AVAILABLE'
                + '; USING KNOWLEDGE VALUES * [1.5,2.25,3]'
                )
          print(msg)
          scal = numpy.array([150e0,150,100])/numpy.array([100e0,66,33])
          return self.kget(KnSig=True) * scal

    if DeEig:                                   ### Delivery eigenvalues
      try   : return self.delivEigVals
      except:
        try   : return self.delivSigmas * self.delivSigmas
        except:
          msg = ('### WARNING:  DELIVERY (CONTROL) VALUES NOT AVAILABLE'
                + '; USING KNOWLEDGE VALUES * [1.5,2,25,3]'
                )
          print(msg)
          scal = numpy.array([150e0,150,100])/numpy.array([100e0,66,33])
          return self.kget(KnEig=True) * scal * scal

    return self.kget(KnSig=True)   ### Default:  return Knowledge sigmas
