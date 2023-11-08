import numpy
import spiceypy as sp


########################################################################
class DIRECTEDPLANE:
  """SPICE plane wrapper that maintains input normal's direction"""
  def __init__(self,normal,plane):
    ###print(dict(n=normal,p=plane,r='DIRECTEDPLANE.__init__'))
    self.init_normal,self.plane = normal,plane

  @staticmethod
  def nvc2pl(normal,konst):
    ###print(dict(n=normal,k=konst,r='DIRECTEDPLANE.nvc2pl'))
    return DIRECTEDPLANE(normal,sp.nvc2pl(normal,konst))
  @staticmethod
  def nvp2pl(normal,point):
    ###print(dict(n=normal,p=point,r='DIRECTEDPLANE.nvp2pl'))
    ###print(dict(plane=sp.nvp2pl(normal,point),r='DIRECTEDPLANE.nvp2pl'))
    return DIRECTEDPLANE(normal,sp.nvp2pl(normal,point))

  def pl2nvc(self):
    normal,konst = sp.pl2nvc(self.plane)
    if sp.vdot(self.init_normal,normal) < 0.0: normal *= -1.0
    return normal,konst
  def pl2nvp(self):
    ###print(self.plane)
    normal,point = sp.pl2nvp(self.plane)
    if sp.vdot(self.init_normal,normal) < 0.0: normal *= -1.0
    return normal,point
  def csppl(self): return sp.nvp2pl(*sp.pl2nvp(self.plane))

########################################################################
### Find surface point of ellipse where surface normal is parallel and
###   opposite in direction to an input plane (argument plane) normal.
###   The A, B & C semi-axes of the ellipse are in argument abcArg.
###
###   The return structure contains the ellipse surface point and the
###   plane height above that surface point.  The sign of the height
###   is the same as the sign of the dot product of the ellipse surface
###   normal and a vector to the plane from the ellipse surface point.
###   Put another way, the height is zero if the plane is tangent to
###   the surface, it is positive if the plane is above (in the
###   direction of the surface normal from) the ellipse, otherwise
###   it is negative.
########################################################################
class EllSurfTangent:

  def __init__(self):
    self.converged = False                         ### True if converged
    self.status = 'EllSurfTangent:  uninitialized' ### Text status
    self.msg = ''                                  ### Error message

  def calculate(self, abcArg, dirplane
               ,epsilon=1e-7
               ,maxIter=20
               ):
    try:
      self.converged = False
      self.abc = numpy.abs(numpy.array(abcArg,dtype=numpy.float64))
      self.dirplane = dirplane
      self.epsilon = epsilon
      self.maxIter = maxIter

      ### N.B. uvPlane normal points inward i.e. towards negative height
      uvPlane,pPlane = dirplane.pl2nvp()

      outvec = -10e3 * self.abc.max() * uvPlane
      eps2 = self.epsilon * self.epsilon
      self.iter = 0
      surfpt = outvec
      err2 = (eps2 * 2e0) + 2e0
      iwct = 3
      while err2 > eps2 and self.iter < self.maxIter and iwct != 0:
        self.iter += 1
        lastSurfpt = surfpt
        surfpt,alt = sp.nearpt(sp.vadd(lastSurfpt,outvec), *self.abc)
        uvNorm = sp.surfnm(*self.abc,surfpt)
        crossProd = sp.vcrss(uvPlane, uvNorm)
        err2 = sp.vdot(crossProd,crossProd)
        iwct = len(numpy.where(surfpt != lastSurfpt)[0])

      self.surfpt = surfpt
      self.height = sp.vdot(uvPlane,sp.vsub(self.surfpt,pPlane))
      self.errHgt = self.height - sp.vdot(uvNorm,sp.vsub(pPlane,self.surfpt))
      self.errAng = numpy.sqrt(err2)

      if err2 >eps2:
        self.status = 'Failed to converg'
        raise Exception('FAILED TO CONVERGE')

      self.converged = True
      self.status = 'Converged successfully'
      self.msg = 'OK'

    except:
      self.status = 'EllSurfTangent:  Failed'
      import traceback
      self.status = 'Failed to converg'
      self.msg = traceback.format_exc()

    ###print(vars(self))
    return self

"""

if n_elements(uvPlaneNormal) ne 3L then cspice_vhat,-[1,5d2,1],uvPlaneNormal
if n_elements(cPlane) ne 1L then cPlane = 500d0
if n_elements(abc) ne 3L then abc=[14d3,132,66]
btcdirpln_nvc2pl, uvPlaneNormal, cPlane, plane
btcdirpln_nvc2pl, uvPlaneNormal, -cPlane, plane2
eps=0d0
#eps=1d-5
rtn = ellsurftangent( abc, plane, eps=eps, debug=debug)
rtn2 = ellsurftangent( abc, plane2, eps=eps, debug=debug)
help,rtn
help,/st,rtn
print,rtn.surfpt,rtn.abc
help,rtn2
help,/st,rtn2
print,rtn2.surfpt,rtn2.abc

end
"""
