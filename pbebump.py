import sys
import numpy
import spiceypy as sp
from ellsurftangent import DIRECTEDPLANE,EllSurfTangent

try: dpr
except:
  dpr = sp.dpr()
  v0 = sp.vpack(0,0,0)
########################################################################
### Find plane passing through observer position (argument obsPos[3])
###   that passes at a known distance (argument targRadius) to the
###   nearest point on an ellipse.  The orientation of the plane is such
###   that its normal is in a plane with obsPos[3] and a ray (argument
###   scanRay[3]) passing through the center of the ellipse.  The dot
###   product of scanRay[3] and the plane normal will be negative, so
###   the scanRay vector will pass through the plane.  The A, B, & C
###   semi-axes of the ellipse are in argument abc[3].
###
###   Output structure members are 1) the angle between the plane and a
###   vector from the observer to the ellipse center.
########################################################################
class PBEBUMP:
  def __init__(self
,scanRay=None
,obsPos=None
,abc=None
,targRadius=0.0
,epsilon=1e-7
,maxIter=20
,sample=False
):
    if sample:
      self.scanRay = sp.vequ(v0)
      self.scanAxis = sp.vequ(v0)
      self.planeNorm = sp.vequ(v0)
      self.aimPt = sp.vequ(v0)
      self.nom2Tangent_Deg = self.rangeAlongScanRay = 0e0
      self.iter = 0
      return

    ####################################################################
    ### Convention:  shown as projections into nominal FOV plane normal
    ###   to argument obsPos[3]
    ###
    ### - Character @  at center of ellipse represents argument obsPos[3],
    ###     nominal vector from center of ellipse to observer, coming out
    ###     of the screen
    ###
    ### - Argument scanRay[3] is line of scan# it need not be in FOV
    ###
    ### - VscanRayN are the projection of argument scanRay[3] into the
    ###     FOV plane normal to obsPos[3].  Ellipse need not be symmetrical
    ###     about VscanRayN.
    ###
    ### - planeN with normal plNormN will be the plane associated with VscanRayN
    ###
    ### - obsPos[3] cross VscanRayN will be VscanAxisN# plane will rotate
    ###     around VscanAxisN to satisfy targRadius constraint.  VscanAxisN
    ###     will not change with plane rotating out from center
    ###
    ### - aimPt[3] will be aimpoint through which plane passes, marked by +
    ###
    ### - VscanAxisN cross [obsPos[3]-aimpt] will be plNorm1
    ###
    ###                                   ^
    ###                                   |VscanAxis1
    ### |                                 |                                 |
    ### |                               __|__                         plane1|
    ### |                            _--  |  --_                            |
    ### |                           /     |     \                           |
    ### +-->plNorm2  VscanRay2<----|------@------|---->VscanRay1  plNorm1<--+
    ### |                           \_    |    _/                           |
    ### |                             --__|__--                             |
    ### |plane2                           |                                 |
    ### |                       VscanAxis2|
    ###                                   |
    ###                                   v
    ####################################################################

    self.scanRay = scanRay
    self.scanAxis = sp.ucrss(obsPos, self.scanRay)   ### Get scan axis

    itercount = 0
    aimPt = sp.vpack(0,0,0)
    maxIter = 20
    while itercount < maxIter:
      itercount += 1
      planeNorm = sp.ucrss(sp.vsub(obsPos,aimPt), self.scanAxis)
      plane = DIRECTEDPLANE.nvp2pl(planeNorm, aimPt)
      tangInfo = EllSurfTangent().calculate(abc, plane
                                           ,epsilon=epsilon
                                           ,maxIter=maxIter
                                           )
      try:
        if abs(targRadius - tangInfo.height) < 1e-10: break
      except:
        import pprint
        print(tangInfo)
        pprint.pprint(vars(tangInfo))
        raise

      aimPt = sp.vsub(tangInfo.surfpt, sp.vscl(targRadius,planeNorm))

    self.iter = itercount
    self.planeNorm = planeNorm
    self.aimPt = aimPt

    self.nom2Tangent_Deg = dpr * sp.vsep(obsPos, sp.vsub(obsPos,aimPt))

    ### Intersection of scanRay with plane
    nxpts, scanRayPt = sp.inrypl(v0, scanRay, plane.csppl())
    self.rangeAlongScanRay = sp.vnorm(scanRayPt)

    return

if "__main__" == __name__ and sys.argv[1:] and 'test' == sys.argv[1]:
  obsPos = sp.vpack(3,4,0)
  scanRay = sp.vpack(0,1,0)
  abc = sp.vequ([1.25]*3)
  targRadius = abc[0]
  import pprint
  pprint.pprint(dict(pbPlusY=vars(PBEBUMP(scanRay,obsPos,abc,targRadius))
                    ,pbMinusY=vars(PBEBUMP(sp.vsub(v0,scanRay),obsPos,abc,targRadius))
                    ))
"""
### function pbebump, scanRay, obsPos, abcArg, targRadius

### special case to test code:  sphere i.e. A=B=C

obsPos = [3.,4.,0]         ### Obs position in XY plane
scanRay = [0d0,1,0]        ### scan ray in XY plane along Y axis
abc = [1.25d0,1.25,1.25]   ### Sphere of radius 1.25
targRadius = 1.25d0        ### Target radius is 1.25

### plane should pass 1.25 + 1.25 = 2.5km from center, observer at 5km from
### center, angle should be arcsin(2.5/5) = 30 degrees
###
### pbPlusY range should be   2.718528 = 4 - 3 * tan(arctan(4/3)-PI/6))
### pbMinusY range should be 20.900346 = 3 * tan(arctan(4/3)+PI/6)) - 4

pbPlusY = pbebump( scanRay, obsPos, abc, targRadius, tangInfo=tiPlus)
pbMinusY = pbebump( -scanRay, obsPos, abc, targRadius, tangInfo=tiMinus)
help,/st,pbPlusY,pbMinusY

end
"""
