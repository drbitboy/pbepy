import spiceypy as sp
from ellsurftangent import DIRECTEDPLANE,EllSurfTangent


########################################################################
class VERTEX:
  """.xyz is a point at an altitude above an elllipsoid surface point
.normal is the unit normal vector at that surface point
.lon is normal's East longitude (RA)
.lat is normal's latitude (DEC)
- .lat and .lon refer to surface normal in radians

"""
  def __init__(self,idx,abc,dtSc,utSc,altitude):
    """idx is the index number of this vertex
abc has the three radii of the ellipsoid
dtSc is the downtrack scaling of the ellipsoid
utSc is the uptrac scaling of the ellipsoid
r is the altitude above the surface point

"""
    self.idx,self.abc,self.ready = idx,abc,False
    self.dtSc,self.utSc,self.altitude = map(float,(dtSc,utSc,altitude,))

  def set_normal(self,xLat,xLon):
    """xLat is the latitude (DEC) of the surface normal
xLon is the East longitude of the surface normal

"""
    self.lat,self.lon = xLat,xLon
    self.norm = sp.radrec(1e0,self.lon,self.lat)
    return self

  def calculate(self):
    """Calculate surface point based on normal, abc, scaling
- Use a plane passing through the center of the ellipsoid, with the
  plane's normal opposite the normal defined by .lat/.lon in __init__
"""
    dirplane = DIRECTEDPLANE.nvc2pl(-self.norm, 0e0)
    scale = self.dtSc if self.norm[0] >= 0e0 else self.utSc
    tangInfo = EllSurfTangent().calculate(self.abc * scale, dirplane)
    self.xyz = tangInfo.surfpt + (self.norm * self.altitude)
    self.ready = True
    ##print(vars(self))
