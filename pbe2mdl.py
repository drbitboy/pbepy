#!/usr/bin/env python3

import os
import sys
import numpy
import collections
import spiceypy as sp
from vertex import VERTEX

########################################################################
def v2s(vector,sep=', '):
  """Vector to string"""
  return '[' + sep.join(map(str,vector)) + ']'

def m2ss(matrix,sep='\n',ffmt=':19.15f',prefix='     '):
  """Matrix to multi-line string"""
  fmt = ffmt.join('{0 }'.split())
  return '\n'.join(
         [prefix + '[' + ', '.join([fmt.format(v) for v in row]) + ']'
          for row in matrix
         ]
         )

########################################################################
### create a vector from three floats
makevec = lambda seq: sp.vequ(list(map(float,seq[:3])))


########################################################################
### Build a data label
def getDataLbl(fn,nLbl,nNameLbl,userLabelArg,vArr):
  textureFn = fn
  return nLbl and 'DataTx' or 'Data',nLbl and f"""
    Texture
      RGB {textureFn}
      Parm AA
      Parm Transp
    EndTexture
""".lstrip('\r\n') or ''


########################################################################
def pbe2mdl( abcArg                       ### Ellipsoid semi-axes data*
           , r                            ### Radius bump
           , fn                           ### output MDL filename
           , iWhichPbeArrArg=0            ### Index into abcArg.pbeArr
           , downtrackScaleArg=1.0        ### Downtrack scaling
           , uptrackScaleArg=1.0          ### Uptrack scaling
           , resolutionArg=5.0            ### Nom. MDL resolution, deg
           , mtxJ2k2UncertArg=sp.ident()  ### Ellipsoid orientation, J2k
           , userCommentsArg=None         ### User comments for MDL
           , nameLabelArg=False           ### Write fn to MDL comments?
           , userLabelArg=None
           , wrlArg=False
           , debugArg=None
           ):
  """Write STK MDL of Plumped Bumped Ellipse (PBE)"""

  ### * May be an instance of the PBESTRUCT class; see pbecalcs.py

  ### Misc constants
  dpr,halfpi,onepi,twopi = sp.dpr(),sp.halfpi(),sp.pi(),sp.twopi()
  res = float(resolutionArg)
  nlon0 = round(360e0 / res)
  nlat0 = round(180e0 / res) + 1

  ### Initialize provenance strings
  sigmaMultiple = baseSigmas = sUNKNOWN = 'UNKNOWN'
  satRTN = pluRTN = pluRTNSrc = xKinetxCall = sUNKNOWN
  pbecalcSuccess = False

  try:
    ### 3-element ABC from input argument sequence
    abc = makevec(abcArg)
    noArgYet = False

  except:

    ####################################################################
    ### Test for output of pbecalcs i.e object pbeStr or pbeStr.pbeArr
    ###
    ### - If .pbeArr is present
    ###   - Check .SUCCESS
    ###   - ABC[3] from .pbeArr[I].plumpedABC
    ###   - mtxJ2u[3,3] from .Kinetx.mJ2u (if present)
    ###
    ### - else check for .plumpedABC
    ###   - use abcArg[I].plumpedABC if .PBEARR not present
    ###
    ### - use iWhichPbeArrArg to select which pbeArr to use; default = 0

    iWhichPbeArr = int(iWhichPbeArrArg)

    ### Get provenance string from KinetX object
    try: xKinetxCall = abcArg.kinetx.KinetxCall
    except: pass

    ### Get success string from KinetX object
    try: pbecalcSuccess = abcArg.success
    except: pass

    try:
      assert pbecalcSuccess

      ### Get plumped ABC and J2k-to-ABC xform from PBECALCS.PRO return value

      (abc,sigmaMultiple,baseSigmas,satRTN,pluRTN,pluRTNsrc,mtxJ2uLcl
      ,) = (makevec(abcArg.pbeArr[0][iWhichPbeArr].plumpedABC[:3])
           ,str(abcArg.sigmaMultiple).strip()
           ,v2s(abcArg.baseSigmas)
           ,v2s(abcArg.satEphemUncert)
           ,v2s(abcArg.pluEphemUncert)
           ,str(abcArg.pluEphemUncertSource).strip()
           ,abcArg.kinetx.mJ2u
           ,)
    except:
      try:
        abc = makevec(abcArg[0][iWhichPbeArr].plumpedABC[:3])
      except:
        assert False,f'Unable to get ABC values from argument [{abcArg}]'

  ######################################################################
  ### Override any previous mtxJ2u with argument:

  if not (mtxJ2k2UncertArg is None): mtxJ2uLcl = mtxJ2k2UncertArg

  nNameLbl = nameLabelArg and True or False
  nUserLbl = userLabelArg and len(userLabelArg) or 0
  nLbl = nUserLbl + (nNameLbl and 1 or 0)

  ######################################################################
  ### Downtrack (X>=0) & Uptrack (X<0) scaling

  dtSc = float(downtrackScaleArg)
  utSc = float(uptrackScaleArg)

  ######################################################################
  ### Build parallels array of latitudes

  latsArr = (numpy.arange(nlat0,dtype=numpy.float64) * sp.pi() / (nlat0-1)) - halfpi

  ### Build array with number of longitude vertices per parallel

  nlonsArr = numpy.round(nlon0 * numpy.cos(latsArr)).astype(numpy.int32)

  ### Ensure there are at least six vertices per parallel, except at the
  ### two poles, which have one each

  iw = numpy.where(nlonsArr < 6)
  if iw[0].shape: nlonsArr[iw]=6
  nLon = nlonsArr[0] = nlonsArr[-1] = 1

  ### Build the starting and ending vertices' indices at each parallel

  nextLonIdx = numpy.cumsum(nlonsArr)
  firstLonIdx = numpy.insert(nextLonIdx[:-1],0,0)

  ######################################################################
  ### Calculate vertices' Cartesian coordinates

  nV = nextLonIdx[-1]
  vArr = [VERTEX(idx,abc,dtSc,utSc,r) for idx in range(nV)]
  iLat,xLat = 0,-halfpi
  for v in vArr:
    if v.idx == nextLonIdx[iLat]:
      iLat += 1
      nLon = nlonsArr[iLat]
      xLat = (iLat * onepi / (nlat0-1)) - halfpi
    xLon = (v.idx - firstLonIdx[iLat]) * twopi / nLon
    v.set_normal(xLat,xLon).calculate()

  ######################################################################
  ### Build polygon (triangle) list; each polygon ha 3 vertex indices

  nC = (2 * nV) - 4

  class CONN:
    def __init__(self,idx): self.idx = idx
    def set_verts(self,verts): self.verts = verts

  connArr = [CONN(idx) for idx in range(nC)]

  ### Initialize quantities for loop

  botLat = 0
  topLat = botLat + 1
  nLonBot = nlonsArr[botLat]
  nLonTop = nlonsArr[topLat]
  offsetBot = offsetTop = sumBot = sumTop = 0

  ### Loop over polygons

  for conn in connArr:

    ### Calculate
    ###   ((offsetBot+1)/nLonBot) * (nLonBot*nLonTop)
    ### and
    ###   ((offsetTop+1)/nLonTop) * (nLonBot*nLonTop)
    ### to select which parallel, bottom or top, should be advanced on
    ### this pass; if unequal, the parallel with the lower longitude,
    ### the proxy for which is the ratio (offsetXxx/nLonXxx), will be
    ### advanced
    nextSumBot = sumBot + nLonTop
    nextSumTop = sumTop + nLonBot

    ### Advance by one longitude increment along the bottom parallel if
    ### - EITHER the top parallel's offset has come back to its start
    ### - OR there is one longitude in the top parallel (north pole)
    ### - OR
    ###   - BOTH the bottom parallel is not the south pole
    ###   - AND
    ###     - EITHER the next bottom longitude would be < the top's
    ###     - OR
    ###       - BOTH the next bottom and top longitudes would be equal
    ###       - AND the bottom increment in longitude would be smaller
    incBot = (offsetTop == nLonTop
             ) or (nLonTop == 1
             ) or ((nLonBot > 1)
                   and ((nextSumBot < nextSumTop)
                        or
                        (nextSumBot == nextSumTop and sumBot < sumTop)
                       )
                  )

    ### The first and last vertex are the bottom and top current offsets
    triple = [firstLonIdx[botLat] + (offsetBot % nLonBot)
             ,firstLonIdx[topLat] + (offsetTop % nLonTop)
             ]

    ### The middle vertex will be the one that is incremented
    if incBot:
      offsetBot += 1
      middleIdx = firstLonIdx[botLat] + (offsetBot % nLonBot)
      sumBot = nextSumBot
    else:
      offsetTop += 1
      middleIdx = firstLonIdx[topLat] + (offsetTop % nLonTop)
      sumTop = nextSumTop

    ### Insert the middle vertex, write it to the connection
    triple.insert(1,middleIdx)
    #print((botLat,topLat,offsetBot,offsetTop,triple,))
    conn.set_verts(triple)

    ### Continue with the current parallel rows if they are incomplete
    if (offsetBot < nLonBot and nLonBot > 1) or (offsetTop < nLonTop):
      continue

    ### Increment the parallel rows, reset the rows' parameters
    botLat = topLat
    topLat += 1
    nLonBot = nlonsArr[botLat]
    nLonTop = nlonsArr[topLat]
    offsetBot = offsetTop = sumBot = sumTop = 0

    ### End of loop

  ######################################################################
  ### Various dimensional constants

  maxR2 = max([sp.vdot(v.xyz,v.xyz) for v in vArr])
  maxR = numpy.sqrt(maxR2)
  boundRad = str(maxR * 1040e0)  ### BoundRadius:  scale by 1000, add 4%

  ######################################################################
  ### Build PBE description and provenancew as MDL comments

  commaj = ','.join
  fromsrc = lambda src:f' from {src}' if src else ''

  def slist(lbl,seq,prefix='  '):
    if isinstance(seq,str): return slist(lbl,[seq])
    if isinstance(seq,collections.Sequence):
      return [f'{prefix}{lbl}:']+[f'{prefix}  {item}' for item in seq]
    return []
  def sslist(lbl,seq,prefix='  '):
    return '\n'.join(slist(lbl,seq,prefix=prefix))

  def list_kernels(prefix='  # - '):
    return '\n'.join([f'{prefix}{fil} type={typ}{fromsrc(src)}'
                      for fil,typ,src,handle in
                      [sp.kdata(i,'all')
                       for i in range(sp.ktotal('all'))
                      ]])

  ### Description and provenance string

  descr = f"""
ModelDesc
  Plumped Bumped Ellipse:
    A,B,C = {commaj(map(str,abc))}
    - includes nSigma scaling = {sigmaMultiple}
    - includes base (ToF,Bnorm,Bmag) sigmas = {baseSigmas}
    - includes Target ephem (RTN) sigmas = {satRTN}
    - includes Pluto ephem (RTN) sigmas = {pluRTN}
      - Source = {pluRTNSrc}
    - includes kinetxNNN = {xKinetxCall}
    Scaling downtrack,uptrack = {dtSc},{utSc}
    Radius bump = {r}
    Xform matrix, J2k to ABC:
{m2ss(mtxJ2uLcl)}
  #  Spice Kernels:
{list_kernels()}
{sslist('User-supplied comments',userCommentsArg)}
{sslist('Filename label',fn if nNameLbl else None)}
{sslist('User-supplied Label',userLabelArg if nUserLbl else None)}
EndModelDesc""".replace('\n\n\n','\n'  ### Remove empty lines ...
).replace('\n\n','\n'                  ### ... and again
).replace('\n','\n### '                ### Make each line a comment
).lstrip()                             ### Remove initial newline
  ###print(descr)

  with open(fn,'w') as fout:
    if wrlArg: pass

    ### Build vertex and vertices list strings
    vtxs = [v2s(vtx.xyz,sep=" ").strip("[]") for vtx in vArr]
    conns = [v2s(conn.verts,sep=' ').strip('[]') for conn in connArr]
    dataLbl,texture = getDataLbl(fn,nLbl,nNameLbl,userLabelArg,vArr)

    rot = sp.m2eul(mtxJ2uLcl, 1, 2, 3)  ### ax, ay, az
    rotStr = v2s(sp.vscl(dpr,rot),sep=' ').strip('[]')

    shin = 76e0

    axisStr = [s.strip().split('/') for s in """
Xaxis/255000000/0   0  0
Yaxis/000255000/0   0 90
Zaxis/000000255/0 -90  0
""".strip().split('\n')]

    axes = ''
    for i in range(3):
      axes += f"""
Component {axisStr[i][0]}
  Translucency 0.67
  FaceColor %{axisStr[i][1]}
  Scale {r+abc[i]} 1 1
  Rotate {axisStr[i][2]}
  Refer
    Component BaseXaxis
  EndRefer
EndComponent
"""

    fout.write(f"""
# STK/VO 3D Model file

{descr}

Component BaseXaxis
  Revolve
    StartAngle 0.0
    EndAngle 360.0
    NumRevolve 10
    NumVerts 2
    Data
      0 10 0
      1.03 10 0
  EndRevolve
EndComponent

{axes}

Component XYZaxes
  Refer
    Component Xaxis
  EndRefer
  Refer
    Component Yaxis
  EndRefer
  Refer
    Component Zaxis
  EndRefer
EndComponent

Component PBE
  PolygonMesh
    FaceColor %255255255
    SmoothShading No
    BackfaceCullable Yes
    Translucency 0.5
    Specularity {shin/128e0}
    Shininess {shin}

{texture}

    NumVerts {nV}
    {dataLbl}
""".lstrip('\r\n'))

    for vtx in vtxs: fout.write(f'      {vtx}\n')

    fout.write(f"""
    NumPolys {nC}
    Polys
""".lstrip('\r\n'))

    for conn in conns: fout.write(f'      3 {conn}\n')

    fout.write(f"""
  EndPolygonMesh
EndComponent

Component PBE_Rotated
  Rotate {rotStr}
  Refer
    Component PBE
  EndRefer
  Refer
    Component XYZaxes
  EndRefer
EndComponent

Component PBE_ROOT
  Root
  UniformScale 1000.0
  Refer
    Component PBE_Rotated
  EndRefer
  #Refer
  #  Component XYZaxes
  #EndRefer
EndComponent
BoundRadius {boundRad}
""".lstrip('\r\n'))

########################################################################
if "__main__" == __name__:
  try:
    L = len(sys.argv)
    target,observer,kerns,mdlpath = sys.argv[1:5]
    kerns = kerns.split(',')
    radi = L > 5 and float(sys.argv[5]) or 0.0
    sigm = L > 8 and list(map(float,sys.argv[6:9])) or None
    satE = L > 11 and list(map(float,sys.argv[9:12])) or 0.0
    pStE = L > 14 and list(map(float,sys.argv[12:15])) or None

    """Instantiate  PBE structure"""
    import pprint
    pprint.pprint(dict(pbe2mdl_main=locals()))
    from pbecalcs import PBESTRUCT
    pbec = PBESTRUCT(targArg=target,obsArg=observer,kern=kerns
                    ,radi=radi,sigm=sigm,satE=satE,pStE=pStE
                    )
    """Write MDL file from that PBE structure"""
    pbe2mdl(pbec,radi,mdlpath,mtxJ2k2UncertArg=pbec.mtx_j2k2Uncert)
  except:
    if 'DEBUG' in os.environ:
      import traceback as tb
      tb.print_exc()
    print("""
Usage:

  pbe2mdl.py target observer kernel[,kernel[,...]] MDLpath \\
             [radius_bump \\
             [scTOF,km      scBnorm,km    scBmag,km \\
             [targetR,km    targetT,km    targetN,km   \\
             [primaryR,km   primaryT,km   primaryN,km]]]]
""")
