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
  fmt = ffmt.join('{0 }'.split())
  return '\n'.join(
         [prefix + '[' + ', '.join([fmt.format(v) for v in row]) + ']'
          for row in matrix
         ]
         )

########################################################################
makevec = lambda seq: sp.vequ(list(map(float,seq[:3])))


########################################################################
def pbe2mdl( abcArg
           , r
           , fn
           , iWhichPbeArrArg=0
           , downtrackScaleArg=1.0
           , uptrackScaleArg=1.0
           , resolutionArg=5.0
           , mtxJ2k2UncertArg=sp.ident()
           , userCommentsArg=None
           , nameLabelArg=None
           , userLabelArg=None
           , wrlArg=False
           , debugArg=None
           ):

  dpr,halfpi,onepi,twopi = sp.dpr(),sp.halfpi(),sp.pi(),sp.twopi()
  res = float(resolutionArg)
  nlon0 = round(360e0 / res)
  nlat0 = round(180e0 / res) + 1

  sigmaMultiple = baseSigmas = sUNKNOWN = 'UNKNOWN'
  satRTN = pluRTN = pluRTNSrc = xKinetxCall = sUNKNOWN
  pbecalcSuccess = False

  abc = None

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

    try: xKinetxCall = abcArg.kinetx.KinetxCall
    except: pass

    try: pbecalcSuccess = abcArg.success
    except: pass

    try:
      assert pbecalcSuccess

      ### Get plumped ABC and J2k-to-ABC xform from PBECALCS.PRO return value

      (abc,sigmaMultiple,baseSigmas,satRTN,pluRTN,pluRTNsrc,mtxJ2uLcl
      ,) = (makevec(abcArg.pbeArr[iWhichPbeArr].plumpedABC[:3])
           ,str(abcArg.sigmaMultiple).strip()
           ,v2s(abcArg.baseSigmas)
           ,v2s(abcArg.satEphemUncert)
           ,v2s(abcArg.pluEphemUncert)
           ,str(abcArg.pluEphemUncertSource).strip()
           ,abcArg.kinetx.mJ2u
           ,)
    except:
      try:
        abc = makevec(abcArg[iWhichPbeArr].plumpedABC[:3])
      except:
        assert False,'Unable to get ABC values from argument [{abcArg}]'

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

  latsArr = (numpy.arange(nlat0,dtype=numpy.float64) * sp.pi() / (nlat0-1)) - halfpi

  nlonsArr = numpy.round(nlon0 * numpy.cos(latsArr)).astype(numpy.int)
  iw = numpy.where(nlonsArr < 6)
  if iw[0].shape: nlonsArr[iw]=6
  nLon = nlonsArr[0] = nlonsArr[-1] = 1
  nextLonIdx = numpy.cumsum(nlonsArr)
  firstLonIdx = numpy.insert(nextLonIdx[:-1],0,0)

  #print(nlonsArr)
  #print(firstLonIdx)
  #print(nextLonIdx)

  ######################################################################
  ### vertices:  index, lon, lat, normal to ellipse => surface point +R => XYZ

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
  ### connection list

  nC = (2 * nV) - 4

  class CONN:
    def __init__(self,idx): self.idx = idx
    def set_verts(self,verts): self.verts = verts

  connArr = [CONN(idx) for idx in range(nC)]

  botLat = 0
  topLat = botLat + 1
  nLonBot = nlonsArr[botLat]
  nLonTop = nlonsArr[topLat]
  offsetBot = offsetTop = sumBot = sumTop = 0

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

    ### restart the loop

  """
  topLat = 1
  while topLat < nlat0:
    botLat = topLat - 1
    offsetBot = offsetTop = sumBot = sumTop = 0
    nLonBot = nlonsArr[botLat]
    nLonTop = nlonsArr[topLat]
    while (offsetBot < nLonBot and nLonBot > 1
          ) or (offsetTop < nLonTop and nLonTop > 1):
      nextSumBot = sumBot + nLonTop
      nextSumTop = sumTop + nLonBot
      incBot = (offsetTop == nLonTop or nLonTop == 1
               ) or (nextSumBot < nextSumTop
               ) or (nextSumBot == nextSumTop and sumBot < sumTop
               )
      triple = [firstLonIdx[botLat] + (offsetBot % nLonBot)
               ,firstLonIdx[topLat] + (offsetTop % nLonTop)
               ]
      if incBot:
        offsetBot += 1
        insertIdx = firstLonIdx[botLat] + (offsetBot % nLonBot)
        sumBot = nextSumBot
      else:
        offsetTop += 1
        insertIdx = firstLonIdx[topLat] + (offsetTop % nLonTop)
        sumTop = nextSumTop
      triple.insert(1,insertIdx)
    topLat += 1
  """


  ######################################################################
  ### Various constants

  maxR2 = max([sp.vdot(v.xyz,v.xyz) for v in vArr])
  maxR = numpy.sqrt(maxR2)
  boundRad = str(maxR * 1040e0)  ### BoundRadius:  scale by 1000, add 4%

  ######################################################################
  ### Description / comments

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
  print(descr)

  with open(fn,'w') as fout:
    if wrlArg: pass

    vts = [v2s(c.verts,' ') for c in connArr]
    dataLbl = getDataLbl(fn,nNameLbl,userLabelArg,vArr)

    rot = sp.m2eul(mtxJ2uLcl, 1, 2, 3)  ### ax, ay, az
    rotStr = v2s(sp.vscl(dpr,rot))

    shin = 76e0

    axisStr = [s.strip().split('/') for s in """
Xaxis/255000000/0   0  0
Yaxis/000255000/0   0 90
Zaxis/000000255/0 -90  0
""".strip().split('\n')]

    axes = ['']
    for i in range(3):
      axes.extend("""
Component {axisStr[i][0]
  Translucency 0.67
  FaceColor %{axisStr[i][1]
  Scale {r+abc[i]} 1 1
  Rotate {axisStr[i][2]}
  Refer
    Component BaseXaxis
  EndRefer
EndComponent

""".split('\n'))

"""

    printf,lun, f='(a)'
    , '# STK/VO 3D Model file'
    , ''
    , '### ' + [ 'ModelDesc', '  ' + descr , 'EndModelDesc']
    , ''
    , 'Component BaseXaxis'
    , '  Revolve'
    , '    StartAngle 0.0'
    , '    EndAngle 360.0'
    , '    NumRevolve 10'
    , '    NumVerts 2'
    , '    Data'
    , '      0 10 0'
    , '      1.03 10 0'
    , '  EndRevolve'
    , 'EndComponent'
    , axes
    , ''
    , 'Component XYZaxes'
    , '  Refer'
    , '    Component Xaxis'
    , '  EndRefer'
    , '  Refer'
    , '    Component Yaxis'
    , '  EndRefer'
    , '  Refer'
    , '    Component Zaxis'
    , '  EndRefer'
    , 'EndComponent'
    , ''
    , 'Component PBE'
    , '  PolygonMesh'
    , '    FaceColor %255255255'
    , '    SmoothShading No'
    , '    BackfaceCullable Yes'
    , '    Translucency 0.5'
    , '    Specularity ' + strtrim(shin/128e0,2)
    , '    Shininess ' + strtrim(shin,2)
    , (nLbl eq 0L) ? '' : '    ' +
      [ 'Texture'
      , '  RGB ' + textureFn
      , '  Parm ' + ( 1?['AA','Transp']:['AA','Transp','Mipmap','ClampS','ClampT'])
      , 'EndTexture'
      ]
    , '    NumVerts ' + strtrim(n_elements(vArr),2)
    , '    ' + dataLbl
    , '      ' + pts
    , '    NumPolys ' + strtrim(n_elements(cArr),2)
    , '    Polys'
    , '      3 ' + vts
    , '  EndPolygonMesh'
    , 'EndComponent'
    , ''
    , 'Component PBE_Rotated'
    , '  Rotate ' + rotStr
    , '  Refer'
    , '    Component PBE'
    , '  EndRefer'
    , '  Refer'
    , '    Component XYZaxes'
    , '  EndRefer'
    , 'EndComponent'
    , ''
    , 'Component PBE_ROOT'
    , '  Root'
    , '  UniformScale 1000.0' $       ### Convert km to m
    , '  Refer'
    , '    Component PBE_Rotated'
    , '  EndRefer'
    , '  #Refer' $                    ### J2000 axes removed
    , '  #  Component XYZaxes'
    , '  #EndRefer'
    , 'EndComponent'
    , 'BoundRadius ' + boundRad

  endif else begin

    ####################################################################
    ### - WRL

    pts = strjoin( strtrim(vArr.xyz * 5e0 / maxR ,2), ' ')
    vts = strjoin( strtrim(cArr.verts,2), ' ')
    printf,lun, f='(a)'
    , '#VRML V2.0 utf8'
    , '### ' + [ '', 'Model description:', '  ' + descr, '']
    , 'Shape {'
    , '  geometry IndexedFaceSet {'
    , '    ccw TRUE  ### East longitude'
    , '    coord Coordinate {'
    , '      point [ ' + pts[0]
    , '            , ' + pts[1:*]
    , '            ]'
    , '    }'
    , '    coordIndex [ ' + vts[0] + ' -1'
    , '               , ' + vts[1:*] + ' -1'
    , '               ]'
    , '  }'
    , '  appearance Appearance {'
    , '    material Material { '
    , '      diffuseColor .25 .25 .25 '
    , '      emissiveColor .25 .25 .25 '
    , '    }'
    , '  }'
    , '}'
  endelse
  free_lun,lun

  return
  message,'OK'
end
"""

"""
    ########################################################################
    ########################################################################
    ### Test code

      cspice_furnsh,'v_od059b.tm'
      k=kinetxcurrent(target=999)

      ###sigABC=[1378e0,66e0,33e0]
      sigABC=kget(k,/KnSig)
      sigABC2=sigABC * 2e0

      pbe2mdl, sigABC2, 91.0e0, 'pbetesthydra_2sigma.mdl', res=5e0, mtxJ2k=k.mJ2u
      pbe2mdl, sigABC2, 91.0e0, 'pbetesthydra_2sigma.wrl', res=5e0, /wrl

      pbe2mdl, sigABC2, 76.8e0, 'pbetestnix_2sigma.mdl', res=5e0, mtxJ2k=k.mJ2u
      pbe2mdl, sigABC2, 76.8e0, 'pbetestnix_2sigma.wrl', res=5e0, /wrl

      sigTmp=[100e0,35,30]
      pbe2mdl, sigTmp, 91.0e0, 'pbetestasymm.mdl', up=2.5e0, down=1.5e0, mtxJ2k=k.mJ2u
      pbe2mdl, sigTmp, 91.0e0, 'pbetestasymm.wrl', up=2.5e0, down=1.5e0, /wrl

      hRTN=[ 64e0, 1512, 102]    ### Hydra
      dEts=[-100e0,200,500]
      pbecalcs,'903',dEts=dEts,satE=hRTN, debug=debug, pbeStr=pbe903

      pbe2mdl, pbe903, 91e0, 'pbetesthydra_pbecalc_asymm.mdl', up=2.5, down=1.5
      pbe2mdl, pbe903, 91e0, 'pbetesthydra_pbecalc_asymm.wrl', up=2.5, down=1.5, /wrl

      pbe2mdl, pbe903, 91e0, 'pbetest_texture00.mdl'
             , up=2.5, down=1.5
             , userLabel='Test User Label 00 '+['0','1']

      pbe2mdl, pbe903, 91e0, 'pbetest_texture01.mdl'
             , up=2.5, down=1.5
             , /nameLabel

      pbe2mdl, pbe903, 91e0, 'pbetest_texture02.mdl'
             , up=2.5, down=1.5
             , /nameLabel
             , userLabel='Test User Label 02 '+['0','1']
end

"""
def getDataLbl(fn,nNameLbl,userLabelArg,vArr):
  return 'Data'
  """

  ######################################################################
  ### Choose WRL or MDL

  if keyword_set(wrlArg) eq 0b then begin

    ####################################################################
    ### - MDL

    pts = strjoin( strtrim(vArr.xyz,2), ' ')

    if nLbl gt 0L then begin
      textureFn = fn
      posSlash=strpos(textureFn,'/',/reverse_search)
      posDot=strpos(textureFn,'.',/reverse_search)
      if posDot gt (posSlash+1L) then textureFn = strmid(textureFn,0,posDot)
      textureFn = textureFn + '.png'

      lclLbls= nNameLbl ? strtrim(['',fn],2) : ['']
      if nUserLbl gt 0L then lclLbls=[lclLbls,strtrim(userLabelArg[*],2)]
      lclLbls = lclLbls[1:*]

      u=transpose(vArr.lon/twopi)
      v=transpose((vArr.lat+!dpi/2)/!dpi)
      uv=strjoin(strtrim([u,v],2),' ')
      pts = pts + ' ' + uv
      dataLbl = 'DataTx'
      write_png,textureFn, pbetext(lclLbls)
    endif else begin
      dataLbl = 'Data'
    endelse
"""
