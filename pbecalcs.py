import math
import numpy
import collections
import spiceypy as sp
from pbebump import PBEBUMP
from kinetxcurrent import kinetxcurrent
from findsatephuncabc import FINDSATEPHUNCABC

spd,dpr = sp.spd(),sp.dpr()

########################################################################
### Plumped, Bumped Ellipsoid calculations
###
### Conventions:
###
###   sigm[3]:  A, B, C of ellipse; A is more or less TOF
###
### mEig[3,3]:  [ [ Xbpl, Ybpl, Zbpl ] $  of sigm[0] PBE axis
###             , [ Xbpl, Ybpl, Zbpl ] $  of sigm[1] PBE axis
###             , [ Xbpl, Ybpl, Zbpl ] ]  of sigm[2] PBE axis
###
class PBESTRUCT:

  def __init__(self
, targArg=2486958### STRING, default = '2486958' (MU69)
, obsArg=-98     ### STRING, default = '-98' (New Horizons)
, kern=[]        ### STRING[], SPICE kernels to load and unload
, ukoe=False     ### /ukoe = Unload SPICE Kernels On Exit
, radi=None      ### DOUBLE, Target radius for bumping, km
, UTCsArg=[]     ### STRING[], UTCs for which to calculate info
, dEtsArg=[]     ### DOUBLE, delta ET to add to UTCs for times at which to calc
, sigm=None      ### DOUBLE[3], ell sigmas, km; dflt=>use kinetxcurrent(target=targArg)
, mJ2u=None      ### DOUBLE[3,3], J2k to T],[LTOF]], dflt=> findvinf()
, nSig=2e0       ### DOUBLE, Sigma multiple to use, dflt=2
, satE=0e0       ### DOUBLE, Satellite ephemeris error, km, dflt=[0,0,0]
, pStE=None      ### DOUBLE, Pluto err, NIX & HYDRA only, km, dflt=findsatephuncabc
, mdlF=None      ### STRING, => MDL file name to which to write PBE for STK
, datF=None      ### STRING, => Flat ASCII filename to which to write info
, DelivSigmas=False ### Misc keywords
):
    unloadKernelsOnExit = ukoe and True or False
    self.success = False
    self.status = 'PBEcalcs SKIPPED'
    self.msg = ''

    try:
      self.status = 'PBEcalcs FAILED to load kernels'
      for k in kern: sp.furnsh(k)
      self.status = 'PBEcalcs FAILED'

      self.TargName = str(targArg)
      self.ObsName = str(obsArg)

      self.TargID = sp.bods2c(self.TargName)
      self.ObsID = sp.bods2c(self.ObsName)

      try   : self.TargName = sp.bodc2n(self.TargID,99)
      except: pass

      try   : self.ObsName = sp.bodc2n(self.ObsID,99)
      except: pass

      if radi is None:
        try   : t = numpy.max(sp.bodvcd(self.TargID,'RADII',4)[-1])
        except: t = 0e0
      else    : t = radi
      self.TargRadius = float(t)

      self.kinetx = kinetx = kinetxcurrent(target=targArg)

      ### KGET defaults to Knowledge Sigmas

      if sigm is None: self.baseSigmas = kinetx.kget(DeSig=DelivSigmas)
      else           : self.baseSigmas = list(map(float,sigm))

      if mJ2u is None: self.mtx_j2k2Uncert  = kinetx.mJ2u
      else           : self.mtx_j2k2Uncert = mJ2u

      self.sigmaMultiple = float(nSig)
      self.mdlFilename = mdlF
      self.dataFilename = datF

      try:
        L = len(satE)
        if L == 3: self.satEphemUncert = sp.vequ(*satE)
        else     : self.satEphemUncert = sp.vpack(0e0, satE[0], 0e0)
      except     : self.satEphemUncert = sp.vpack(0e0, satE, 0e0)

      self.pluEphemUncert = numpy.zeros(3,dtype=numpy.float64)
      self.pluEphemUncertSource = ''

      if isinstance(UTCsArg,str): self.UTCs = [UTCsArg]
      else:
        assert isinstance(UTCsArg,collections.Sequence)
        self.UTCs = list(map(str,UTCsArg))
      nUTCs = len(self.UTCs)
      if nUTCs == 0:
        self.UTCs = [sp.et2utc(kinetx.Tca.et,'ISOC',3)]
        nUTCs = 1

      try: self.dEts = [float(dEtsArg)]
      except: self.dEts = list(map(float,dEtsArg))
      ndEts = len(self.dEts)
      if ndEts == 0:
        self.dEts = [0e0]
        ndEts = 1

    except:
      import traceback
      self.message = traceback.format_exc()
      raise

    ### Create 2-D array PBEs (below) at various base UTCs and deltaETs
    lpbe = lambda iutc,idet: PBE(iutc,idet,self,pStE=pStE)

    arr =         [[lpbe(iUTC,idEt) for idEt in range(ndEts)]
                   for iUTC in range(nUTCs)
                  ]

    self.pbeArr = numpy.array(arr,dtype=object)
    """
    self.pbeArr = numpy.array(
                  [[lpbe(iUTC,idEt) for idEt in range(ndEts)]
                   for iUTC in range(nUTCs)
                  ])
    """

    self.success = True
    self.status = 'PBEcalcs OK'
    self.msg = 'OK'

    for k in kern: sp.unload(k)


########################################################################
class PBE:
  def __init__(self,iUTC,idEt,pbeStr,pStE=None):

    abc0Sq = pbeStr.baseSigmas*pbeStr.baseSigmas  ### Apply sigma multiplier later

    j2u = pbeStr.mtx_j2k2Uncert

    et0 = sp.utc2et(pbeStr.UTCs[iUTC])

    ### Time
    self.et = et = et0 + pbeStr.dEts[idEt]
    self.utc = utc = sp.et2utc(et, 'ISOC',3)
    self.ETminusTCA = et - pbeStr.kinetx.Tca.et

    ### Ephemeris

    stJ2k,ltim = sp.spkezr(pbeStr.TargName, et, 'J2000', 'LT', pbeStr.ObsName)
    self.nomVec_J2k = stJ2k[:3]
    r,ra,dec = sp.recrad(stJ2k[:3])
    self.nomRRaDec_J2k = [r, ra*dpr, dec*dpr]

    ### Perpendiculars to projection into FOV of Vinfinty & of Vtarget lines
    ### * stJ2k[:3] is vector to nominal target position
    ### * .mtx_j2k2Uncert[*,0] is Vinfinity, pointing downtrack
    ### * stJ2k[3:] is target velocity, Vtarg, pointing uptrack-ish
    ### * Orientation of FOV (in plane normal to nominal target position vec):
    ###   * Vinf or Vtarg is horizontal, Downtrack is Left
    ###   * VInfPerp and VTargPerp are Up

    vInf = pbeStr.mtx_j2k2Uncert[0].flatten()
    self.nomVinfPerp_J2k = vInfPerp = sp.ucrss(stJ2k[:3], vInf)
    vTargPerp = sp.ucrss(stJ2k[3:], stJ2k[:3])
    self.projVinfVtargRot_Deg = sp.vsep(vInfPerp,vTargPerp) * dpr

    ### Convert vectors to Uncertainty frame

    pTargU     = sp.mxv(j2u, stJ2k[:3])     ### Target position
    vTargU     = sp.mxv(j2u, stJ2k[3:])     ### Target velocity, uptrack
    vInfU      = sp.mxv(j2u, vInf)          ### Vinfinity, downtrack
    vInfPerpU  = sp.mxv(j2u, vInfPerp)
    vTargPerpU = sp.mxv(j2u, vTargPerp)

    ### Calculate offset btw Vinf & Vtarg at 200s in FOV in degrees

    tipTargDt = pTargU - (vTargU * 200e0)    ### 200s downtrack along vTarg
    tipTargUt = pTargU + (vTargU * 200e0)    ### 200s uptrack along vTarg
    l = sp.vnorm(vTargU) * 200e0
    tipInfDt = pTargU + (vInfU * l)          ### 200s downtrack along vInfU
    tipInfUt = pTargU - (vInfU * l)          ### 200s uptrack along vInfU

    tipDiffDt = tipInfDt - tipTargDt
    tipUiffUt = tipInfUt - tipTargUt

    tipVPerpDt = sp.vperp(tipDiffDt, tipTargDt)
    tipVPerpUt = sp.vperp(tipUiffUt, tipTargUt)
    self.VinfVtarg200sDiff_deg = (dpr
      * math.atan( max( [ sp.vnorm(tipVPerpDt) / sp.vnorm(tipTargDt)
                        , sp.vnorm(tipVPerpUt) / sp.vnorm(tipTargUt)
                        ] ) ) )

    ##################################################################
    ### Plump ellipse:
    ### * Use function FINDSATEPHUNCABC

    kxttrg = pbeStr.kinetx.Tca.Target

    ### - Square of Target uncertainties wrt Pluto barycenter

    abcTargStr = FINDSATEPHUNCABC(pbeStr.TargName
                                 ,sigmaRTN=pbeStr.satEphemUncert
                                 ,observer=pbeStr.ObsName
                                 ,tcaTarget=kxttrg
                                 ,etRtnOffset=self.ETminusTCA
                                 )
    abcTarg = abcTargStr.sigmaABC

    ### - adjustment to convert to uncertainty wrt Pluto

    ### - For multiple-body MU69 KEM cases, when TargID is 2486958nn,
    ###   use 2486958 as targIDclass to find correct clause here
    ### - Otherwise use targID as targIDclass

    targID = pbeStr.TargID

    if targID >= 248695800 and targID < 248695900: targIDclass = 2486958
    else                                         : targIDclass = targID

    if targIDclass == 2486958:
      abcPlu = sp.vpack(0,0,0)  ### MU69 single or multiple body:  no other required
      if iUTC == 0 and idEt == 0:
        pbeStr.pluEphemUncert = [0e0,0,0]
        pbeStr.pluEphemUncertSource = 'N/A'
        if targID != targIDclass and sp.vnorm(abcTarg) == 0e0:
          print('WARNING:  MU69 satellite barycenter-relative uncertainty is zero; consider using keyword argument satE')

    elif targIDclass == 999:
      abcPlu = sp.vpack(0,0,0)             ### Pluto:  no other required
      if iUTC == 0 and idEt == 0:
        pbeStr.pluEphemUncert = [0e0,0,0]
        pbeStr.pluEphemUncertSource = 'N/A'

    elif targIDclass == 901:
      ### Charon:  add 10% after RSS, 21% after squaring:
      ###   1.1^2 = 1.21 = 1^2 + sqrt(.21)^2
      abcPlu = abcTarg * numpy.sqrt(.21e0)
      if iUTC == 0 and idEt == 0:
        pbeStr.pluEphemUncert = abcTargStr.sigmaRTN / 10e0
        pbeStr.pluEphemUncertSource = 'sqrt(.21) of Target (Charon)'

    elif targIDclass > 901 and targIDclass < 999: ### Nix & Hydra:  add Pluto
      abcPlutoStr = FINDSATEPHUNCABC( '999'
                                    , observer=pbeStr.ObsName
                                    , sigmaRTN=pStE
                                    , tcaTarget=kxttrg
                                    , etRtnOffset=self.ETminusTCA
                                    )
      abcPlu = abcPlutoStr.sigmaABC
      if iUTC == 0 and idEt == 0:
        pbeStr.pluEphemUncert = abcPlutoStr.sigmaRTN
        pbeStr.pluEphemUncertSource = abcPlutoStr.sigmaRTNSource

    else:
      abcPlu = sp.vpack(0,0,0)  ### Single body?  No other required
      if iUTC == 0 and idEt == 0:
        pbeStr.pluEphemUncert = [0e0,0,0]
        pbeStr.pluEphemUncertSource = 'N/A'
        if targID != targIDclass and sp.vnorm(abcTarg) == 0e0:
          print('WARNING:  ???? satellite barycenter-relative uncertainty is zero; consider using keyword argument satE')


    ### RSS Plumped:  sqRt of Sum Squared; scaled by sigma multiple

    abcP = numpy.sqrt(abc0Sq + abcTarg*abcTarg + abcPlu*abcPlu
                     ) * pbeStr.sigmaMultiple

    self.plumpedABC = abcP

    self.scanVinfUptrack = PBEBUMP( -vInfU, -pTargU, abcP, pbeStr.TargRadius)
    self.scanVinfDowntrack = PBEBUMP( vInfU, -pTargU, abcP, pbeStr.TargRadius)
    self.scanVinfOvertrack = PBEBUMP( vInfPerpU, -pTargU, abcP, pbeStr.TargRadius)
    self.scanVinfBelowtrack = PBEBUMP( -vInfPerpU, -pTargU, abcP, pbeStr.TargRadius)

    self.scanVtargUptrack = PBEBUMP( vTargU, -pTargU, abcP, pbeStr.TargRadius)
    self.scanVtargDowntrack = PBEBUMP( -vTargU, -pTargU, abcP,pbeStr.TargRadius)
    self.scanVtargOvertrack = PBEBUMP( vTargPerpU, -pTargU, abcP, pbeStr.TargRadius)
    self.scanVtargBelowtrack = PBEBUMP( -vTargPerpU, -pTargU, abcP, pbeStr.TargRadius)

#######################################################################
### Test code

"""
!quiet=1

####################################
### Generate delta-ETs with exponential
###   steps between them as they move
###   away from the base UTC time
### - these delta-ETs will be the times
###     wrt the UTC base time (see below)
###     at which the calculations are
###     performed

dEts=(numpy.arange(3)-1)*spd
dEts=(numpy.arange(41)-20)*36e2
halfDEts = exp( dindgen(20) * alog(cspice_spd()) / 19e0)
dEts = [ -reverse(halfDEts), 0e0, halfDEts ]

pk='v_od059b.tm'       ### Pluto SPICE meta-kernel

mk= 'mu69altikore.tm'  ### MU69 SPICE meta kernel

####################################
### Call the PBECALCS.PRO routine
### - the base UTC time will default to the Pluto TCA
###   - this can be modified by specifying UTC
###     = e.g. , UTC='2008-07-14T10:20:30'
####################################

cRTN=[ 62e0,   22,   6]    ### Charon RTN ephemeris uncertainties wrt barycenter
nRTN=[109e0, 1629, 209]    ### Nix
hRTN=[ 64e0, 1512, 102]    ### Hydra

pRTN=cRTN / 10e0           ### Pluto <= 10% of Charon

mRTN=[ 64e0, 1512, 102]    ### Mu69

sMU69 = '2486958'

###        +-- INPUT Target body, STRING
###        |
###        |       +-- INPUT Kernels to load, STRING[]
###        |       |
###        |       |      +-- INPUT Delta ETs from base time, seconds, DOUBLE[]
###        |       |      |
###        |       |      |         +-- INPUT Sat ephem uncert, km, DOUBLE
###        |       |      |         |
###        |       |      |         |                        +-- OUTPUT
###        |       |      |         |                        |     Structure
###        |       |      |         |                        |     that contains
###        |       |      |         |                        |     calculated
###        |       |      |         |                        |     values
###        |       |      |         |                        |
###        |       |      |         |                        |          +-- Unload kernels
###        |       |      |         |                        |          |     on exit
###        |       |      |         |                        |          |
###        V       V      V         V                        V          V
pbecalcs,sMU69,ker=mk,dEts=dEts,satE=mRTN, debug=debug, pbeStr=pbmu69, /ukoe
pbecalcs,'999',ker=pk,dEts=dEts,satE=pRTN, debug=debug, pbeStr=pbe999, /ukoe
pbecalcs,'901',ker=pk,dEts=dEts,satE=cRTN, debug=debug, pbeStr=pbe901, /ukoe
pbecalcs,'902',ker=pk,dEts=dEts,satE=nRTN, debug=debug, pbeStr=pbe902, /ukoe
pbecalcs,'903',ker=pk,dEts=dEts,satE=hRTN, debug=debug, pbeStr=pbe903, /ukoe

####################################
### New case 2009-02-20:  Use Delivery control sigmas

pbecalcs,'903',kern=pk,dEts=dEts,satE=hRTN, debug=debug, pbeStr=pbeHDlv, /ukoe, /DelivSig

####################################
### Save all structures in IDL save file

save,pbmu69,pbe999,pbe901,pbe902,pbe903,pbeHDlv,fil='pbetest.idlsav'

####################################
### Convert all structures to TSV files
###   suitable for eXcel

### - put them in an array of structures

pbeStructs =
[ { name: 'mu69', pbePtr: ptr_new(pbmu69) }
, { name: 'charon', pbePtr: ptr_new(pbe901) }
, { name: 'nix', pbePtr: ptr_new(pbe902) }
, { name: 'hydra', pbePtr: ptr_new(pbe903) }
, { name: 'pluto', pbePtr: ptr_new(pbe999) }
, { name: 'hydra_deliverysigmas', pbePtr: ptr_new(pbeHDlv) }
]
for i=0L,n_elements(pbeStructs)-1L do begin
  ps=pbeStructs[i]
  pbe = *ps.pbePtr
  pbe2ascii,pbe,'pbetest'+ps.name+'.tsv'
  ###pbe2mdl,pbe,'pbetest'+ps.name+'.tsv'
  ptr_free, ps.pbePtr
endfor

ptr_free,ptr_new(pbeStructs,/no_copy),ptr_new(pbe,/no_copy)

####################################
### Again with zero s/c uncertainties for comparing Pluto and Charon

pbecalcs,'999',kern=pk, satE=pRTN, debug=debug, pbeStr=pbe999Zero, sigm=[1d-16,1d-16,1d-16], nSig=1, /ukoe
pbecalcs,'901',kern=pk, satE=cRTN, debug=debug, pbeStr=pbe901Zero, sigm=[1d-16,1d-16,1d-16], nSig=1, /ukoe
save,pbe999Zero,pbe901Zero,fil='pbetestzero.idlsav'

end

function pbecArg, arg, dflt, outVal=outVal, noErase=noErase, usedArg=usedArg
  ###xxx = keyword_set(noErase) ? ptr_new() : ptr_new(outVal,/no_copy)
  ###ptr_free, xxx
  if not keyword_set(noErase) then ptr_free, ptr_new(outVal,/no_copy)
  if n_elements(arg) gt 0L then begin
    outVal = arg
    usedArg = 1b
  endif else begin
    usedArg = 0b
    if n_elements(dflt) gt 0L then outVal = dflt
  endelse
  ### Impress default variable's type onto outVal
  vtyp=size(outVal,/typ)
  ptr_free, ptr_new(xxx,/no_copy)
  ntyp = size(xxx,/typ)              ### null type
  dtyp=size(dflt,/typ)
  if vtyp != dtyp and dtyp != ntyp then outVal = fix(outVal,type=dtyp)
  return, n_elements(outVal)
end
"""
