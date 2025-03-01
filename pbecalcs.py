import os
import math
import numpy
import collections
import spiceypy as sp
from pbebump import PBEBUMP
from kinetxcurrent import kinetxcurrent
from findsatephuncabc import FINDSATEPHUNCABC,PRIM,CORR,SATE,UNKN

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
, targArg='POLYMELE' ### STRING, target, default = 'POLYMELE'
, obsArg='LUCY'  ### STRING, default = '-98' (New Horizons)
, kern=[]        ### STRING[], SPICE kernels to load and unload
, ukoe=False     ### /ukoe = Unload SPICE Kernels On Exit
, radi=None      ### DOUBLE, Target radius for bumping, km
, UTCsArg=[]     ### STRING[], UTCs for which to calculate info
, dEtsArg=[]     ### DOUBLE, delta ET to add to UTCs for times at which to calc
, sigm=None      ### DOUBLE[3], S/C ephemeris error wrt primary, km; dflt=>use kinetxcurrent(target=targArg)
, mJ2u=None      ### DOUBLE[3,3], J2k to T],[LTOF]], dflt=> findvinf()
, nSig=2e0       ### DOUBLE, Sigma multiple to use, dflt=2
, satE=None      ### DOUBLE, target ephemeris error wrt bary, target as non-primary only, km, dflt=[0,0,0]
, pStE=None      ### DOUBLE, primary ephemeris error wrt bary, target as non-primary only, km, dflt=findsatephuncabc
, mdlF=None      ### STRING, => MDL file name to which to write PBE for STK
, datF=None      ### STRING, => Flat ASCII filename to which to write info
, DelivSigmas=False ### Misc keywords
):
    if 'DEBUG' in os.environ:
      import pprint
      pprint.pprint(dict(PBESTRUCT_locals=locals()))
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
      else           : self.baseSigmas = numpy.array(list(map(float,sigm)),dtype=numpy.float64)

      if mJ2u is None: self.mtx_j2k2Uncert  = kinetx.mJ2u
      else           : self.mtx_j2k2Uncert = mJ2u

      self.sigmaMultiple = float(nSig)
      self.mdlFilename = mdlF
      self.dataFilename = datF

      self.satEphemUncertArg = self.satEphemUncert = satE

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

    if unloadKernelsOnExit: [sp.unload(k) for k in kern]


########################################################################
class PBE:
  """
Combine one, or three, flyby uncertainties:
  abc0      - Spacecraft-Primary*; from pbeStr.baseSigmas
  abcPri    - Primary-Barycenter**; RTN RSSed into ABC frame
  abcTarg   - Barycenter-Target**; RTN RSSed into ABC frame

  * always used
  ** only used if Target is not primary

"""
  def __init__(self,iUTC,idEt,pbeStr,pStE=None):

    ####################################################################
    ### Spacecraft-primary uncertainties
    ####################################################################
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
    ### * Use class FINDSATEPHUNCABC

    kxttrg = pbeStr.kinetx.Tca.Target

    ### - Square of Target uncertainties wrt system barycenter

    ####################################################################
    ### Target-Barycenter uncertainties
    ####################################################################
    abcTargStr = FINDSATEPHUNCABC(pbeStr.TargName
                                 ,sigmaRTN=pbeStr.satEphemUncertArg
                                 ,observer=pbeStr.ObsName
                                 ,tcaTarget=kxttrg
                                 ,etRtnOffset=self.ETminusTCA
                                 )
    assert abcTargStr.success
    pbeStr.satEphemUncert = abcTargStr.sigmaRTN
    abcTarg = abcTargStr.sigmaABC

    ### - adjustment to convert to uncertainty wrt Pluto/Primary

    if abcTargStr.bodyFlag == PRIM:
      ### Target is Pluto/Primary; no body uncertainties required
      abcTarg = sp.vpack(0,0,0)
      abcPri = sp.vpack(0,0,0)
      if iUTC == 0 and idEt == 0:
        pbeStr.pluEphemUncert = [0e0,0,0]
        pbeStr.pluEphemUncertSource = 'N/A'

    elif abcTargStr.bodyFlag == SATE:

      ### Uncorrelated target is not Primary:  RSSum uncertainties
      ### Convert Primary RTN uncertainties to B-Plane (ABC) frame

      abcPriStr = FINDSATEPHUNCABC(pbeStr.kinetx.PRIMARY
                                  ,sigmaRTN=pStE
                                  ,observer=pbeStr.ObsName
                                  ,tcaTarget=kxttrg
                                  ,etRtnOffset=self.ETminusTCA
                                  )
      assert abcPriStr.success

      ### Uncorrelated Primary RTN uncertainties in B-Plane frame
      abcPri = abcPriStr.sigmaABC

      if iUTC == 0 and idEt == 0:
        pbeStr.pluEphemUncert = list(abcPriStr.sigmaRTN)
        pbeStr.pluEphemUncertSource = 'Primary uncertainties Root Sum-of-Square with uncorrelated target'

    elif abcTargStr.bodyFlag == CORR:

      ### Correlated Primary RTN uncertainties in B-Plane frame,
      ### directly summed with satellite RTN uncertainties

      abcPriStr = FINDSATEPHUNCABC(pbeStr.kinetx.PRIMARY
                                  ,sigmaRTN=pStE
                                  ,observer=pbeStr.ObsName
                                  ,tcaTarget=kxttrg
                                  ,etRtnOffset=self.ETminusTCA
                                  )
      assert abcPriStr.success

      ### Calculate magnitudes of primary and target uncertanties
      pMag = sp.vnorm(abcPriStr.sigmaRTN)
      tMag = sp.vnorm(abcTargStr.sigmaRTN)

      ### 1) Target uncertainties are zero:  use primary uncertainties
      ### 2) Non-0 correlated target & primary:  sum uncerts via scaling
      ### 3) Else primary uncerts are 0:  use target uncerts as-is
      if tMag == 0.0 : abcTarg = abcPriStr.sigmaABC
      elif pMag > 0.0: abcTarg = sp.vscl(1.0+(pMag/tMag), abcTarg)

      ### Since primary uncertainties have either been folded into the
      ### target uncertainties, or are zero, use zero for the local
      ### primary uncertainties below
      abcPri = sp.vpack(0,0,0)

      if iUTC == 0 and idEt == 0:
        pbeStr.pluEphemUncert = list(abcPriStr.sigmaRTN)
        pbeStr.pluEphemUncertSource = 'Primary uncertainties direct-summed with correlated target, not Root-Sum-of-Squares'

    else: ### abcTargStr.bodyFlag == UNKN:
      ### Unknown primary/satellite/barycenter configuration
      assert abcTargStr.bodyFlag == UNKN
      abcTarg = sp.vpack(0,0,0)
      abcPri = sp.vpack(0,0,0)
      if iUTC == 0 and idEt == 0:
        pbeStr.pluEphemUncert = [0e0,0,0]
        pbeStr.pluEphemUncertSource = 'Body uncertainties not used:  unknown primary/satellite/barycenter configuration'

    ### RSS Plumped:  sqRt of Sum Squared; scaled by sigma multiple

    abcP = numpy.sqrt(abc0Sq + abcTarg*abcTarg + abcPri*abcPri
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
