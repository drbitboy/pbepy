import collections
import spiceypy as sp
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
, DelivSigmas=0  ### Misc keywords
):
    unloadKernelsOnExit = ukoe and True or False
    self.success = False
    self.status = 'PBEcalcs SKIPPED'
    self.msg = ''

    try:
      pbeStr.status = 'PBEcalcs FAILED to load kernels'
      for k in kern: sp.furnsh(k)
      pbeStr.status = 'PBEcalcs FAILED'

      self.TargName = str(targArg)
      self.ObsName = str(obsArg)

      self.TargID = sp.bods2c(targNameArg)
      self.ObsID = sp.bods2c(obsNameArg)

      try   : self.TargName = sp.bodc2n(self.targID,99)
      except: pass

      try   : self.ObsName = sp.bodc2n(self.obsID,99)
      except: pass

      if radi is None:
        try   : t = numpy.max(sp.bodvcd(self.targID,'RADII',4))
        except: t = 0e0
      else    : t = radi
      self.TargRadius = float(t)
 
      self.kinetx = kx = kinetxcurrent(target=targArg)

      ### KGET defaults to Knowledge Sigmas

      if sigm is None: self.baseSigmas = kx.kget(DeSig=DelivSigmas)
      else           : self.baseSigmas = list(map(float,sigm))

      if mJ2u is None: self.mtx_j2k2Uncert  = kx.mJ2u
      else           : self.mtx_j2k2Uncert = mJ2u

      self.sigmaMultiple = float(nSig)
      self.mdlFilename = mdlF
      self.dataFilename = datF

      try:
        L = len(satE)
        if L == 3: self.satEphemUncert = sp.vequ(*satE)
        else     : self.satEphemUncert = sp.vequ(00, satE[0], 0e0)
      except     : self.satEphemUncert = sp.vequ(00, satE, 0e0)

      self.pluEphemUncert = numpy.zeros(3,dtype=numpy.float64)
      self.pluEphemUncertSource = ''

      if isinstance(UTCsArg,str): UTCs = [UTCsArg]
      else:
        assert isinstance(UTCsArg,collections.Sequence)
        UTCs = UTCsArg
      nUTCs = len(UTCs)

      if isinstance(dEtsArg,str): UTCs = [UTCsArg]
      else:
        assert isinstance(UTCsArg,collections.Sequence)
        UTCs = UTCsArg
      nUTCs = len(UTCs)
      if nUTCs == 0:
        UTCs = [sp.et2utc(kinetx.tca.et,'ISOC',3)]
        nUTCs = 1

      try: dEts = [float(dEtsArg)]
      except: dEts = list(map(float,dEtsArg)))
      ndEts = len(dEts)
      if ndEts == 0:
        dEts = [0e0]
        ndEts = 1

    except:
      raise
"""

  pbeScanSample = pbebump(/sample)

  pbeBase =
  { UTC: ''
  , ETminusTCA: 0d0
  , nomVec_J2k: dblarr(3)
  , nomRRaDec_J2k: dblarr(3)
  , nomVinfPerp_J2k:  dblarr(3)
  , projVinfVtargRot_Deg: 0d0
  , plumpedABC: dblarr(3)
  , $  ### Scans up & down along Vinfinity...
    scanVinfUptrack: pbeScanSample
  , scanVinfDowntrack: pbeScanSample
  , $  ### Scans over & below across Vinfinity...
    scanVinfOvertrack: pbeScanSample
  , scanVinfBelowtrack: pbeScanSample
  , $  ### Scans up & down along Vtarget...
    scanVtargUptrack: pbeScanSample
  , scanVtargDowntrack: pbeScanSample
  , $  ### Scans over & below across Vtarget ...
    scanVtargOvertrack: pbeScanSample
  , scanVtargBelowtrack: pbeScanSample
  , VinfVtarg200sDiff_deg: 0d0 $      ### Offset @ 200s btw Vinf & Vtarg in FOV
  }

  pbeArr = replicate(pbeBase, nUTCs, ndEts)

  dpr = cspice_dpr()
  abc0Sq = pbeStr.baseSigmas*pbeStr.baseSigmas  ### Apply sigma multiplier later

  j2u = pbeStr.mtx_j2k2Uncert

  for iUTC=0L,nUTCs-1L do begin
    cspice_utc2et, UTCs[iUTC], et0
    for idEt=0L,ndEts-1L do begin
      pbe = pbeBase

      ### Time
      et = et0 + dEts[idEt]
      cspice_et2utc,et, 'ISOC',3L,utc
      pbe.UTC = utc
      pbe.ETminusTCA = et - kinetx.tca.et

      ### Ephemeris

      cspice_spkezr, targName, et, 'J2000', 'LT', obsName, stJ2k, ltim
      pbe.nomVec_J2k = stJ2k[0:2]
      cspice_recrad, stJ2k[0:2], r,ra,dec
      pbe.nomRRaDec_J2k = [ r, [ra,dec] * dpr]

      ### Perpendiculars to projection into FOV of Vinfinty & of Vtarget lines
      ### * stJ2k[0:2] is vector to nominal target position
      ### * .mtx_j2k2Uncert[*,0] is Vinfinity, pointing downtrack
      ### * stJ2k[3:5] is target velocity, Vtarg, pointing uptrack-ish
      ### * Orientation of FOV (in plane normal to nominal target position vec):
      ###   * Vinf or Vtarg is horizontal, Downtrack is Left
      ###   * VInfPerp and VTargPerp are Up

      vInf = pbeStr.mtx_j2k2Uncert[*,0]
      cspice_ucrss, stJ2k[0:2], vInf, vInfPerp
      pbe.nomVinfPerp_J2k = vInfPerp
      cspice_ucrss, stJ2k[3:5], stJ2k[0:2], vTargPerp
      pbe.projVinfVtargRot_Deg = cspice_vsep(vInfPerp,vTargPerp) * dpr

      ### Convert vectors to Uncertainty frame

      cspice_mxv, j2u, stJ2k[0:2], pTargU   ### Target position
      cspice_mxv, j2u, stJ2k[3:5], vTargU   ### Target velocity, uptrack
      cspice_mxv, j2u, vInf, vInfU          ### Vinfinity, downtrack
      cspice_mxv, j2u, vInfPerp, vInfPerpU
      cspice_mxv, j2u, vTargPerp, vTargPerpU

      ### Calculate offset btw Vinf & Vtarg at 200s in FOV in degrees

      tipTargDt = pTargU - (vTargU * 200d0)    ### 200s downtrack along vTarg
      tipTargUt = pTargU + (vTargU * 200d0)    ### 200s uptrack along vTarg
      l = cspice_vnorm(vTargU) * 200d0
      tipInfDt = pTargU + (vInfU * l)          ### 200s downtrack along vInfU
      tipInfUt = pTargU - (vInfU * l)          ### 200s uptrack along vInfU

      tipDiffDt = tipInfDt - tipTargDt
      tipUiffUt = tipInfUt - tipTargUt

      cspice_vperp, tipDiffDt, tipTargDt, tipVPerpDt
      cspice_vperp, tipUiffUt, tipTargUt, tipVPerpUt
      pbe.VinfVtarg200sDiff_deg = dpr
      * atan( max( [ cspice_vnorm(tipVPerpDt) / cspice_vnorm(tipTargDt)
                   , cspice_vnorm(tipVPerpUt) / cspice_vnorm(tipTargUt)
                   ] ) )

      ##################################################################
      ### OLD CODE:
      ##################################################################

      ##################################################################
      ### Plump ellipse:
      ### * Satellite ephem uncertainty times satellite velocity unit vector
      ### * RSS satellite ephemeris uncertainties
      ###     with ellipse uncertainties in Uncertainty frame

      ###cspice_vhat, vTargU, uvVTargU
      ###abcAdd = uvVTargU * pbeStr.satEphemUncert
      ##################################################################
      ### END OLD CODE
      ##################################################################

      ##################################################################
      ### Plump ellipse:
      ### * Use function findsatephuncabc

      seu = pbeStr.satEphemUncert
      ttrg = kinetx.tca.target


      ### - Square of Target uncertainties wrt Pluto barycenter

      abcTargStr = findsatephuncabc( targName
                   , sigmaRTN=pbeStr.satEphemUncert, observer=obsName
                   , tcaTarget=kinetx.tca.target, etRtnOffset=pbe.ETminusTCA
                   , debug=debug)
      abcTarg = abcTargStr.sigmaABC

      ### - adjustment to convert to uncertainty wrt Pluto

      ### - For multiple-body MU69 KEM cases, whenre targID is 2486958nn,
      ###   use 2486958 as targIDclass to find correct clause here
      ### - Otherwise use targID as targIDclass

      targIDclass = (targID ge 248695800 and targID lt 248695900) ? 2486958 : targID

      case targIDclass of
      2486958: begin
             abcPlu = [0d0,0,0]  ### MU69 single or multiple body:  no other required
             if iUtc eq 0 and idEt eq 0 then begin
               pbeStr.pluEphemUncert = [0d0,0,0]
               pbeStr.pluEphemUncertSource = 'N/A'
               if targID ne targIDclass and cspice_vnorm(abcTarg) eq 0d0 then begin
                 message,/continue,'WARNING:  MU69 satellite barycenter-relative uncertainty is zero; consider using keyword argument satEArg'
               endif
             endif
           end
      999: begin
             abcPlu = [0d0,0,0]                   ### Pluto:  no other required
             if iUtc eq 0 and idEt eq 0 then begin
               pbeStr.pluEphemUncert = [0d0,0,0]
               pbeStr.pluEphemUncertSource = 'N/A'
             endif
           end
      901: begin
             abcPlu = abcTarg / 10d0              ### Charon:  add 10%
             if iUtc eq 0 and idEt eq 0 then begin
               pbeStr.pluEphemUncert = abcTargStr.sigmaRTN / 10d0
               pbeStr.pluEphemUncertSource = '10% of Target (Charon)'
             endif
             ### Previous code is wrong and gets overridden unless
             ###   environment variable WRONGTENPERCENT is non-empty string
             ### Charon:  add 10% after RSS, 21% after squaring:
             ###   1.1^2 = 1.21 = 1^2 + sqrt(.21)^2
             if getenv("WRONGTENPERCENT") eq '' then begin
               abcPlu = abcTarg * sqrt(.21d0)
               if iUtc eq 0 and idEt eq 0 then begin
                 pbeStr.pluEphemUncert = abcTargStr.sigmaRTN / 10d0
                 pbeStr.pluEphemUncertSource = 'sqrt(.21) of Target (Charon)'
               endif
             endif
           end
      else: begin                               ### Nix & Hydra:  add Pluto
              abcPlutoStr = findsatephuncabc( '999'
                , observer=obsName
                , sigmaRTN=pStEArg
                , tcaTarget=kinetx.tca.target, etRtnOffset=pbe.ETminusTCA
                , debug=debug)
              abcPlu = abcPlutoStr.sigmaABC
              if iUtc eq 0 and idEt eq 0 then begin
                pbeStr.pluEphemUncert = abcPlutoStr.sigmaRTN
                pbeStr.pluEphemUncertSource = abcPlutoStr.sigmaRTNSource
              endif
            end
      endcase

      ### RSS Plumped:  sqRt of Sum Squared; scaled by sigma multiple

      abcP = sqrt(abc0Sq + abcTarg^2 + abcPlu^2) * sigmaMultiple


      pbe.plumpedABC = abcP

      pbe.scanVinfUptrack = pbebump( -vInfU, -pTargU, abcP, targRadius)
      pbe.scanVinfDowntrack = pbebump( vInfU, -pTargU, abcP, targRadius)
      pbe.scanVinfOvertrack = pbebump( vInfPerpU, -pTargU, abcP, targRadius)
      pbe.scanVinfBelowtrack = pbebump( -vInfPerpU, -pTargU, abcP, targRadius)

      pbe.scanVtargUptrack = pbebump( vTargU, -pTargU, abcP, targRadius)
      pbe.scanVtargDowntrack = pbebump( -vTargU, -pTargU, abcP,targRadius)
      pbe.scanVtargOvertrack = pbebump( vTargPerpU, -pTargU, abcP, targRadius)
      pbe.scanVtargBelowtrack = pbebump( -vTargPerpU, -pTargU, abcP, targRadius)

      ### Put results into array

      pbeArr[iUTC,idEt] = pbe
    endfor ### idEt
  endfor ### iUTC

  pbeStr = create_struct( pbeStr, 'pbeArr', pbeArr)

  pbeStr.success = 1b
  pbeStr.status = 'PBEcalcs OK'
  message,'OK'
end

#######################################################################
### Test code

!quiet=1

####################################
### Generate delta-ETs with exponential
###   steps between them as they move
###   away from the base UTC time
### - these delta-ETs will be the times
###     wrt the UTC base time (see below)
###     at which the calculations are
###     performed

dEts=(lindgen(3)-1L)*cspice_spd()
dEts=(lindgen(41)-20L)*36d2
halfDEts = exp( dindgen(20) * alog(cspice_spd()) / 19d0)
dEts = [ -reverse(halfDEts), 0d0, halfDEts ]

pk='v_od059b.tm'       ### Pluto SPICE meta-kernel

mk= 'mu69altikore.tm'  ### MU69 SPICE meta kernel

####################################
### Call the PBECALCS.PRO routine
### - the base UTC time will default to the Pluto TCA
###   - this can be modified by specifying UTC
###     = e.g. , UTC='2008-07-14T10:20:30'
####################################

cRTN=[ 62d0,   22,   6]    ### Charon RTN ephemeris uncertainties wrt barycenter
nRTN=[109d0, 1629, 209]    ### Nix
hRTN=[ 64d0, 1512, 102]    ### Hydra

pRTN=cRTN / 10d0           ### Pluto <= 10% of Charon

mRTN=[ 64d0, 1512, 102]    ### Mu69

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
  if vtyp ne dtyp and dtyp ne ntyp then outVal = fix(outVal,type=dtyp)
  return, n_elements(outVal)
end
"""
