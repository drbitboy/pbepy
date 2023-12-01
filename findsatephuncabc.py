import numpy
import traceback
import collections
import spiceypy as sp
from findtca import FINDTCA,ftca_furnsh
from findvinf import FINDVINF
from kinetxcurrent import kinetxcurrent
########################################################################
### findsatephuncabc.py
###
### Find satellite ephemeris uncertainties in ABC (ellipsoid) frame
###   starting with uncertainties in orbit RTN frame (Radial, Trajectory, 
###   orbit plane Normal)
###
### ARGUMENTS
###
###   Target       Satellite to use
###
### KEYWORDS
###
###   SigmaRTN     DOUBLE[3 or 1] Uncertainties in RTN frame
###                     - one value => along track only, [0,T,0]
###   Observer     STRING Observer to use, default to that for FINDTCA
###   TcaTarget    STRING Target to use as Base TCA to determine time of calc
###                     - defaults to Target argument
###                     - offset via etRTNoffset below
###   etRTNoffset  DOUBLE seconds of offset to add to Base TCA
###                     - positive => later
###                     - negative => earlier
###   Kernels      STRING[], kernels to load and unload
###   /debug       Turn on logging, turn off exception handler (CATCH)
########################################################################
class FINDSATEPHUNCABC:
  def __init__(self
              , target
              , sigmaRTN=None
              , observer=None
              , tcaTarget=None
              , etRtnOffset=0e0
              , kernels=[]
              ):

    ### Target to determine TCA to use

    tcaTarget = str(target if tcaTarget is None else tcaTarget).strip()

    ftca_furnsh(kernels)      ### Load kernels

    ### Exception handling

    self.success = False
    self.status = 'FAILED'
    self.msg = ''
    self.extended_msg = ''
    """
    er=0L
    if nodebug then catch,er
    if er ne 0L then begin
      catch,/cancel
      r.msg = !error_state.msg
      if not r.success then message,/cont, !error_state.msg
      message,/reset
      ftca_furnsh(kernels, unload=True)
      return, r
    endif
    """

    ### Parse sigmaRTN argument i.e. uncertainties
    ### - three values => [sigmaR,sigmaT,sigamN]
    ### - one value => sigmaT => [0e0, sigmaT, 0e0]
    ### - else use target to lookup [sigmaR,sigmaT,sigmaN] from default values 
    ### - also store where the value came from

    if not (sigmaRTN is None):
      sRTN = numpy.float64([sigmaRTN]).flatten()
      nRtn = len(sRTN)
      if nRtn == 3:
        source = 'User input three values:  [sigmaR, sigmaT, sigmaN]'
      elif nRtn == 1:
        sRTN = sp.vpack(0,sRTN[0],0)
        source = 'User input one value (along-track element):  [0e0, sigmaT, 0e0]'
      else:
        self.msg = 'User sigmaRTN input is not 1 or 3 numeric values'
        return

    else:

      ### Default values:  [sigmaR, sigmaT, sigmaN]

      rtnArr = { 901: [  62e0,   22e0,   6e0]
               , 902: [ 109e0, 1625e0, 209e0]
               , 903: [  64e0, 1512e0, 102e0]

      ### MU69
      ### - For all cases, S/C uncertainty is relative to 2486958
      ### - For single body:
      ###   - 2486958 is SPICE ID for MU69
      ###   - 2486958 (MU69) positional uncertainty is included in S/C
      ###     uncertainty wrt itself, so it needs no additional
      ###     uncertainty
      ### - For multiple bodies
      ###   - 2486958 is SPICE ID for MU69 system barycenter, 
      ###   - 2486958 (MU69 barycenter) positional uncertainty is included
      ###     in S/C uncertainty wrt itself, so it needs no additional
      ###     uncertainty

               , 2486958: [   0e0,   0e0,    0e0]

      ### - Multiple bodies'SPICE IDs are 2486958nn
      ### - We provide no nominal uncertainties for bodies 2486958nn here,
      ###   as they must be be provided via sigmaRTN argument
      ###   - sigmaRTN typically comes from sigEArg argument in PBECALCS
      ###   - If sigmaRTN is not provided, this function will throw an
      ###      exception via the MESSAGE call below:
      ###
      ###        message, 'No default RTN uncertainties available ...'
               }

      ### No default values for Pluto, use one tenth of Charon values

      rtnArr[999] = numpy.array(rtnArr[901]) /10e0   ### Pluto~Charon/10

      ### - convert target argument string to target ID

      try: sRTN = rtnArr[sp.bods2c(target)]
      except:
        self.msg = f'No default RTN uncertainties available for target [{target}]]'
        self.extended_msg = traceback.format_exc()
        return

      source = 'Default value selected from target'

    ### Add sigmaRTN values and source to return structure

    self.sigmaRTN = numpy.array(sRTN)
    self.sigmaRTNSource = source

    ### Call kinetxcurrent to get xform from J2000 to ABC (uncertainty) frame
    ### Call FINDTCA to get TCA of tcaTarget
    ### Call FINDVINF to get xform from RTN frame to J2000 frame

    k = kinetxcurrent(target=tcaTarget)
    tcaOut = FINDTCA(tcaTarget, observer)
    vinfOut = FINDVINF(tcaOut, target, etRtnOffset=etRtnOffset)

    self.kinetx = k
    self.TCA = tcaOut
    self.Vinf = vinfOut

    ### Matrix multiply to get xform from RTN to ABC frame

    self.mtx_RTN2Abc = sp.mxm(k.mJ2u,vinfOut.mtx_RTNTarg2j)

    ### RSS ABC elements from each RTN axis, put results into return structure

    v3x3 = self.mtx_RTN2Abc * numpy.vstack([self.sigmaRTN]*3)
    self.sigmaABC = numpy.sqrt( (v3x3 * v3x3).sum(axis=1) )

    ### Indicate success

    self.success = True
    self.status = self.msg = 'OK'

"""
end

########################################################################
### Test code

  !quiet=1b

  k='v_od059b.tm'   ### New Horizons meta-kernel

  ### Make array of RTN structures for Pluto system bodies

  tts=['999','901','902','903']
  
  dt = -120L ### * 0L

  for itt=0L,n_elements(tts)-1L do begin

  cdt = '+(' + strtrim(dt,2) + 's)'

  tt=tts[itt]
  rtnAbcArr = [ findsatephuncabc('901', tcaTarget=tt, k=k, etRtnOffset=dt) $
              , findsatephuncabc('902', tcaTarget=tt, k=k, etRtnOffset=dt) $
              , findsatephuncabc('903', tcaTarget=tt, k=k, etRtnOffset=dt) $
              , findsatephuncabc('999', tcaTarget=tt, k=k, etRtnOffset=dt) $
              ]

  if itt eq 0L then begin
    help,/st,rtnAbcArr
  endif
  print,f='(a)' $
       ,'','RTN Sigmas in ABC frame, Charon, Nix, Hydra, Pluto at TCA@'+tt+cdt
  print,rtnAbcArr.sigmaABC, f='(3f15.3)'

  ### RSS SQRT(KINETX eigenvalues), Body & Pluto ABC uncertainties

  RSSedABCus = sqrt( rtnAbcArr[0:2].sigmaABC^2 $
               + ( (rtnAbcARR[3].sigmaABC^2) xxx [1,1,1] ) $   ### Pluto
               + ( kget(rtnAbcArr[0].kinetx,/KnEig) xxx [1,1,1] ) $ ###Eig=>Squared
               )

  print,f='(a)', 'RSSed sigmas in ABC from Charon, Nix, Hydra at ' $
               , 'TCA@' + tt + cdt + '=' + rtnAbcArr[0].vInf.tdbRTNCaldate
  print,RSSedABCus, f='(3f15.3)'

  dt = dt * 2L

  endfor
end
"""
