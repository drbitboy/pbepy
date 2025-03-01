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
###   starting from uncertainties in orbit RTN frame (Radial,
###   Transverse, Normal to orbit plane)
###
### ARGUMENTS
###
###   Target       Satellite to use
###
### KEYWORDS
###
###   sigmaRTN     DOUBLE[3 or 1 or optional] Uncertainties in RTN frame
###                     - one value => along track only, [0,T,0]
###                     - defaults to values from RTNs dict
###   bFlag        Target type in local system; must be one of these:
###                - findsatephuncabc.PRIM - target is primary
###                - findsatephuncabc.CORR - correlated satellite
###                - findsatephuncabc.SATE - uncorrelated satellite
###                - None - defaults to value from RTNs dict
###   observer     STRING Observer to use, default to that for FINDTCA
###   tcaTarget    STRING Target to use as Base TCA to determine time of calc
###                     - defaults to Target argument
###                     - offset via etRTNoffset below
###   etRTNoffset  DOUBLE seconds of offset to add to Base TCA
###                     - positive => later
###                     - negative => earlier
###   Kernels      STRING[], kernels to load and unload
###   /debug       Turn on logging, turn off exception handler (CATCH)
########################################################################

########################################################################
### Dictionary of RTN uncertainties and Primary/Satellite flags
### - First three elements are default RTN values, km
### - Last element is body flag (PRIM, CORR, SATE)
########################################################################

### Body flags:  Primary; Correlated satellite; Uncorrelated satelite
PRIM,CORR,SATE,UNKN = 'PRIMARY CORRELATED UNCORRELATED UNKNOWN'.split()

### Pluto encounter
### - For Pluto, use one tenth of Charon values

RTNs = { 999: [  62e-1,   22e-1,   6e-1, PRIM]   ### Pluto
       , 901: [  62e00,   22e00,   6e00, CORR]   ### Charon
       , 902: [ 109e00, 1625e00, 209e00, SATE]   ### Nix
       , 903: [  64e00, 1512e00, 102e00, SATE]   ### Hydra

### MU69
### - For all cases, S/C uncertainty is relative to 2486958
### - For single body:
###   - 2486958 is SPICE ID for MU69
###   - 2486958 (MU69) positional uncertainty is included in S/C
###     uncertainty wrt itself, so it needs no additional
###     uncertainty
### - Multiple body cases have been dropped

         , 2486958: [   0e00,   0e00,    0e00, PRIM]  ### MU69/Arrokoth

#################################################################
### Lucy
### - Placeholder values ca. 2025-02-18:
###   - [10m, 100m, 20m] for primaries of binaries
###   - [1km, 10km, 2km] for secondaries of binaries
###   - [  0,   0,    0] for non-binaries
#################################################################
         , 920000617: [   1e-2,  10e00,    2e-2, PRIM]  ### Patroclus
         , 120000617: [   1e00,  10e00,    2e00, CORR]  ### Menoetius

         , 920003548: [   1e-2,  10e-2,    2e-2, PRIM]  ### Eurybates
         , 120003548: [   1e00,  10e00,    2e00, CORR]  ### Queta

         ,  20011351: [   0e00,   0e00,    0e00, PRIM]  ### Leucus

         , 920015094: [   1e-2,  10e-2,    2e-2, PRIM]  ### Polymele
         , 120015094: [   1e00,  10e00,    2e00, CORR]  ### Shaun

         ,  20021900: [   0e00,   0e00,    0e00, PRIM]  ### Orus

         ,  20052246: [   0e00,   0e00,    0e00, PRIM]  ### Donaldjohanson

         , 920152830: [   1e-2,  10e-2,    2e-2, PRIM]  ### Dinkinesh
         , 120152083: [   1e00,  10e00,    2e00, CORR]  ### Selam
         }

### - Multiple bodies' SPICE IDs are 2486958nn
### - We provide no nominal uncertainties for bodies 2486958nn here,
###   as they must be be provided via sigmaRTN argument
###   - sigmaRTN typically comes from sigEArg argument in PBECALCS
###   - If sigmaRTN is not provided, this function will return
###      a value fo False in the return object's .success attribute

class FINDSATEPHUNCABC:
  def __init__(self
              , target
              , bFlag=None
              , sigmaRTN=None
              , observer=None
              , tcaTarget=None
              , etRtnOffset=0e0
              , kernels=[]
              ):

    ### Target to determine TCA to use

    tcaTarget = str(target if tcaTarget is None else tcaTarget).strip()
    if tcaTarget is None: tcaTarget = str(target).strip()

    ftca_furnsh(kernels)      ### Load kernels

    ### Assume failure

    self.success = False
    self.status = 'FAILED'
    self.msg = ''
    self.extended_msg = ''
    self.bodyFlag = UNKN

    if not bFlag in (PRIM,CORR,SATE,None,):
      self.msg = 'Invalid body flag (bFlag keyword) value ({bFlag})'
      return

    ### Parse sigmaRTN argument i.e. uncertainties
    ### - three values => [sigmaR,sigmaT,sigamN]
    ### - one value => sigmaT => [0e0, sigmaT, 0e0]
    ### - else use target to lookup [sigmaR,sigmaT,sigmaN] from defaults
    ### - also store where the value came from

    if not (sigmaRTN is None):

      ### sigmaRTN argument is not the default of None:
      ### - Convert scalar (T) or sequence ([T] or [R,T,N]) to array

      sRTN = list(numpy.float64([sigmaRTN]).flatten())
      nRtn = len(sRTN)
      if nRtn == 3:
        source = 'User input three values:  [sigmaR, sigmaT, sigmaN]'
      elif nRtn == 1:
        sRTN = [0.,sRTN[0],0.]
        source = 'User input one value (transverse element):  [0e0, sigmaT, 0e0]'
      else:
        self.msg = 'User sigmaRTN input is not 1 or 3 numeric values'
        return

      ### Get target body flag one of (PRIM, CORR, SATE)

      bodyFlag = bFlag 
      if bFlag is None:
        ### Attempt to get default body flag from target name and RTNs
        try   : bodyFlag = RTNs[sp.bods2c(target)][-1]
        except: pass
      if not bodyFlag in (PRIM,CORR,SATE,):
        self.msg = 'Cannot determine body flag for [{target}]; bFlag keyword is [{bFlag}]'
        return

      sRTN.append(bodyFlag)

    else:

      ### sigmaRTN keyword argument is None (default)
      ### - Use target as lookup to retrieve default values
      ###   [sigmaR, sigmaT, sigmaN] from dict RTNs

      ### - convert target argument string to target ID integer for
      ###   lookup into RTNs

      try: sRTN = RTNs[sp.bods2c(target)]
      except:
        self.msg = f'No default RTN uncertainties available for target [{target}]]'
        self.extended_msg = traceback.format_exc()
        return

      source = 'Default RTN value(s) selected from target'

    ####################################################################
    ### Assign flag, sigmaRTN values and source attributes

    self.bodyFlag = sRTN[3] if bFlag is None else bFlag
    self.sigmaRTN = numpy.array(sRTN[:3])
    self.sigmaRTNSource = source

    ### Get xform from J2000 to ABC (uncertainty) frame at tcaTarget
    ### Get TCA of tcaTarget

    self.kinetx = k = kinetxcurrent(target=tcaTarget)
    self.TCA = tcaOut = FINDTCA(tcaTarget, observer)

    ### Get xform from RTN frame of target to J2000 frame
    ### Calculate xform from RTN frame to ABC frame

    self.Vinf = vinfOut = FINDVINF(tcaOut, target, etRtnOffset=etRtnOffset)
    self.mtx_RTN2Abc = sp.mxm(k.mJ2u,vinfOut.mtx_RTNTarg2j)

    ### Calculate Root-Sum-Squared (RSS) contribution from each RTN axis
    ### into ABC elements, put results into return structure
    ###
    ### .mtx_RTNTarg2Abc matrix xforms vectors from RTN frame
    ### to ABC frame.  Each row of that matrix comprises the RTN frame's
    ### XYZ vector components for one basis unit vector of the ABC
    ### frame; top row is X basis unit vector of ABC frame; middle row
    ### is Y; bottom row is Z.
    ###
    ### The products, v3x3, of each row's (each ABC basis vector's)
    ### elements of .mtx_RTNTarg2Abc with the RTN magnitudes yields the
    ### contributions of R, T, & N uncertainties along that basis vector

    ###     .mtx_RTN2Abc          vstack([.sigmaRTN]*3)
    ### [ XabcR  XabcT  XabcN ] . [ Rmag  Tmag  Nmag ]   => [XabcRMag  XabcTMag XabcNMag]
    ### [ YabcR  YabcT  YabcN ] . [ Rmag  Tmag  Nmag ]   => [YabcRMag  YabcTMag YabcNMag]
    ### [ ZabcR  ZabcT  ZabcN ] . [ Rmag  Tmag  Nmag ]   => [ZabcRMag  ZabcTMag ZabcNMag]

    v3x3 = self.mtx_RTN2Abc * numpy.vstack([self.sigmaRTN]*3)

    ### Squaring those v3x3 elements, then summing along the rows, and
    ### finally taking the square root of those sums, yields the RSS
    ### contributions of the RTN magnitudes along the ABC frame's basis
    ### vectors

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
