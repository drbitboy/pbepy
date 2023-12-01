import numpy
import spiceypy as sp
from findtca import FINDTCA
from findvinf import FINDVINF
from kinetxbase import KINETXBASE

class KINETX_MU69_20180626(KINETXBASE):
  """
Uncertainties class for MU69 flyby

 .provenance     ### Like it says
 .mJ2b           ### xform from J2000 to Bplane
 .mEig           ### xform from Bplane to Uncertainty frame
 .mJ2u           ### xform from J2000 to Uncertainty frame
 .sigmas         ### Knowledge:  Semi-major axes
 .eigVals        ### Knowledge:  Square of semi-major axes
 .delivSigmas    ### Delivery:  Semi-major axes
 .delivEigVals   ### Delivery:  Square of semi-major axes
 .Tca            ### TCA structure from FINDTCA.PRO
 .Vinf           ### Vinfinity structure from FINDVINF.PRO
 .KinetxCall     ### Name of this routine

"""
  def __init__(self):

    self.provenance = [ '20180626:  placeholder'
                      , ' - no matrix available yet'
                      , ' - Used Vinf to set A-hat, B-hat, C-hat'
                      , ' - Set delivery sigma to 1.5A, 2B, 4C'
                      ]

    ### Nominal 1-sigma values, NH S/C wrt MU69 or MU69_barycenter, km
    ### - These will typically be overidden by the [sigm=] keyword to
    ###   pbecalcs, so only the matrix will be used here
    ###                   TOF, s; Bnorm, km; Bmag, km
    self.sigmas = numpy.array([7215e0, 96e0,      105e0])

    ####################################################################
    ### Placeholders ca. 2017-11-08
    ### - We should not need to call FINDTCA/FINDVINF when we get explicit
    ###     values for the vectors for rawAHat, rawBHat and rawCHat
    ### - findving return struct vinfOut provides, in rows of mtx_j2b,
    ###   the following:
    ###
    ###                                                            |    kinetx*.pro
    ###                                          findvinf.pro      |    [ToF,Bnorm,Bmag]
    ###                                          B-plane (bpl)  <- | -> [A,B,C]
    ###       -----------------------------------------------------+------------------
    ###       mtx_j2b[*,2] = +Zbpl = [S/C vel wrt targ] = uvToF    |   uvToF   =  A-hat
    ###       mtx_j2b[*,1] = +Ybpl = [vTarget X uvToF]  = uvBdotT  |  -uvBnorm = -B-hat
    ###       mtx_j2b[*,0] = +Xbpl = [uvBdotR X uvToF]  = uvBdotR  |   uvBmag  =  C-hat

    ###   values B-plane A-hat, B-hat, C-hat

    TcaOut = FINDTCA( '2486958', '-98', utcEstArg='2019-01-01T12:00:00')
    vinfOut = FINDVINF( TcaOut, '2486958')

    ### ToF Axes of ellipse in J2000 frame
    rawCHat, rawBHat, rawAHat = vinfOut.mtx_j2b
    rawBHat -= 1.0
    ### +Zbpl = ToF(findvinf)   = ToF   = A-hat
    ### -Ybpl = BdotT(findvinf) = Bnorm = B-hat
    ### +Xbpl = BdotR(findving) = Bmag  = C-hat

    ### - Delivery 1-sigma values; [A, B, C]; units are s, km, km
    ###   - A, B and C delivery 1-sigma values scaled from knowledge
    ###     sigmas by 1.5, 2 and 4
    ###   - Division by TcaOut.speed converts A_km to A_s
    ### *** N.B. these are placholders

    delivTimeKmKm = self.sigmas * numpy.array([1.5/TcaOut.speed,2e0,4e0])

    ### End placeholder
    ######################################################################

    ######################################################################
    ### Future:  B-plane matrix and delivery denominator come from KinetX
    ###          document, typically an email, stored in array of strings
    ###          as the variable provenance, above
    ###
    ### It will likely look look similar to this in provenance:
    ###   ...
    ###   , '' $
    ###   , '  1-sigma requirement numbers:' $
    ###   , '' $
    ###   , '   |A| [km] = l*Vinf  1378	' $
    ###   , '   |B| [km] = n       66' $
    ###   , '   |C| [km] = m       33' $
    ###   , '   A-hat (x,y,z)   0.0113423      -0.9678945       -0.2511008' $
    ###   , '   B-hat (x,y,z)  -0.4532379       0.21886297      -0.8641033' $
    ###   , '   C-hat (x,y,z)   0.89131745      0.12360932      -0.436204' $
    ###   , '' $
    ###   ...
    ###
    ### rawZHat = [   0.0113423e0,   -0.9678945e0,   -0.2511008e0 ]  ### Axes of ellipse
    ### rawBHat = [  -0.4532379e0,    0.21886297e0,  -0.8641033e0 ]  ###   in J2000
    ### rawCHat = [   0.89131745e0,   0.12360932e0,  -0.436204e0  ]  ###     frame
    ###
    ### Delivery uncertainties are in Time (for A) and km (for B & C)
    ### - Will be sourced from provenance document
    ###
    ### delivTimeKmKm = [A_seconds, B_km, C_km]
    ######################################################################


    ######################################################################
    ######################################################################
    ### From here, all versions of kinetx_mu69_*.py will perform the
    ### same calculations

    ### Compile raw*Hat vectors into matrix

    self.mJ2u = sp.twovec(rawAHat,1,sp.ucrss(rawCHat, rawAHat),2)

    ### Knowledge sigmas' eigenvalues

    self.eigVals = self.sigmas*self.sigmas

    ### Alignment is along Vinfinity of spacecraft (-98) either wrt MU69
    ### (SPICE ID 2486958), or wrt MU69_barycenter (also SPICE ID 2486958),
    ### on 2019-01-01
    ### - N.B. this sets TcaOut.target to '2486958' either way, which is
    ###        used elsewhere, specifically in PBECALCS, which passes
    ###        kinetx.Tca.target as the value keyword argument tcaTarget
    ###        in calls to FINDSATEPHUNCABC

    ### May duplicate calls above
    self.Tca = FINDTCA( '2486958', '-98', utcEstArg='2019-01-01T12:00:00')
    self.Vinf = FINDVINF( self.Tca, '2486958')

    ### mEig => B-plane to Uncertainty:
    ###                                 T
    ### [B2u]= [J2u] [B2j] = [J2u] [J2b]
    self.mEig = sp.mxmt(self.mJ2u, self.Vinf.mtx_j2b)

    ### Calculate delivery (control) sigmas in (km, km, km) from
    ### delivery sigmas in (s, km, km)

    self.delivSigmas = delivTimeKmKm * numpy.array([TcaOut.speed,1e0,1e0])
    self.delivEigVals = self.delivSigmas * self.delivSigmas

    self.mJ2b = self.Vinf.mtx_j2b

    self.KinetxCall = __name__

"""
!quiet=1

### Test code
if n_elements(kinetx001kernels) eq 0L then begin
  kinetx001kernels = 'mu69altikore.tm'
  print,'Loading: ' + strjoin([kinetx001kernels],',')
  for i=0L,n_elements(kinetx001kernels)-1L do cspice_furnsh,kinetx001kernels[i]
endif

k=kinetx_mu69_20180626()
help,/st,k
print,'J2b:'
print,k.mJ2b
print,'Eigenvectors:'
print,k.mEig
print,'J2k to Uncertainty'
print,k.mJ2u
print,f='(a,3(/3f16.11,a))', 'J2u RA,DEC,R:' $
,cv_coord(from_r=k.mJ2u[*,0],/to_s,/deg),' = J2u[*,0]' $
,cv_coord(from_r=k.mJ2u[*,1],/to_s,/deg),' = J2u[*,1]' $
,cv_coord(from_r=k.mJ2u[*,2],/to_s,/deg),' = J2u[*,2]'
cspice_m2eul, k.mJ2u, 3L, 2L, 1L, az, ay, ax
print,f='(a,/3f16.11)', 'J2u Euler angles, X, Y, Z', [ax, ay, az] * cspice_dpr()
help,/st,k.tca

end
"""
