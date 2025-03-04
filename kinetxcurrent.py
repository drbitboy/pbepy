import spiceypy as sp
from kinetx_mu69_20250210 import KINETX_MU69
from kinetx_orus_20250210 import KINETX_ORUS
from kinetx_leucus_20250210 import KINETX_LEUCUS
from kinetx_polymele_20250210 import KINETX_POLYMELE
from kinetx_patroclus_20250210 import KINETX_PATROCLUS
from kinetx_dinkinesh_20250210 import KINETX_DINKINESH
from kinetx_eurybates_20250210 import KINETX_EURYBATES
from kinetx_donaldjohanson_20250210 import KINETX_DONALDJOHANSON
########################################################################
### Wrapper for current info from Kinetx
########################################################################
### 20180626:  Removed Pluto
def kinetxcurrent(target=None,initialize=True):

  sTarget = str(2486958 if target is None else target).strip()

  try: iTarget = sp.bods2c(sTarget)
  except:
    import traceback
    traceback.print_exc()
    assert False,f'Could not understand target [{sTarget}]'

  ######################################################################
  ### Lucy targets
  ######################################################################

  KMODULE = False

  ### ORUS                => SPICE ID 20021900
  if (iTarget % 100000000) == 20021900:
    KMODULE = KINETX_ORUS

  ### Leucus              => SPICE ID 20011351
  if (iTarget % 100000000) == 20011351:
    KMODULE = KINETX_LEUCUS

  ### Polymele            => SPICE ID 920015094
  ### Polymele_Barycenter => SPICE ID  20015094
  ### Shaun(?)            => SPICE ID 120015094
  if (iTarget % 100000000) == 20015094:
    KMODULE = KINETX_POLYMELE

  ### Patroclus            => SPICE ID 920000617
  ### Patroclus_Barycenter => SPICE ID  20000617
  ### Menoetius            => SPICE ID 120000617
  if (iTarget % 100000000) == 20000617:
    KMODULE = KINETX_PATROCLUS

  ### Dinkinesh            => SPICE ID 920152830
  ### Dinkinesh_Barycenter => SPICE ID  20152830
  ### Selam                => SPICE ID 120152830
  if (iTarget % 100000000) == 20152830:
    KMODULE = KINETX_DINKINESH

  ### DonaldJohanson => SPICE ID 20052246
  if (iTarget % 100000000) == 20052246:
    KMODULE = KINETX_DONALDJOHANSON

  ### Eurybates            => SPICE ID 920003458
  ### Eurybates_Barycenter => SPICE ID  20003458
  ### Queta                => SPICE ID 120003458
  if (iTarget % 100000000) == 20003548:
    KMODULE = KINETX_EURYBATES

  ######################################################################
  ### New Horizon targets
  ######################################################################

  ### MU69            => SPICE ID 2486958
  ### MU69_BARYCENTER => SPICE ID 2486959
  ### - Multiple bodies SPICE IDs are 218695801, 148695802, etc.
  if iTarget == 2486958 or (iTarget//100) == 2486958:
    KMODULE = KINETX_MU69

  assert KMODULE, f'Target {sTarget} is unrelated to a known target'
  if initialize: return KMODULE()
  return KMODULE
