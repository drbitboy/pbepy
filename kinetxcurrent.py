import spiceypy as sp
from kinetx_mu69_20180626 import KINETX_MU69_20180626
from kinetx_eurybates_20231211 import KINETX_EURYBATES_20231211
from kinetx_donaldjohanson_20250210 import KINETX_DONALDJOHANSON_20250210
########################################################################
### Wrapper for current info from Kinetx
########################################################################
### 20180626:  Removed Pluto
def kinetxcurrent(target=None):

  sTarget = str(2486958 if target is None else target).strip()

  try: iTarget = sp.bods2c(sTarget)
  except:
    import traceback
    traceback.print_exc()
    assert False,f'Could not understand target [{sTarget}]'

  ### DonaldJohanson => SPICE ID 20052246
  if iTarget == 20052246:
    return KINETX_DONALDJOHANSON_20250210()

  ### Eurybates             => SPICE ID  20003458
  ### Eurybaters_Barycenter => SPICE ID 920003458
  ### Queta                 => SPICE ID 120003458
  if (iTarget % 100000000) == 20003548:
    return KINETX_EURYBATES_20250210()

  ### Polymele            => SPICE ID  20015094
  ### Polymele_Barycenter => SPICE ID 920015094
  ### Shaun(?)            => SPICE ID 120015094
  if (iTarget % 100000000) == 20015094:
    return KINETX_POLYMELE_20250210()

  ### MU69            => SPICE ID 2486958
  ### MU69_BARYCENTER => SPICE ID 2486959
  ### - Multiple bodies SPICE IDs are 218695801, 148695802, etc.
  if iTarget == 2486958 or (iTarget//100) == 2486958:
    return KINETX_MU69_20250210()    ### 2025-02-10

  assert False, f'Target {sTarget} is unrelated to MU69'
