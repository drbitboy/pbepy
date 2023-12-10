import os
import pprint
import pbecalcs

os.system("kernels/meta_kernel.tm --get-kernels")

pc = pbecalcs.PBESTRUCT(kern=["kernels/meta_kernel.tm"])
pprint.pprint(dict(pbecalcs=vars(pc)
                  ,pbecalcs_kinetx=vars(pc.kinetx)
                  ,pbecalcs_pbeArr00=vars(pc.pbeArr[0,0])
                  ,yw_scanVinfUptrack=vars(pc.pbeArr[0,0].scanVinfUptrack)
                  ,yx_scanVinfDowntrack=vars(pc.pbeArr[0,0].scanVinfDowntrack)
                  ,yy_scanVinfOvertrack=vars(pc.pbeArr[0,0].scanVinfOvertrack)
                  ,yz_scanVinfBelowtrack=vars(pc.pbeArr[0,0].scanVinfBelowtrack)
                  ,zw_scanVtargUptrack=vars(pc.pbeArr[0,0].scanVtargUptrack)
                  ,zx_scanVtargDowntrack=vars(pc.pbeArr[0,0].scanVtargDowntrack)
                  ,zy_scanVtargOvertrack=vars(pc.pbeArr[0,0].scanVtargOvertrack)
                  ,zz_scanVtargBelowtrack=vars(pc.pbeArr[0,0].scanVtargBelowtrack)
                  ,))
