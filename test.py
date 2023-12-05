import os
import pprint
import pbecalcs

os.system("kernels/meta_kernel.tm --get-kernels")

pc = pbecalcs.PBESTRUCT(kern=["kernels/meta_kernel.tm"])
pprint.pprint(dict(pbecalcs=vars(pc)
                  ,pbecalcs_kinetx=vars(pc.kinetx)
                  ,pbecalcs_pbeArr00=vars(pc.pbeArr[0,0])
                  ,pbecalcs_pbeArr00_scanVinfUptrack=vars(pc.pbeArr[0,0].scanVinfUptrack)
                  ,))
