from cc3d import CompuCellSetup
from Model2022Steppables import Model2022Steppable

CompuCellSetup.register_steppable(steppable=Model2022Steppable(frequency=1))
CompuCellSetup.run()
