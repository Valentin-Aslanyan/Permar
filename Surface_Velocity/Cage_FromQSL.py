

data_dir="Data"
QSL_file="./QSL/0000/qslLare.qsl"


import sys
sys.path[:0]=['/Change/This/Path']
from Permar_Functions_Python import *


z_actual,y,x,Q=parse_QSL_Larebinfile(QSL_file)

x_out=np.array(x.flatten(),dtype=np.float64)
y_out=np.array(y.flatten(),dtype=np.float64)
write_surface_cage_file(x_out,y_out,os.path.join(data_dir,"SurfaceCage0.bin"))



