

data_dir="Data"
cage_x=[ 0.0, 0.0,0.0,0.0,0.0]
cage_y=[-0.7,-0.3,0.0,0.3,0.7]


import sys
sys.path[:0]=['/Change/This/Path']
from Permar_Functions_Python import *


x_out=np.array(cage_x,dtype=np.float64)
y_out=np.array(cage_y,dtype=np.float64)
write_surface_cage_file(x_out,y_out,os.path.join(data_dir,"SurfaceCage0.bin"))



