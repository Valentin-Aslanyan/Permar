

input_directory='./Data'
output_directory='./QSL'

zero_z=True
write_axis_txt=True
write_axis_bin=True
write_bfld_txt=False
write_bfld_bin=True


import sys
sys.path[:0]=['/Change/This/Path']
from Permar_Functions_Python import *


infiles = [f for f in os.listdir(input_directory) if os.path.isfile(os.path.join(input_directory, f)) and '.sdf' in f]

for fname in infiles:
	lare3d_file_properties,lare3d_blocks=read_cfsasdf_file(os.path.join(input_directory,fname),target_quantities=["grid","Bx","By","Bz"])
	grid_x,grid_y,grid_z,B=regularize_Bfield(lare3d_blocks)
	if zero_z:
		grid_z=grid_z-grid_z[0]
	print(np.shape(B))

	Bx = B[:,:,:,0].transpose(2,1,0).flatten(order='F').astype('float64')
	By = B[:,:,:,1].transpose(2,1,0).flatten(order='F').astype('float64')
	Bz = B[:,:,:,2].transpose(2,1,0).flatten(order='F').astype('float64')

	pathnm = os.path.join(output_directory,fname.replace('.sdf',''))
	if os.path.exists(pathnm):
		print('path exists for ',pathnm)
	else:
		print('generating path ',pathnm)
		os.mkdir(pathnm)
	if write_axis_bin:
		grid_x.tofile(os.path.join(pathnm,'xs0.dat'))
		grid_y.tofile(os.path.join(pathnm,'ys0.dat'))
		grid_z.tofile(os.path.join(pathnm,'zs0.dat'))

	if write_axis_txt:
		np.savetxt(os.path.join(pathnm,'xs0_txt.dat'), grid_x)
		np.savetxt(os.path.join(pathnm,'ys0_txt.dat'), grid_y)
		np.savetxt(os.path.join(pathnm,'zs0_txt.dat'), grid_z)

	if write_bfld_bin:
		Bx.tofile(os.path.join(pathnm,'bx0.dat'))
		By.tofile(os.path.join(pathnm,'by0.dat'))
		Bz.tofile(os.path.join(pathnm,'bz0.dat'))

	if write_bfld_txt:
		np.savetxt(os.path.join(pathnm,'bx0_txt.dat'), Bx)
		np.savetxt(os.path.join(pathnm,'by0_txt.dat'), By)
		np.savetxt(os.path.join(pathnm,'bz0_txt.dat'), Bz)

	np.savetxt(os.path.join(pathnm,'dim.txt'), np.array((len(grid_z), len(grid_y), len(grid_x))).astype('int64'), fmt='%i')
	



