

show_nulls=True
lare3d_filenames=[
"./Data/0000.sdf",
"./Data/0010.sdf",
"./Data/0020.sdf",
"./Data/0030.sdf",
"./Data/0040.sdf",
"./Data/0050.sdf",
"./Data/0060.sdf",
"./Data/0070.sdf",
"./Data/0080.sdf"]
lare3d_datfilename="./Data/lare3d.dat"


fieldline_starts=[
[0.269,0.0,0.1],
[0.269,0.1,0.1],
[0.269,-0.1,0.1],
[0.369,0.0,0.1],
[0.169,0.0,0.1],
[0.274,0.05,0.1],
[0.264,0.05,0.1],
[0.274,-0.05,0.1],
[0.264,-0.05,0.1],
[0.269,0.0,0.2],
[0.269,0.1,0.2],
[0.269,-0.1,0.2],
[0.369,0.0,0.2],
[0.169,0.0,0.2],
[0.274,0.05,0.2],
[0.264,0.05,0.2],
[0.274,-0.05,0.2],
[0.264,-0.05,0.2],
[0.8,0.0,0.8]]

frames_per_step=10
frames_per_sec=10
pad_start_frames=3
pad_end_frames=3


import sys
sys.path[:0]=['/Change/This/Path']
from Permar_Functions_Python import *
if show_nulls:
	import HQVseg
import matplotlib.pyplot as plt
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
from mayavi import mlab
from subprocess import call


mlab.options.offscreen = True

times=np.zeros((len(lare3d_filenames)))
fieldlines=[]
if show_nulls:
	nulls=[]

lare3d_file_properties,lare3d_blocks=read_cfsasdf_file(lare3d_filenames[0],target_quantities=["grid","Bx","By","Bz"])
lare3d_dat=read_lare3d_dat(lare3d_datfilename)
times[0]=lare3d_file_properties.time*lare3d_dat.time_norm
grid_x,grid_y,grid_z,B=regularize_Bfield(lare3d_blocks)
grid_z=grid_z-grid_z[0]

grid_x_mg,grid_y_mg=np.meshgrid(grid_x,grid_y)
grid_z_mg=np.zeros(np.shape(grid_x_mg))
B_slices=np.zeros((len(lare3d_filenames),np.shape(B)[1],np.shape(B)[2]))
B_slices[0,:,:]=B[0,:,:,2]
fieldline_current=[]
for idx in range(len(fieldline_starts)):
	fieldline_current.append(shift_fieldline(trace_fieldline(fieldline_starts[idx],grid_x,grid_y,grid_z,B,max_steps=10000,dt=0.01,periodic_X=True,periodic_Y=True),grid_x,grid_y,grid_z))
fieldlines.append(fieldline_current)
if show_nulls:
	nulls_current=HQVseg.null_finder(len(grid_x),len(grid_y),len(grid_z),grid_x,grid_y,grid_z,np.transpose(B[:,:,:,0]),np.transpose(B[:,:,:,1]),np.transpose(B[:,:,:,2]))
	nulls_arr=np.zeros((len(nulls_current),3))
	for idx_n in range(len(nulls_current)):
		nulls_arr[idx_n,0]=nulls_current[idx_n]['x']
		nulls_arr[idx_n,1]=nulls_current[idx_n]['y']
		nulls_arr[idx_n,2]=nulls_current[idx_n]['z']
	nulls.append(nulls_arr)

for idx_f in range(1,len(lare3d_filenames)):
	lare3d_file_properties,lare3d_blocks=read_cfsasdf_file(lare3d_filenames[idx_f],target_quantities=["grid","Bx","By","Bz"])
	times[idx_f]=lare3d_file_properties.time*lare3d_dat.time_norm
	grid_x,grid_y,grid_z,B=regularize_Bfield(lare3d_blocks)
	B_slices[idx_f,:,:]=B[0,:,:,2]
	fieldline_current=[]
	for idx in range(len(fieldline_starts)):
		fieldline_current.append(shift_fieldline(trace_fieldline(fieldline_starts[idx],grid_x,grid_y,grid_z,B,max_steps=10000,dt=0.01,periodic_X=True,periodic_Y=True),grid_x,grid_y,grid_z))
	fieldlines.append(fieldline_current)
	if show_nulls:
		nulls_current=HQVseg.null_finder(len(grid_x),len(grid_y),len(grid_z),grid_x,grid_y,grid_z,np.transpose(B[:,:,:,0]),np.transpose(B[:,:,:,1]),np.transpose(B[:,:,:,2]))
		nulls_arr=np.zeros((len(nulls_current),3))
		for idx_n in range(len(nulls_current)):
			nulls_arr[idx_n,0]=nulls_current[idx_n]['x']
			nulls_arr[idx_n,1]=nulls_current[idx_n]['y']
			nulls_arr[idx_n,2]=nulls_current[idx_n]['z']
		nulls.append(nulls_arr)


R_max=1.2

call_result=call(["mkdir","./anim_temp"])


mlab.figure(bgcolor=(1.0,1.0,1.0),size=(1200,1000))

mlab.clf()

for idx in range(len(fieldlines[0])):
	mlab.plot3d(fieldlines[0][idx][:,0],fieldlines[0][idx][:,1],fieldlines[0][idx][:,2],line_width=0.01,color=(0/255,0/255,0/255),tube_radius=0.002)

plane_mesh=mlab.mesh(grid_x_mg,grid_y_mg,grid_z_mg,scalars=B_slices[idx_f,:,:],colormap='PiYG',vmin=-max(abs(B[0,:,:,2].flatten())),vmax=max(abs(B[0,:,:,2].flatten())))

mlab.points3d(1.1*np.array([R_max,R_max,R_max,R_max,-R_max,-R_max,-R_max,-R_max]), 1.1*np.array([R_max,R_max,-R_max,-R_max,R_max,R_max,-R_max,-R_max]), 1.1*np.array([R_max,-R_max,R_max,-R_max,R_max,-R_max,R_max,-R_max]),opacity=0.0)

mlab.view(azimuth=90.0, elevation=70.0, distance=4.0, focalpoint=(0.0,0.0,0.4))#, roll=None, reset_roll=True, figure=None)

mlab.savefig("./anim_temp/img{:03d}.png".format(0))

fig=mlab.gcf()
mlab.savefig("./anim_temp/img{:03d}.png".format(0))

plot_idx=0
for idx_f in range(len(lare3d_filenames)):
	mlab.clf()
	fig=mlab.gcf()

	if show_nulls:
		null_points=mlab.points3d(nulls[idx_f][:,0],nulls[idx_f][:,1],nulls[idx_f][:,2],scale_factor=0.025,color=(1.0, 0.0, 0.0))

	for idx in range(len(fieldlines[idx_f])):
		mlab.plot3d(fieldlines[idx_f][idx][:,0],fieldlines[idx_f][idx][:,1],fieldlines[idx_f][idx][:,2],line_width=0.01,color=(0/255,0/255,0/255),tube_radius=0.002)


	plane_mesh=mlab.mesh(grid_x_mg,grid_y_mg,grid_z_mg,scalars=B_slices[idx_f,:,:],colormap='PiYG',vmin=-max(abs(B[0,:,:,2].flatten())),vmax=max(abs(B[0,:,:,2].flatten())))

	mlab.points3d(1.1*np.array([R_max,R_max,R_max,R_max,-R_max,-R_max,-R_max,-R_max]), 1.1*np.array([R_max,R_max,-R_max,-R_max,R_max,R_max,-R_max,-R_max]), 1.1*np.array([R_max,-R_max,R_max,-R_max,R_max,-R_max,R_max,-R_max]),opacity=0.0)

	mlab.view(azimuth=90.0, elevation=70.0, distance=4.0, focalpoint=(0.0,0.0,0.4))#, roll=None, reset_roll=True, figure=None)

	mlab.savefig("./anim_temp/img{:03d}.png".format(plot_idx))
	plot_idx+=1
	if idx_f==0:
		for idx in range(pad_start_frames):
			call_result=call(["cp","./anim_temp/img{:03d}.png".format(plot_idx-1),"./anim_temp/img{:03d}.png".format(plot_idx)])
			plot_idx+=1
	for idx in range(1,frames_per_step):
		call_result=call(["cp","./anim_temp/img{:03d}.png".format(plot_idx-1),"./anim_temp/img{:03d}.png".format(plot_idx)])
		plot_idx+=1
	if idx_f==len(lare3d_filenames)-1:
		for idx in range(pad_end_frames):
			call_result=call(["cp","./anim_temp/img{:03d}.png".format(plot_idx-1),"./anim_temp/img{:03d}.png".format(plot_idx)])
			plot_idx+=1

call_result=call(['ffmpeg -framerate '+str(frames_per_sec)+' -i ./anim_temp/img%03d.png -c:v libx264 -r 30 -pix_fmt yuv420p -vf "pad=ceil(iw/2)*2:ceil(ih/2)*2" ./anim_temp/anim.mp4'],shell=True)
call_result=call(["cp","./anim_temp/anim.mp4","./Animate_Lare3d_Fieldlines3D.mp4"])
call_result=call(["rm","-r","./anim_temp/"])





