
data_directory="./Data"
remove_split_files=False	#removes raw surface velocity files spat out by Lare directly
overwrite_outputs=True	#if True will create combined time and velocity files from scratch, otherwise will append


import sys
sys.path[:0]=['/Change/This/Path']
from Permar_Functions_Python import *


merge_split_surface_velocity_coords_to_file(data_directory)
coordinate_timesteps,vx_timesteps,vy_timesteps,vz_timesteps=list_split_surface_velocity_files(data_directory)

if len(vx_timesteps)>0:
	x_grid,y_grid,processors,x_indices,y_indices=get_split_surface_velocity_coords(data_directory,vx_timesteps[0])
	merge_split_surface_velocity_to_file(data_directory,'x',processors,vx_timesteps[0],[len(y_grid),len(x_grid)],x_indices,y_indices,overwrite_outputs,ignore_ghost_cells=True)
for idx in range(1,len(vx_timesteps)):
	print(idx)
	x_grid,y_grid,processors,x_indices,y_indices=get_split_surface_velocity_coords(data_directory,vx_timesteps[idx])
	merge_split_surface_velocity_to_file(data_directory,'x',processors,vx_timesteps[idx],[len(y_grid),len(x_grid)],x_indices,y_indices,False,ignore_ghost_cells=True)
if len(vy_timesteps)>0:
	x_grid,y_grid,processors,x_indices,y_indices=get_split_surface_velocity_coords(data_directory,vy_timesteps[0])
	merge_split_surface_velocity_to_file(data_directory,'y',processors,vy_timesteps[0],[len(y_grid),len(x_grid)],x_indices,y_indices,overwrite_outputs,ignore_ghost_cells=True)
for idx in range(1,len(vy_timesteps)):
	x_grid,y_grid,processors,x_indices,y_indices=get_split_surface_velocity_coords(data_directory,vy_timesteps[idx])
	merge_split_surface_velocity_to_file(data_directory,'y',processors,vy_timesteps[idx],[len(y_grid),len(x_grid)],x_indices,y_indices,False,ignore_ghost_cells=True)
if len(vz_timesteps)>0:
	x_grid,y_grid,processors,x_indices,y_indices=get_split_surface_velocity_coords(data_directory,vz_timesteps[0])
	merge_split_surface_velocity_to_file(data_directory,'z',processors,vz_timesteps[0],[len(y_grid),len(x_grid)],x_indices,y_indices,overwrite_outputs,ignore_ghost_cells=True)
for idx in range(1,len(vz_timesteps)):
	x_grid,y_grid,processors,x_indices,y_indices=get_split_surface_velocity_coords(data_directory,vz_timesteps[idx])
	merge_split_surface_velocity_to_file(data_directory,'z',processors,vz_timesteps[idx],[len(y_grid),len(x_grid)],x_indices,y_indices,False,ignore_ghost_cells=True)

if remove_split_files:
	for f in os.listdir(data_directory):
		if os.path.isfile(os.path.join(data_directory, f)) and ('SurfV' in f or 'Surfx' in f or 'Surfy' in f) and f.count('_')==2:
			os.remove(os.path.join(data_directory, f))
		if os.path.isfile(os.path.join(data_directory, f)) and ('Surft' in f) and f.count('_')==1:
			os.remove(os.path.join(data_directory, f))


