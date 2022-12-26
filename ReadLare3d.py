"""
Read in and plot basics from .sdf file output by Lare3d
partial_load=True will load in target_quantities and therefore conserve memory

"""


lare3d_filename="./Data/0000.sdf"
lare3d_datfilename="./Data/lare3d.dat"

partial_load=False
target_quantities=["grid","Bx","By","Bz","Rho","Temperature"]


import sys
sys.path[:0]=['/Change/This/Path']
from Permar_Functions_Python import *
import matplotlib.pyplot as plt
plt.rc('text', usetex=True)
plt.rc('font', family='serif')

if partial_load:
	lare3d_file_properties,lare3d_blocks=read_cfsasdf_file(lare3d_filename,target_quantities=target_quantities)
else:
	lare3d_file_properties,lare3d_blocks=read_cfsasdf_file(lare3d_filename)
lare3d_dat=read_lare3d_dat(lare3d_datfilename)

#Use regularize_Bfield(lare3d_blocks) for magnetic field and quantity_from_blocks(lare3d_blocks,...) for everything else
grid_x,grid_y,grid_z,B=regularize_Bfield(lare3d_blocks)
grid_x,grid_y,grid_z,Rho=quantity_from_blocks(lare3d_blocks,"Rho")
grid_x,grid_y,grid_z,temperature=quantity_from_blocks(lare3d_blocks,"Temperature")

print("Available data quantities    |    sizes")
for idx in range(len(lare3d_blocks)):
	print(lare3d_blocks[idx].block_id,len(lare3d_blocks[idx].data))
print()

#This will undo the normalization to give SI quantities e.g. [m], [T]
grid_x_SI=grid_x*lare3d_dat.L_norm
grid_y_SI=grid_y*lare3d_dat.L_norm
grid_z_SI=grid_z*lare3d_dat.L_norm
B_SI=B*lare3d_dat.B_norm
Rho_SI=Rho*lare3d_dat.Rho_norm
temperature_SI=temperature*lare3d_dat.temp_norm


#Plot and save figures
#Slight normalization to turn [m] -> [Mm] and [T] -> [G]
fig1=plt.figure("Density, Temperature lines",figsize=(12,10))
ax1=plt.gca()
plt.tick_params(axis='both', which='major',labelsize=20,direction='in',bottom=True, top=True, left=True, right=True)
plt.ylabel(r"$T$ [K]",fontsize=24,color="red")
plt.xlabel(r"$Z$ [Mm]",fontsize=24)
plt.plot(grid_z_SI*1E-6,temperature_SI[:,0,0],color="red",linewidth=3)
plt.yscale('log')
ax2=ax1.twinx()
plt.tick_params(axis='both', which='major',labelsize=20,direction='in',bottom=True, top=True, left=True, right=True)
plt.ylabel(r"$\rho$ [kg m$^{-3}$]",fontsize=24,color="blue")
plt.xlabel(r"$Z$ [Mm]",fontsize=24)
plt.plot(grid_z_SI*1E-6,Rho_SI[:,0,0],'--',color="blue",linewidth=3)
plt.yscale('log')
plt.savefig("Rho_T.pdf", format="pdf", dpi=100,bbox_inches='tight',pad_inches=0.1)

fig2=plt.figure("Magnetic field plane",figsize=(12,10))
ax2=plt.gca()
color_plot=plt.pcolormesh(grid_x_SI*1E-6,grid_y_SI*1E-6,np.sqrt(B_SI[0,:,:,0]**2+B_SI[0,:,:,1]**2+B_SI[0,:,:,2]**2)*1E4,rasterized=True,cmap="Oranges")
plt.tick_params(axis='both', which='major',labelsize=20,direction='in',bottom=True, top=True, left=True, right=True)
plt.xlabel(r"$X$ [Mm]",fontsize=24)
plt.ylabel(r"$Y$ [Mm]",fontsize=24)

cbar=fig2.colorbar(color_plot)
cbar.ax.tick_params(labelsize=20,direction='in', left=True, right=True)
cbar.set_label(label=r"$|B|$ [G]",fontsize=24)
plt.savefig("mod_B.pdf", format="pdf", dpi=100,bbox_inches='tight',pad_inches=0.1)

plt.show()





