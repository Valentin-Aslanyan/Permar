
"""
Displays a connectivity map, see Connectivity_Map_Legend.py for colors
"""


connectivity_filename="./Data/Connectivity_0_1.bin"
title=r"$z=0\quad$ $t=\;1000$ s"


import sys
sys.path[:0]=['/Change/This/Path']
from Permar_Functions_Python import *
import matplotlib
import matplotlib.pyplot as plt
plt.rc('text', usetex=True)
plt.rc('font', family='serif')


connection_colors=[(171/255,210/255,229/255),(244/255,166/255,131/255),(102/255,102/255,102/255),(143/255,145/255,0/255),(1.0,128/255,0)]
connection_cmap=matplotlib.colors.ListedColormap(((171/255,210/255,229/255),(244/255,166/255,131/255),(102/255,102/255,102/255),(143/255,145/255,0/255),(1.0,128/255,0)))


y_grid,x_grid,connectivity_map=load_connectivity_map(connectivity_filename)
connectivity_map=np.floor(connectivity_map/2)

y_spacing=y_grid[0,1]-y_grid[0,0]
x_spacing=x_grid[1,0]-x_grid[0,0]

fig1=plt.figure("Connectivity map",figsize=(10,10))
ax1_1=fig1.gca()

color_plot1_1=plt.pcolormesh(x_grid-0.5*x_spacing,y_grid-0.5*y_spacing,connectivity_map+0.5,cmap=connection_cmap,vmin=0,vmax=5,rasterized=True)
plt.title(title,fontsize=46)
plt.tick_params(axis='both', which='major',labelsize=40,direction='in',bottom=True, top=True, left=True, right=True)
plt.ylabel(r"$y$ [a.u.]",fontsize=42)
plt.xlabel(r"$x$ [a.u.]",fontsize=42)
#plt.xlim(x_limits)
#plt.ylim(y_limits)

plt.savefig("Connectivity_Map.pdf", format="pdf", dpi=100,bbox_inches='tight',pad_inches=0.1)
plt.show()








