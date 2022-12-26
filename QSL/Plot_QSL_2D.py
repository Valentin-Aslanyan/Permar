
"""

"""

x_shift=850	#grid points
filename="./QSL/0000/qslLare_Pad.qsl"


import sys
sys.path[:0]=['/Change/This/Path']
from Permar_Functions_Python import *
import matplotlib.pyplot as plt
plt.rc('text', usetex=True)
plt.rc('font', family='serif')


z_actual,y_grid,x_grid,Q=parse_QSL_Larebinfile(filename)
max_Q=max(abs(Q).flatten())
Q_grid=np.copy(Q)
Q_grid[x_shift:,:]=Q[:np.shape(Q)[0]-x_shift,:]
Q_grid[:x_shift,:]=Q[np.shape(Q)[0]-x_shift:,:]
Q_grid=np.sign(Q_grid)*np.log(np.clip(abs(Q_grid),2.0,max_Q))/np.log(10.0)
#Q_grid=np.sign(Q)*np.log(abs(Q))/np.log(10.0)
Q_grid[np.isinf(Q_grid)]=np.nan


fig=plt.figure("Q 2D",figsize=(12,10))
ax=fig.gca()
plt.title(r"$z=0.0$",fontsize=20)
color_plot=plt.pcolormesh(x_grid,y_grid,Q_grid,cmap='RdBu_r',vmin=-5,vmax=5,rasterized=True)


plt.tick_params(axis='both', which='major',labelsize=19,direction='in',bottom=True, top=True, left=True, right=True)
plt.ylabel(r"$y$",fontsize=20)
plt.xlabel(r"$x$",fontsize=20)
cbar=fig.colorbar(color_plot)#,ticks=[-4,-3,-2,-1,0])
cbar.ax.tick_params(labelsize=19,direction='in', left=True, right=True)
cbar.set_label(label=r"$\mathrm{slog}(Q)$",fontsize=20)

plt.savefig("QSL.pdf", format="pdf", dpi=100,bbox_inches='tight',pad_inches=0.1)

plt.show()



