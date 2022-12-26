
"""
Animate successive Q maps from QSL_squasher (not included in this repository)

Requires ffmpeg; written for Linux machines (possibly MacOS), needs altering for Windows
pad_start_frames and pad_end_frames will repeat first/last frames
"""


target_files=[
"./QSL/0000/qslLare.qsl",
"./QSL/0001/qslLare.qsl",
"./QSL/0002/qslLare.qsl",
"./QSL/0003/qslLare.qsl",
"./QSL/0004/qslLare.qsl",
"./QSL/0005/qslLare.qsl",
"./QSL/0006/qslLare.qsl",
"./QSL/0007/qslLare.qsl",
"./QSL/0008/qslLare.qsl",
"./QSL/0009/qslLare.qsl",
"./QSL/0010/qslLare.qsl",
"./QSL/0011/qslLare.qsl",
"./QSL/0012/qslLare.qsl",
"./QSL/0013/qslLare.qsl",
"./QSL/0014/qslLare.qsl",
"./QSL/0015/qslLare.qsl",
"./QSL/0016/qslLare.qsl",
"./QSL/0017/qslLare.qsl",
"./QSL/0018/qslLare.qsl",
"./QSL/0019/qslLare.qsl",
"./QSL/0020/qslLare.qsl"]

frames_per_step=3
frames_per_sec=10
pad_start_frames=3
pad_end_frames=3


import sys
sys.path[:0]=['/Change/This/Path']
from Permar_Functions_Python import *
import matplotlib.pyplot as plt
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
from subprocess import call



call_result=call(["mkdir","./anim_temp"])


plot_idx=0
for idx_f in range(len(target_files)):
	z_actual,y_grid,x_grid,Q=parse_QSL_Larebinfile(target_files[idx_f])
	max_Q=max(abs(Q).flatten())
	Q_grid=np.sign(Q)*np.log(np.clip(abs(Q),2.0,max_Q))/np.log(10.0)
	Q_grid[np.isinf(Q_grid)]=np.nan	

	plt.clf()
	fig=plt.figure("Output",figsize=(12,10))

	color_plot=plt.pcolormesh(x_grid,y_grid,Q_grid,cmap='RdBu_r',vmin=-5,vmax=5,rasterized=True)
	plt.tick_params(axis='both', which='major',labelsize=19,direction='in',bottom=True, top=True, left=True, right=True)
	plt.ylabel(r"$y$",fontsize=20)
	plt.xlabel(r"$x$",fontsize=20)
	cbar=fig.colorbar(color_plot)#,ticks=[-4,-3,-2,-1,0])
	cbar.ax.tick_params(labelsize=19,direction='in', left=True, right=True)
	cbar.set_label(label=r"$\mathrm{slog}(Q)$",fontsize=20)

	plt.savefig("./anim_temp/img{:03d}.png".format(plot_idx), format="png", dpi=100,bbox_inches='tight',pad_inches=0.1)
	plt.close(fig)
	plot_idx+=1
	if idx_f==0:
		for idx in range(pad_start_frames):
			call_result=call(["cp","./anim_temp/img{:03d}.png".format(plot_idx-1),"./anim_temp/img{:03d}.png".format(plot_idx)])
			plot_idx+=1
	for idx in range(1,frames_per_step):
		call_result=call(["cp","./anim_temp/img{:03d}.png".format(plot_idx-1),"./anim_temp/img{:03d}.png".format(plot_idx)])
		plot_idx+=1
	if idx_f==len(target_files)-1:
		for idx in range(pad_end_frames):
			call_result=call(["cp","./anim_temp/img{:03d}.png".format(plot_idx-1),"./anim_temp/img{:03d}.png".format(plot_idx)])
			plot_idx+=1


call_result=call(['ffmpeg -framerate '+str(frames_per_sec)+' -i ./anim_temp/img%03d.png -c:v libx264 -r 30 -pix_fmt yuv420p -vf "pad=ceil(iw/2)*2:ceil(ih/2)*2" ./anim_temp/anim.mp4'],shell=True)
call_result=call(["cp","./anim_temp/anim.mp4","./Animate_QSL2D.mp4"])
call_result=call(["rm","-r","./anim_temp/"])



