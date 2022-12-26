
"""
Animate successive connectivity maps, see Connectivity_Map_Legend.py for colors

Requires ffmpeg; written for Linux machines (possibly MacOS), needs altering for Windows
pad_start_frames and pad_end_frames will repeat first/last frames
"""


connectivity_filenames=[
"./Data/Connectivity_0_0.bin",
"./Data/Connectivity_0_1.bin",
"./Data/Connectivity_0_2.bin",
"./Data/Connectivity_0_3.bin",
"./Data/Connectivity_0_4.bin",
"./Data/Connectivity_0_5.bin",
"./Data/Connectivity_0_6.bin",
"./Data/Connectivity_0_7.bin",
"./Data/Connectivity_0_8.bin",
"./Data/Connectivity_0_9.bin",
"./Data/Connectivity_0_10.bin",
"./Data/Connectivity_0_11.bin",
"./Data/Connectivity_0_12.bin",
"./Data/Connectivity_0_13.bin",
"./Data/Connectivity_0_14.bin",
"./Data/Connectivity_0_15.bin",
"./Data/Connectivity_0_16.bin",
"./Data/Connectivity_0_17.bin"]

titles=[
r"$z=0\quad$ $t=\;\;\;\;0$ s",
r"$z=0\quad$ $t=\;1000$ s",
r"$z=0\quad$ $t=\;2000$ s",
r"$z=0\quad$ $t=\;3000$ s",
r"$z=0\quad$ $t=\;4000$ s",
r"$z=0\quad$ $t=\;5000$ s",
r"$z=0\quad$ $t=\;6000$ s",
r"$z=0\quad$ $t=\;7000$ s",
r"$z=0\quad$ $t=\;8000$ s",
r"$z=0\quad$ $t=\;9000$ s",
r"$z=0\quad$ $t=10000$ s",
r"$z=0\quad$ $t=11000$ s",
r"$z=0\quad$ $t=12000$ s",
r"$z=0\quad$ $t=13000$ s",
r"$z=0\quad$ $t=14000$ s",
r"$z=0\quad$ $t=15000$ s",
r"$z=0\quad$ $t=16000$ s",
r"$z=0\quad$ $t=17000$ s"]


frames_per_step=2
frames_per_sec=2
pad_start_frames=2
pad_end_frames=2


import sys
sys.path[:0]=['/Change/This/Path']
from Permar_Functions_Python import *
import matplotlib
import matplotlib.pyplot as plt
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
from subprocess import call


connection_colors=[(171/255,210/255,229/255),(244/255,166/255,131/255),(102/255,102/255,102/255),(143/255,145/255,0/255),(1.0,128/255,0)]
connection_cmap=matplotlib.colors.ListedColormap(((171/255,210/255,229/255),(244/255,166/255,131/255),(102/255,102/255,102/255),(143/255,145/255,0/255),(1.0,128/255,0)))


call_result=call(["mkdir","./anim_temp"])

plot_idx=0
for idx_f in range(len(connectivity_filenames)):

	y_grid,x_grid,connectivity_map=load_connectivity_map(connectivity_filenames[idx_f])
	connectivity_map=np.floor(connectivity_map/2)

	y_spacing=y_grid[0,1]-y_grid[0,0]
	x_spacing=x_grid[1,0]-x_grid[0,0]

	plt.clf()
	fig1=plt.figure("Connectivity map",figsize=(10,10))
	ax1_1=fig1.gca()
	color_plot1_1=plt.pcolormesh(x_grid-0.5*x_spacing,y_grid-0.5*y_spacing,connectivity_map+0.5,cmap=connection_cmap,vmin=0,vmax=5,rasterized=True)
	plt.title(titles[idx_f],fontsize=46)
	plt.tick_params(axis='both', which='major',labelsize=40,direction='in',bottom=True, top=True, left=True, right=True)
	plt.ylabel(r"$y$ [a.u.]",fontsize=42)
	plt.xlabel(r"$x$ [a.u.]",fontsize=42)
	#plt.xlim(x_limits)
	#plt.ylim(y_limits)

	plt.savefig("./anim_temp/img{:03d}.png".format(plot_idx), format="png", dpi=100,bbox_inches='tight',pad_inches=0.1)
	plt.close(fig1)
	plot_idx+=1
	if idx_f==0:
		for idx in range(pad_start_frames):
			call_result=call(["cp","./anim_temp/img{:03d}.png".format(plot_idx-1),"./anim_temp/img{:03d}.png".format(plot_idx)])
			plot_idx+=1
	for idx in range(frames_per_step):
		call_result=call(["cp","./anim_temp/img{:03d}.png".format(plot_idx-1),"./anim_temp/img{:03d}.png".format(plot_idx)])
		plot_idx+=1
	if idx_f==len(connectivity_filenames)-1:
		for idx in range(pad_end_frames):
			call_result=call(["cp","./anim_temp/img{:03d}.png".format(plot_idx-1),"./anim_temp/img{:03d}.png".format(plot_idx)])
			plot_idx+=1


call_result=call(['ffmpeg -framerate '+str(frames_per_sec)+' -i ./anim_temp/img%03d.png -c:v libx264 -r 30 -pix_fmt yuv420p -vf "pad=ceil(iw/2)*2:ceil(ih/2)*2" ./anim_temp/anim.mp4'],shell=True)
call_result=call(["cp","./anim_temp/anim.mp4","./Animate_Connectivity_Map.mp4"])
call_result=call(["rm","-r","./anim_temp/"])







