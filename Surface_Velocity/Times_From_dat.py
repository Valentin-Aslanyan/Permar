

data_directory="./Data"
target_filename="lare3d.dat"


import os
import numpy as np

fortran_max_line=132
indices_raw=[]
times_raw=[]
infile=open(os.path.join(data_directory,target_filename),"r")
line=infile.readline()
while line:
	line=infile.readline()
	if 'dumping' in line.lower():
		indices_raw.append(line.split()[1])
		times_raw.append(line.split()[4])
infile.close()

indices_nooverlap=[int(indices_raw[-1])]
times_nooverlap=[float(times_raw[-1])]
for idx in range(len(indices_raw)-2,-1,-1):
	index_curr=int(indices_raw[idx])
	if index_curr<indices_nooverlap[-1]:
		indices_nooverlap.append(index_curr)
		times_nooverlap.append(float(times_raw[idx]))

indices_nooverlap=np.array(indices_nooverlap[::-1],dtype=np.int32)
times_nooverlap=np.array(times_nooverlap[::-1],dtype=np.float64)

sort_idx=np.argsort(times_nooverlap)

indices=indices_nooverlap[sort_idx]
times=times_nooverlap[sort_idx]

print("      start_time = "+str(times[0])+"_num")
print("")
print("Note, first element in the array below is also start_time, so delete if necessary")
print("")

line_curr="      output_times = (/ "
for idx in range(len(times)-1):
	str_curr=str(times[idx])+"_num,"
	if len(line_curr)+len(str_curr)<fortran_max_line-2:
		line_curr=line_curr+str_curr
	else:
		line_curr=line_curr+' &'
		print(line_curr)
		line_curr="      & "+str_curr
str_curr=str(times[-1])+"_num /)"
if len(line_curr)+len(str_curr)<fortran_max_line-2:
	line_curr=line_curr+str_curr
else:
	line_curr=line_curr+' &'
	print(line_curr)
	line_curr="      & "+str_curr
print(line_curr)



