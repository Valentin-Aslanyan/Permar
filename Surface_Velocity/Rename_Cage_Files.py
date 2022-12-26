

data_directory="Data"
new_name_stem="ConnectivityCage"


import os

cage_files = [f for f in os.listdir(data_directory) if os.path.isfile(os.path.join(data_directory, f)) and 'SurfaceCage' in f]
cage_files.sort()

for f in cage_files:
	os.rename(os.path.join(data_directory, f),os.path.join(data_directory, f.replace('SurfaceCage',new_name_stem)))

