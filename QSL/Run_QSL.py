
#Keep as integers
steps_to_process=[0]

binaries_to_run=["qslSquasher_Lare1"]
#Do not include .qsl extension
outfile_names=["qslLare"]

from os import listdir
from os.path import isfile, join
from subprocess import call

existing_folders = [f for f in listdir("./") if not isfile(f) and f.isnumeric()]
existing_folders.sort()

existing_files = [f for f in listdir("../../") if isfile(join("../../", f))]
existing_files.sort()

binaries_ready=True
for binary in binaries_to_run:
	if binary not in existing_files:
		print("Require "+binary)
		binaries_ready=False

if binaries_ready:
	for current_step in steps_to_process:
		step_folder="{:04d}".format(current_step)
		print(step_folder)
		if step_folder not in existing_folders:
			print("Error: Folder "+step_folder+" does not exist")
		else:
			for idx in range(len(binaries_to_run)):
				print(binaries_to_run[idx])
				call(["cp","../../"+binaries_to_run[idx],step_folder])
				call(["./"+binaries_to_run[idx]],cwd=step_folder)
				call(["mv","qsl.bin",outfile_names[idx]+".qsl"],cwd=step_folder)

