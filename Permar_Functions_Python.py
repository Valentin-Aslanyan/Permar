"""
Functions for reading and processing of .sdf files written by the Lare code
NOTE: these are different to more common types of .sdf file (e.g. Scientific Data Format) for which Python routines already exist
"""

import numpy as np
import struct
import os


class cfsasdf_file_properties:#TODO - None for unset properties?
	def __init__(self):
		self.endianness=0
		self.sdf_version=0
		self.sdf_revision=0
		self.code_name=""
		self.nblocks=0
		self.step=0
		self.time=0.0
		self.jobid1=0
		self.jobid2=0
		self.code_io_version=0
		self.restart_flag=0
		self.subdomain_file=0
		

class cfsasdf_block:#TODO: remaining block types, None for unset properties?
	def __init__(self,blocktype):
		self.blocktype=blocktype
		self.ndims=0
		self.block_id=""
		self.block_name=""
		self.data=0.0
		if blocktype==1:#c_blocktype_plain_mesh
			self.mults=0.0
			self.labels=[]
			self.units=[]
			self.geometry_type=0
			self.minval=0.0
			self.maxval=0.0
			self.dims=0
		elif blocktype==2:#c_blocktype_point_mesh
			self.mults=0.0
			self.labels=[]
			self.units=[]
			self.geometry_type=0
			self.minval=0.0
			self.maxval=0.0
			self.num_points=0
		elif blocktype==3:#c_blocktype_plain_variable
			self.mult=0.0
			self.units=""
			self.mesh_id=""
			self.dims=0
			self.stagger=0
		elif blocktype==4:#c_blocktype_point_variable
			self.mult=0.0
			self.units=""
			self.mesh_id=""
			self.num_points=0


class lare3d_dat_quantities:#TODO - None for unset properties?
	def __init__(self):
		self.nprocx=0
		self.nprocy=0
		self.nprocz=0
		self.nx=0
		self.ny=0
		self.nz=0
		self.length_x=0.0
		self.length_y=0.0
		self.length_z=0.0
		self.visc1=0.0
		self.visc2=0.0
		self.j_max=0.0
		self.eta0=0.0
		self.eta_background=0.0
		self.mass_fraction=0.0
		self.B_norm=0.0
		self.L_norm=0.0
		self.Rho_norm=0.0
		self.v_norm=0.0
		self.time_norm=0.0
		self.temp_norm=0.0


#If target_quantities==None or target_quantities=='all', return all, otherwise supply list
def read_cfsasdf_file(filename,target_quantities=None):#TODO: checks, ignore bad blocks, switch based on endianness, comments, dict?
	sdf_datatypes=[[1,'c','byte'],[4,'i','int32'],[8,'q','int64'],[4,'f','float32'],[8,'d','float64'],[16,'d','float64'],[1,'c','byte'],[1,'c','byte'],[1,'c','byte']]

	current_file=cfsasdf_file_properties()

	infile=open(filename,"rb")
	sdf1_str=infile.read(4)
	current_file.endianness=struct.unpack('i', infile.read(4))[0]
	current_file.sdf_version=struct.unpack('i', infile.read(4))[0]
	current_file.sdf_revision=struct.unpack('i', infile.read(4))[0]
	current_file.code_name=infile.read(32).decode('ascii')
	first_block_location=struct.unpack('q', infile.read(8))[0]
	next_block_location=first_block_location
	summary_location=struct.unpack('q', infile.read(8))[0]
	summary_size=struct.unpack('i', infile.read(4))[0]
	current_file.nblocks=struct.unpack('i', infile.read(4))[0]
	block_header_length=struct.unpack('i', infile.read(4))[0]
	current_file.step=struct.unpack('i', infile.read(4))[0]
	current_file.time=struct.unpack('d', infile.read(8))[0]
	current_file.jobid1=struct.unpack('i', infile.read(4))[0]
	current_file.jobid2=struct.unpack('i', infile.read(4))[0]
	string_length=struct.unpack('i', infile.read(4))[0]
	current_file.code_io_version=struct.unpack('i', infile.read(4))[0]
	current_file.restart_flag=struct.unpack('b', infile.read(1))[0]
	current_file.subdomain_file=struct.unpack('b', infile.read(1))[0]

	block_data=[]
	for idx_b in range(current_file.nblocks):
		infile.seek(next_block_location)
		current_block_location=next_block_location
		next_block_location=struct.unpack('q', infile.read(8))[0]
		data_location=struct.unpack('q', infile.read(8))[0]
		block_id=infile.read(32).decode('ascii')
		data_length=struct.unpack('q', infile.read(8))[0]
		blocktype=struct.unpack('i', infile.read(4))[0]
		datatype=struct.unpack('i', infile.read(4))[0]
		ndims=struct.unpack('i', infile.read(4))[0]
		block_name=infile.read(string_length).decode('ascii')
		block_info_length=struct.unpack('i', infile.read(4))[0]
		current_block=cfsasdf_block(blocktype)
		current_block.block_name=block_name
		current_block.ndims=ndims
		current_block.block_id=block_id
		if blocktype==1 or blocktype==2:
			infile.seek(current_block_location+block_header_length)
			current_block.mults=struct.unpack('d'*ndims, infile.read(8*ndims))
			for idx_l in range(ndims):
				current_block.labels.append(infile.read(32).decode('ascii'))
			for idx_l in range(ndims):
				current_block.units.append(infile.read(32).decode('ascii'))
			current_block.geometry_type=struct.unpack('i', infile.read(4))[0]
			current_block.minval=struct.unpack('d'*ndims, infile.read(8*ndims))
			current_block.maxval=struct.unpack('d'*ndims, infile.read(8*ndims))
			if blocktype==1:
				current_block.dims=struct.unpack('i'*ndims, infile.read(4*ndims))
			elif blocktype==2:
				current_block.num_points=struct.unpack('q', infile.read(8))[0]
		elif blocktype==3 or blocktype==4:
			infile.seek(current_block_location+block_header_length)
			current_block.mult=struct.unpack('d', infile.read(8))
			current_block.units=infile.read(32).decode('ascii')
			current_block.mesh_id=infile.read(32).decode('ascii')
			if blocktype==3:
				current_block.dims=struct.unpack('i'*ndims, infile.read(4*ndims))
				current_block.stagger=struct.unpack('i', infile.read(4))[0]
			if blocktype==4:
				current_block.num_points=struct.unpack('q', infile.read(8))[0]
		infile.seek(data_location)
		data_num_total=data_length//sdf_datatypes[datatype][0]

		if target_quantities==None:
			append_block=True
		elif type(target_quantities)==str:
			if target_quantities.lower()=='all' or target_quantities in current_block.block_id:
				append_block=True
			else:
				append_block=False
		elif hasattr(target_quantities, '__iter__'):
			append_block=False
			for elmt in target_quantities:
				if elmt in current_block.block_id:
					append_block=True
		else:
			append_block=False
		if append_block:
			current_block.data=np.array(struct.unpack(sdf_datatypes[datatype][1]*data_num_total, infile.read(data_length)),dtype=sdf_datatypes[datatype][2])
			block_data.append(current_block)
		else:
			infile.seek(data_length,1)
	infile.close()
	return current_file,block_data


def read_lare3d_dat(filename):

	current_data=lare3d_dat_quantities()

	infile=open(filename,"r")
	line=infile.readline()
	line=infile.readline()
	line=infile.readline().split()
	current_data.nprocx=int(line[4])
	current_data.nprocy=int(line[5])
	current_data.nprocz=int(line[6])
	line=infile.readline().split()
	current_data.nx=int(line[4])
	current_data.ny=int(line[5])
	current_data.nz=int(line[6])
	line=infile.readline()
	line=infile.readline().split()
	current_data.length_x=float(line[2])
	line=infile.readline().split()
	current_data.length_y=float(line[2])
	line=infile.readline().split()
	current_data.length_z=float(line[2])
	line=infile.readline()
	line=infile.readline()
	line=infile.readline().split()
	current_data.visc1=float(line[4])
	line=infile.readline().split()
	current_data.visc2=float(line[4])
	line=infile.readline().split()
	current_data.j_max=float(line[2])
	line=infile.readline().split()
	current_data.eta0=float(line[2])
	line=infile.readline().split()
	current_data.eta_background=float(line[2])
	line=infile.readline().split()
	current_data.kappa=float(line[2])
	line=infile.readline()
	line=infile.readline().split()
	current_data.mass_fraction=float(line[2])
	line=infile.readline().split()
	current_data.B_norm=float(line[3])
	line=infile.readline().split()
	current_data.L_norm=float(line[3])
	line=infile.readline().split()
	current_data.Rho_norm=float(line[3])
	line=infile.readline().split()
	current_data.v_norm=float(line[3])
	line=infile.readline().split()
	current_data.time_norm=float(line[3])
	line=infile.readline().split()
	current_data.temp_norm=float(line[3])
	infile.close()

	return current_data


def unstagger_reshape(grid_x,grid_y,grid_z,current_data,current_stagger):
	nx_actual=len(grid_x)
	ny_actual=len(grid_y)
	nz_actual=len(grid_z)
	if current_stagger & 1 == 0:
		nx_actual=nx_actual-1
		grid_x=0.5*(grid_x[1:]+grid_x[:-1])
	if current_stagger & 2 == 0:
		ny_actual=ny_actual-1
		grid_y=0.5*(grid_y[1:]+grid_y[:-1])
	if current_stagger & 4 == 0:
		nz_actual=nz_actual-1
		grid_z=0.5*(grid_z[1:]+grid_z[:-1])
	current_data=current_data.reshape((nz_actual,ny_actual,nx_actual))
	return grid_x,grid_y,grid_z,current_data


def quantity_from_blocks(block_data,target_block_id):
	for idx in range(len(block_data)):
		if "grid" in block_data[idx].block_id:
			grid_dims=block_data[idx].dims
			grid_x=block_data[idx].data[:grid_dims[0]]
			grid_y=block_data[idx].data[grid_dims[0]:grid_dims[0]+grid_dims[1]]
			grid_z=block_data[idx].data[grid_dims[0]+grid_dims[1]:]
		if target_block_id in block_data[idx].block_id:
			current_data=block_data[idx].data
			current_stagger=block_data[idx].stagger

	return unstagger_reshape(grid_x,grid_y,grid_z,current_data,current_stagger)


#Memory-efficient way to read a single quantity from SDF file; TODO: checks, reject bad blocks
def read_cfsasdf_file_lean(filename,target_quantity):
	sdf_datatypes=[[1,'c','byte'],[4,'i','int32'],[8,'q','int64'],[4,'f','float32'],[8,'d','float64'],[16,'d','float64'],[1,'c','byte'],[1,'c','byte'],[1,'c','byte']]

	infile=open(filename,"rb")

	sdf1_str=infile.read(4)
	endianness=struct.unpack('i', infile.read(4))[0]
	sdf_version=struct.unpack('i', infile.read(4))[0]
	sdf_revision=struct.unpack('i', infile.read(4))[0]
	code_name=infile.read(32).decode('ascii')
	first_block_location=struct.unpack('q', infile.read(8))[0]
	next_block_location=first_block_location
	summary_location=struct.unpack('q', infile.read(8))[0]
	summary_size=struct.unpack('i', infile.read(4))[0]
	nblocks=struct.unpack('i', infile.read(4))[0]
	block_header_length=struct.unpack('i', infile.read(4))[0]
	step=struct.unpack('i', infile.read(4))[0]
	time=struct.unpack('d', infile.read(8))[0]
	jobid1=struct.unpack('i', infile.read(4))[0]
	jobid2=struct.unpack('i', infile.read(4))[0]
	string_length=struct.unpack('i', infile.read(4))[0]
	code_io_version=struct.unpack('i', infile.read(4))[0]
	restart_flag=struct.unpack('b', infile.read(1))[0]
	subdomain_file=struct.unpack('b', infile.read(1))[0]

	for idx_b in range(nblocks):
		infile.seek(next_block_location)
		current_block_location=next_block_location
		next_block_location=struct.unpack('q', infile.read(8))[0]
		data_location=struct.unpack('q', infile.read(8))[0]
		block_id=infile.read(32).decode('ascii')
		if "grid" in block_id:
			data_length=struct.unpack('q', infile.read(8))[0]
			blocktype=struct.unpack('i', infile.read(4))[0]
			datatype=struct.unpack('i', infile.read(4))[0]
			ndims=struct.unpack('i', infile.read(4))[0]
			block_name=infile.read(string_length).decode('ascii')
			block_info_length=struct.unpack('i', infile.read(4))[0]
			infile.seek(current_block_location+block_header_length)
			mults=struct.unpack('d'*ndims, infile.read(8*ndims))
			infile.seek(current_block_location+block_header_length+88*ndims+4)
			dims=struct.unpack('i'*ndims, infile.read(4*ndims))
			infile.seek(data_location)
			grid_x=np.array(struct.unpack(sdf_datatypes[datatype][1]*dims[0],infile.read(sdf_datatypes[datatype][0]*dims[0])),dtype=sdf_datatypes[datatype][2])
			grid_y=np.array(struct.unpack(sdf_datatypes[datatype][1]*dims[1],infile.read(sdf_datatypes[datatype][0]*dims[1])),dtype=sdf_datatypes[datatype][2])
			grid_z=np.array(struct.unpack(sdf_datatypes[datatype][1]*dims[2],infile.read(sdf_datatypes[datatype][0]*dims[2])),dtype=sdf_datatypes[datatype][2])
		elif target_quantity in block_id:
			data_length=struct.unpack('q', infile.read(8))[0]
			blocktype=struct.unpack('i', infile.read(4))[0]
			datatype=struct.unpack('i', infile.read(4))[0]
			ndims=struct.unpack('i', infile.read(4))[0]
			block_name=infile.read(string_length).decode('ascii')
			block_info_length=struct.unpack('i', infile.read(4))[0]
			infile.seek(current_block_location+block_header_length)
			mults=struct.unpack('d', infile.read(8))
			infile.seek(current_block_location+block_header_length+72)
			dims=struct.unpack('i'*ndims, infile.read(4*ndims))
			target_array=np.zeros((dims[2],dims[1],dims[0]),dtype=sdf_datatypes[datatype][2])
			infile.seek(data_location)
			for idx_z in range(dims[2]):
				target_array[idx_z,:,:]=np.array(struct.unpack(sdf_datatypes[datatype][1]*dims[0]*dims[1],infile.read(sdf_datatypes[datatype][0]*dims[0]*dims[1])),dtype=sdf_datatypes[datatype][2]).reshape(dims[1],dims[0])
	infile.close()
	return time,grid_x,grid_y,grid_z,target_array


def regularize_Bfield(block_data): #TODO: currently assumes "ordinary" stagger, make use of _stagger in future
	grid_found=False
	Bx_found=False
	By_found=False
	Bz_found=False
	for idx in range(len(block_data)):
		if "grid" in block_data[idx].block_id:
			grid_dims=block_data[idx].dims
			grid_x=block_data[idx].data[:grid_dims[0]]
			grid_y=block_data[idx].data[grid_dims[0]:grid_dims[0]+grid_dims[1]]
			grid_z=block_data[idx].data[grid_dims[0]+grid_dims[1]:]
			grid_found=True
		if "Bx" in block_data[idx].block_id:
			Bx=block_data[idx].data
			Bx_stagger=block_data[idx].stagger
			Bx_found=True
		if "By" in block_data[idx].block_id:
			By=block_data[idx].data
			By_stagger=block_data[idx].stagger
			By_found=True
		if "Bz" in block_data[idx].block_id:
			Bz=block_data[idx].data
			Bz_stagger=block_data[idx].stagger
			Bz_found=True
	if grid_found and Bx_found and By_found and Bz_found:
		Bx=Bx.reshape((len(grid_z)-1,len(grid_y)-1,len(grid_x)))
		By=By.reshape((len(grid_z)-1,len(grid_y),len(grid_x)-1))
		Bz=Bz.reshape((len(grid_z),len(grid_y)-1,len(grid_x)-1))

		grid_x=0.5*(grid_x[1:]+grid_x[:-1])
		grid_y=0.5*(grid_y[1:]+grid_y[:-1])
		grid_z=0.5*(grid_z[1:]+grid_z[:-1])
		B_full=np.zeros((len(grid_z),len(grid_y),len(grid_x),3))
		B_full[:,:,:,0]=0.5*(Bx[:,:,1:]+Bx[:,:,:-1])
		B_full[:,:,:,1]=0.5*(By[:,1:,:]+By[:,:-1,:])
		B_full[:,:,:,2]=0.5*(Bz[1:,:,:]+Bz[:-1,:,:])
		return grid_x,grid_y,grid_z,B_full
	else:
		print("Warning, error in Bfield regularization")
		return np.zeros((1)),np.zeros((1)),np.zeros((1)),np.zeros((1,1,1,3))


def read_Lare_en_file(en_filename):#TODO: switch based on endianness, real size
	infile_size=os.path.getsize(en_filename)
	infile=open(en_filename,"rb")
	magic=str(infile.read(3))
	version=struct.unpack('i', infile.read(4))[0]
	revision=struct.unpack('i', infile.read(4))[0]
	endianness=struct.unpack('i', infile.read(4))[0]
	header_length=struct.unpack('i', infile.read(4))[0]
	num_sz_in=struct.unpack('i', infile.read(4))[0]
	en_nvars_in=struct.unpack('i', infile.read(4))[0]

	variable_names=infile.read(header_length-3-4*6).replace(b'\x00',b'').decode("utf-8").split()
	ndump=(infile_size-header_length)//(num_sz_in*en_nvars_in)
	time=[]
	variable_sum=[]
	for idx in range(ndump):
		time.append(struct.unpack('d', infile.read(8))[0])
		variable_sum.append(struct.unpack('d'*(en_nvars_in-1), infile.read(8*(en_nvars_in-1))))

	infile.close()
	time=np.array(time)
	variable_sum=np.array(variable_sum)

	return variable_names,time,variable_sum


def calculate_dt(block_data,gamma=1.4):#TODO WORK IN PROGRESS; boris,gamma_boris, alpha
	for idx in range(len(block_data)):
		if "grid" in block_data[idx].block_id:
			grid_dims=block_data[idx].dims
			grid_x=block_data[idx].data[:grid_dims[0]]
			grid_y=block_data[idx].data[grid_dims[0]:grid_dims[0]+grid_dims[1]]
			grid_z=block_data[idx].data[grid_dims[0]+grid_dims[1]:]
		if "Bx" in block_data[idx].block_id:
			Bx=block_data[idx].data
			Bx_stagger=block_data[idx].stagger
		if "By" in block_data[idx].block_id:
			By=block_data[idx].data
			By_stagger=block_data[idx].stagger
		if "Bz" in block_data[idx].block_id:
			Bz=block_data[idx].data
			Bz_stagger=block_data[idx].stagger
		if "Vx" in block_data[idx].block_id:
			Vx=block_data[idx].data
			Vx_stagger=block_data[idx].stagger
		if "Vy" in block_data[idx].block_id:
			Vy=block_data[idx].data
			Vy_stagger=block_data[idx].stagger
		if "Vz" in block_data[idx].block_id:
			Vz=block_data[idx].data
			Vz_stagger=block_data[idx].stagger
		if "Rho" in block_data[idx].block_id:
			Rho=block_data[idx].data
			Rho_stagger=block_data[idx].stagger
		if "Pressure" in block_data[idx].block_id:
			p=block_data[idx].data
			p_stagger=block_data[idx].stagger

	dxb=grid_x[1:]-grid_x[:-1]
	dyb=grid_y[1:]-grid_y[:-1]
	dzb=grid_z[1:]-grid_z[:-1]
	dyb_mg,dzb_mg,dxb_mg=np.meshgrid(dyb,dzb,dxb)
	length=np.minimum(dxb_mg,dyb_mg,dzb_mg)

	xtemp,ytemp,ztemp,Bx=unstagger_reshape(grid_x,grid_y,grid_z,Bx,Bx_stagger)
	xtemp,ytemp,ztemp,By=unstagger_reshape(grid_x,grid_y,grid_z,By,By_stagger)
	xtemp,ytemp,ztemp,Bz=unstagger_reshape(grid_x,grid_y,grid_z,Bz,Bz_stagger)
	w1=Bx[:,:,:-1]**2+By[:,:-1,:]**2+Bz[:-1,:,:]**2

	xtemp,ytemp,ztemp,Rho=unstagger_reshape(grid_x,grid_y,grid_z,Rho,Rho_stagger)
	xtemp,ytemp,ztemp,p=unstagger_reshape(grid_x,grid_y,grid_z,p,p_stagger)
	rho0=np.clip(Rho,2.2250738585072014E-308,1E200)
	cs2=(gamma*p+w1)/rho0

	dt=length/np.sqrt(cs2)
	dt_cfl_idx=np.argmin(dt)
	print("min(Rho) =",min(Rho.flatten()))
	print(" dt(CFL) =",dt.flatten()[dt_cfl_idx])
	print()
	#print(dt.flatten()[dt_cfl_idx],length.flatten()[dt_cfl_idx],cs2.flatten()[dt_cfl_idx],p.flatten()[dt_cfl_idx],w1.flatten()[dt_cfl_idx],rho0.flatten()[dt_cfl_idx])


def calculate_heatdt(block_data,gamma=5/3,dt_parab=1.0E10,kappa_0=1.0,power=2.5):#TODO nonideal
	for idx in range(len(block_data)):
		if "grid" in block_data[idx].block_id:
			grid_dims=block_data[idx].dims
			grid_x=block_data[idx].data[:grid_dims[0]]
			grid_y=block_data[idx].data[grid_dims[0]:grid_dims[0]+grid_dims[1]]
			grid_z=block_data[idx].data[grid_dims[0]+grid_dims[1]:]
		if "Temperature" in block_data[idx].block_id:
			temp=block_data[idx].data
			temp_stagger=block_data[idx].stagger
		if "Rho" in block_data[idx].block_id:
			Rho=block_data[idx].data
			Rho_stagger=block_data[idx].stagger

	dxb=grid_x[1:]-grid_x[:-1]
	dyb=grid_y[1:]-grid_y[:-1]
	dzb=grid_z[1:]-grid_z[:-1]
	dyb_mg,dzb_mg,dxb_mg=np.meshgrid(dyb,dzb,dxb)

	xtemp,ytemp,ztemp,Rho=unstagger_reshape(grid_x,grid_y,grid_z,Rho,Rho_stagger)
	xtemp,ytemp,ztemp,temp=unstagger_reshape(grid_x,grid_y,grid_z,temp,temp_stagger)
	temp2=2.0*kappa_0*temp**power
	
	dt1=Rho*dxb_mg**2/temp2
	dt2=Rho*dyb_mg**2/temp2
	dt3=Rho*dzb_mg**2/temp2

	dt1_idx=np.argmin(dt1)
	dt2_idx=np.argmin(dt2)
	dt3_idx=np.argmin(dt3)
	print(dt1.flatten()[dt1_idx],dt2.flatten()[dt2_idx],dt3.flatten()[dt3_idx])
	print(dzb_mg.flatten()[dt1_idx],dzb_mg.flatten()[dt2_idx],dzb_mg.flatten()[dt3_idx])


def Trilinear_interpolation(x,y,z,x0,x1,y0,y1,z0,z1,f000,f001,f010,f011,f100,f101,f110,f111):
	"""
	Must be a cuboid
	f000=f(x0,y0,z0)
	f100=f(x1,y0,z0) etc
	"""
	x_shift=x-x0
	y_shift=y-y0
	z_shift=z-z0
	delta_x=x1-x0
	delta_y=y1-y0
	delta_z=z1-z0
	a1=(f100-f000)/delta_x
	a2=(f010-f000)/delta_y
	a3=(f001-f000)/delta_z
	a4=(f110-f100-f010+f000)/delta_y/delta_x
	a5=(f110-f100-f001+f000)/delta_z/delta_x
	a6=(f011-f001-f010+f000)/delta_z/delta_y
	a7=(f111-f110-f101-f011+f001+f010+f100-f000)/delta_x/delta_y/delta_z
	return f000+a1*x_shift+a2*y_shift+a3*z_shift+a4*x_shift*y_shift+a5*x_shift*z_shift+a6*y_shift*z_shift+a7*x_shift*y_shift*z_shift


def get_grid_index(target_pos,grid,last_index):
	last_index=min(max(last_index,1),len(grid)-3)
	if grid[last_index]<=target_pos and grid[last_index+1]>=target_pos:
		return last_index
	elif grid[last_index]>target_pos:
		idx=last_index-1
		while idx>0 and grid[idx]>target_pos:
			idx-=1
		return idx
	else:
		idx=last_index+1
		while idx<len(grid)-2 and grid[idx+1]<target_pos:
			idx+=1
		return idx


#Has all components of B in same grid, B[i,j,k,0]=Bx etc; pos_in has 3 coordinates of start point
def trace_fieldline(pos_in,X,Y,Z,B,max_steps=10000,dt=0.01,periodic_X=False,periodic_Y=False,periodic_Z=False):
	r_0=np.array(pos_in)
	delta_X=X[-1]-X[0]
	delta_Y=Y[-1]-Y[0]
	delta_Z=Z[-1]-Z[0]

	idx_x=get_grid_index(r_0[0],X,0)
	idx_y=get_grid_index(r_0[1],Y,0)
	idx_z=get_grid_index(r_0[2],Z,0)
	
	n_forward=1
	r_forward=np.zeros((max_steps,3))
	r_forward[0,:]=r_0
	in_bounds=True
	while n_forward<max_steps and in_bounds:
		idx_x=get_grid_index(r_0[0],X,idx_x)
		idx_y=get_grid_index(r_0[1],Y,idx_y)
		idx_z=get_grid_index(r_0[2],Z,idx_z)
		B_1=Trilinear_interpolation(r_0[0],r_0[1],r_0[2],X[idx_x],X[idx_x+1],Y[idx_y],Y[idx_y+1],Z[idx_z],Z[idx_z+1],B[idx_z,idx_y,idx_x,:], B[idx_z+1,idx_y,idx_x,:], B[idx_z,idx_y+1,idx_x,:], B[idx_z+1,idx_y+1,idx_x,:], B[idx_z,idx_y,idx_x+1,:], B[idx_z+1,idx_y,idx_x+1,:], B[idx_z,idx_y+1,idx_x+1,:], B[idx_z+1,idx_y+1,idx_x+1,:])
		k_1=dt*B_1/np.linalg.norm(B_1)

		r_1=r_0+k_1*0.5
		idx_x=get_grid_index(r_1[0],X,idx_x)
		idx_y=get_grid_index(r_1[1],Y,idx_y)
		idx_z=get_grid_index(r_1[2],Z,idx_z)
		B_2=Trilinear_interpolation(r_1[0],r_1[1],r_1[2],X[idx_x],X[idx_x+1],Y[idx_y],Y[idx_y+1],Z[idx_z],Z[idx_z+1],B[idx_z,idx_y,idx_x,:], B[idx_z+1,idx_y,idx_x,:], B[idx_z,idx_y+1,idx_x,:], B[idx_z+1,idx_y+1,idx_x,:], B[idx_z,idx_y,idx_x+1,:], B[idx_z+1,idx_y,idx_x+1,:], B[idx_z,idx_y+1,idx_x+1,:], B[idx_z+1,idx_y+1,idx_x+1,:])
		k_2=dt*B_2/np.linalg.norm(B_2)
		
		r_2=r_0+k_2*0.5
		idx_x=get_grid_index(r_2[0],X,idx_x)
		idx_y=get_grid_index(r_2[1],Y,idx_y)
		idx_z=get_grid_index(r_2[2],Z,idx_z)
		B_3=Trilinear_interpolation(r_2[0],r_2[1],r_2[2],X[idx_x],X[idx_x+1],Y[idx_y],Y[idx_y+1],Z[idx_z],Z[idx_z+1],B[idx_z,idx_y,idx_x,:], B[idx_z+1,idx_y,idx_x,:], B[idx_z,idx_y+1,idx_x,:], B[idx_z+1,idx_y+1,idx_x,:], B[idx_z,idx_y,idx_x+1,:], B[idx_z+1,idx_y,idx_x+1,:], B[idx_z,idx_y+1,idx_x+1,:], B[idx_z+1,idx_y+1,idx_x+1,:])
		k_3=dt*B_3/np.linalg.norm(B_3)

		r_3=r_0+k_3
		idx_x=get_grid_index(r_3[0],X,idx_x)
		idx_y=get_grid_index(r_3[1],Y,idx_y)
		idx_z=get_grid_index(r_3[2],Z,idx_z)
		B_4=Trilinear_interpolation(r_3[0],r_3[1],r_3[2],X[idx_x],X[idx_x+1],Y[idx_y],Y[idx_y+1],Z[idx_z],Z[idx_z+1],B[idx_z,idx_y,idx_x,:], B[idx_z+1,idx_y,idx_x,:], B[idx_z,idx_y+1,idx_x,:], B[idx_z+1,idx_y+1,idx_x,:], B[idx_z,idx_y,idx_x+1,:], B[idx_z+1,idx_y,idx_x+1,:], B[idx_z,idx_y+1,idx_x+1,:], B[idx_z+1,idx_y+1,idx_x+1,:])
		k_4=dt*B_4/np.linalg.norm(B_4)

		r_0=r_0+(k_1+2.0*k_2+2.0*k_3+k_4)/6.0
		r_forward[n_forward,:]=r_0
		n_forward+=1
		if periodic_X:
			r_0[0]=(r_0[0]-X[0]) % delta_X + X[0]
		elif r_0[0]<X[0] or r_0[0]>X[-1]:
			in_bounds=False
		if periodic_Y:
			r_0[1]=(r_0[1]-Y[0]) % delta_Y + Y[0]
		elif r_0[1]<Y[0] or r_0[1]>Y[-1]:
			in_bounds=False
		if periodic_Z:
			r_0[2]=(r_0[2]-Z[0]) % delta_Z + Z[0]
		elif r_0[2]<Z[0] or r_0[2]>Z[-1]:
			in_bounds=False
	r_forward=r_forward[:n_forward]

	r_0=np.array(pos_in)
	n_back=0
	r_back=np.zeros((max_steps-1,3))
	in_bounds=True
	while n_back<max_steps-1 and r_0[0]>=X[0] and r_0[0]<=X[-1] and r_0[1]>=Y[0] and r_0[1]<=Y[-1] and r_0[2]>=Z[0] and r_0[2]<=Z[-1]:
		idx_x=get_grid_index(r_0[0],X,idx_x)
		idx_y=get_grid_index(r_0[1],Y,idx_y)
		idx_z=get_grid_index(r_0[2],Z,idx_z)
		B_1=Trilinear_interpolation(r_0[0],r_0[1],r_0[2],X[idx_x],X[idx_x+1],Y[idx_y],Y[idx_y+1],Z[idx_z],Z[idx_z+1],B[idx_z,idx_y,idx_x,:], B[idx_z+1,idx_y,idx_x,:], B[idx_z,idx_y+1,idx_x,:], B[idx_z+1,idx_y+1,idx_x,:], B[idx_z,idx_y,idx_x+1,:], B[idx_z+1,idx_y,idx_x+1,:], B[idx_z,idx_y+1,idx_x+1,:], B[idx_z+1,idx_y+1,idx_x+1,:])
		k_1=-dt*B_1/np.linalg.norm(B_1)

		r_1=r_0+k_1*0.5
		idx_x=get_grid_index(r_1[0],X,idx_x)
		idx_y=get_grid_index(r_1[1],Y,idx_y)
		idx_z=get_grid_index(r_1[2],Z,idx_z)
		B_2=Trilinear_interpolation(r_1[0],r_1[1],r_1[2],X[idx_x],X[idx_x+1],Y[idx_y],Y[idx_y+1],Z[idx_z],Z[idx_z+1],B[idx_z,idx_y,idx_x,:], B[idx_z+1,idx_y,idx_x,:], B[idx_z,idx_y+1,idx_x,:], B[idx_z+1,idx_y+1,idx_x,:], B[idx_z,idx_y,idx_x+1,:], B[idx_z+1,idx_y,idx_x+1,:], B[idx_z,idx_y+1,idx_x+1,:], B[idx_z+1,idx_y+1,idx_x+1,:])
		k_2=-dt*B_2/np.linalg.norm(B_2)
		
		r_2=r_0+k_2*0.5
		idx_x=get_grid_index(r_2[0],X,idx_x)
		idx_y=get_grid_index(r_2[1],Y,idx_y)
		idx_z=get_grid_index(r_2[2],Z,idx_z)
		B_3=Trilinear_interpolation(r_2[0],r_2[1],r_2[2],X[idx_x],X[idx_x+1],Y[idx_y],Y[idx_y+1],Z[idx_z],Z[idx_z+1],B[idx_z,idx_y,idx_x,:], B[idx_z+1,idx_y,idx_x,:], B[idx_z,idx_y+1,idx_x,:], B[idx_z+1,idx_y+1,idx_x,:], B[idx_z,idx_y,idx_x+1,:], B[idx_z+1,idx_y,idx_x+1,:], B[idx_z,idx_y+1,idx_x+1,:], B[idx_z+1,idx_y+1,idx_x+1,:])
		k_3=-dt*B_3/np.linalg.norm(B_3)

		r_3=r_0+k_3
		idx_x=get_grid_index(r_3[0],X,idx_x)
		idx_y=get_grid_index(r_3[1],Y,idx_y)
		idx_z=get_grid_index(r_3[2],Z,idx_z)
		B_4=Trilinear_interpolation(r_3[0],r_3[1],r_3[2],X[idx_x],X[idx_x+1],Y[idx_y],Y[idx_y+1],Z[idx_z],Z[idx_z+1],B[idx_z,idx_y,idx_x,:], B[idx_z+1,idx_y,idx_x,:], B[idx_z,idx_y+1,idx_x,:], B[idx_z+1,idx_y+1,idx_x,:], B[idx_z,idx_y,idx_x+1,:], B[idx_z+1,idx_y,idx_x+1,:], B[idx_z,idx_y+1,idx_x+1,:], B[idx_z+1,idx_y+1,idx_x+1,:])
		k_4=-dt*B_4/np.linalg.norm(B_4)

		r_0=r_0+(k_1+2.0*k_2+2.0*k_3+k_4)/6.0
		r_back[n_back,:]=r_0
		n_back+=1
		if periodic_X:
			r_0[0]=(r_0[0]-X[0]) % delta_X + X[0]
		elif r_0[0]<X[0] or r_0[0]>X[-1]:
			in_bounds=False
		if periodic_Y:
			r_0[1]=(r_0[1]-Y[0]) % delta_Y + Y[0]
		elif r_0[1]<Y[0] or r_0[1]>Y[-1]:
			in_bounds=False
		if periodic_Z:
			r_0[2]=(r_0[2]-Z[0]) % delta_Z + Z[0]
		elif r_0[2]<Z[0] or r_0[2]>Z[-1]:
			in_bounds=False
	r_back=r_back[:n_back]

	return np.concatenate((r_back[::-1],r_forward))


def split_fieldline(fieldline,X,Y,Z):
	delta_X=X[-1]-X[0]
	delta_Y=Y[-1]-Y[0]
	delta_Z=Z[-1]-Z[0]
	fieldine_list=[]
	
	for idx in range(len(fieldline[:,0])):
		...
	
	return fieldine_list


def shift_fieldline(fieldline,X,Y,Z):
	fieldline_new=np.copy(fieldline)
	delta_X=X[-1]-X[0]
	delta_Y=Y[-1]-Y[0]
	delta_Z=Z[-1]-Z[0]
	
	shift_X=0.0
	shift_Y=0.0
	shift_Z=0.0
	for idx in range(1,len(fieldline_new[:,0])):
		if abs(fieldline_new[idx,0]+shift_X-fieldline_new[idx-1,0])>0.9*delta_X:
			shift_X-=np.sign(fieldline_new[idx,0]+shift_X-fieldline_new[idx-1,0])*delta_X
		fieldline_new[idx,0]+=shift_X
		if abs(fieldline_new[idx,1]+shift_Y-fieldline_new[idx-1,1])>0.9*delta_Y:
			shift_Y-=np.sign(fieldline_new[idx,1]+shift_Y-fieldline_new[idx-1,1])*delta_Y
		fieldline_new[idx,1]+=shift_Y
		if abs(fieldline_new[idx,2]+shift_Z-fieldline_new[idx-1,2])>0.9*delta_Z:
			shift_Z-=np.sign(fieldline_new[idx,2]+shift_Z-fieldline_new[idx-1,2])*delta_Z
		fieldline_new[idx,2]+=shift_Z
	
	return fieldline_new


def parse_QSL_Larebinfile(filename):
	"""
	Load in QSL output file
	Note - uses special version of QSLsquasher adapted to output binary data
	"""
	filesize=os.path.getsize(filename)
	num_points=filesize//(28)
	infile=open(filename,"rb")
	x=np.zeros((num_points))
	y=np.zeros((num_points))
	z=np.zeros((num_points))
	Q=np.zeros((num_points))
	for idx in range(num_points):
		data=struct.unpack('dddf',infile.read(28))
		z[idx]=data[2]
		y[idx]=data[1]
		x[idx]=data[0]
		Q[idx]=data[3]
	infile.close()

	y_dim=1
	idx=1
	while idx<num_points:
		if x[idx]==x[0]:
			y_dim+=1
			idx+=1
		else:
			idx=num_points
	x_dim=num_points//y_dim


	z_actual=z[0]
	y=y.reshape((x_dim,y_dim))
	x=x.reshape((x_dim,y_dim))
	Q=Q.reshape((x_dim,y_dim))
	return z_actual,y,x,Q


def parse_QSL_LarebinfileExpansion(filename):
	"""
	Load in QSL output file with expansion factor calculation
	Note - uses special version of QSLsquasher adapted to output binary data
	"""
	filesize=os.path.getsize(filename)
	num_points=filesize//(44)
	infile=open(filename,"rb")
	z=np.zeros((num_points))
	y=np.zeros((num_points))
	x=np.zeros((num_points))
	Q=np.zeros((num_points))
	Bs=np.zeros((num_points))
	Be=np.zeros((num_points))
	for idx in range(num_points):
		data=struct.unpack('=dddfdd',infile.read(44))
		z[idx]=data[2]
		y[idx]=data[1]
		x[idx]=data[0]
		Q[idx]=data[3]
		Bs[idx]=data[4]
		Be[idx]=data[5]
	infile.close()

	y_dim=1
	idx=1
	while idx<num_points:
		if x[idx]==x[0]:
			y_dim+=1
			idx+=1
		else:
			idx=num_points
	x_dim=num_points//y_dim


	z_actual=z[0]
	y=y.reshape((x_dim,y_dim))
	x=x.reshape((x_dim,y_dim))
	Q=Q.reshape((x_dim,y_dim))
	Bs=Bs.reshape((x_dim,y_dim))
	Be=Be.reshape((x_dim,y_dim))
	return z_actual,y,x,Q,Bs,Be


#####
#Reading parallel surface velocity files
def read_surface_coordinates(filename):
	data_type="d"
	data_size=8

	filesize=os.path.getsize(filename)
	infile=open(filename,"rb")
	raw_data=infile.read(filesize)
	infile.close()
	return np.array(struct.unpack(data_type*(filesize//data_size),raw_data))


#timestep_label must be correct string
def get_split_surface_velocity_coords(data_dir,timestep_label):
	x_files = [f for f in os.listdir(data_dir) if os.path.isfile(os.path.join(data_dir, f)) and 'Surfx_' in f and f.count('_')==2]
	y_files = [f for f in os.listdir(data_dir) if os.path.isfile(os.path.join(data_dir, f)) and 'Surfy_' in f and f.count('_')==2]

	target_x_files=[]
	for idx in range(len(x_files)):
		if x_files[idx].replace('Surfx_','').replace('.bin','').split('_')[1]==timestep_label:
			target_x_files.append(x_files[idx])
	target_x_files.sort()
	
	target_y_files=[]
	for idx in range(len(y_files)):
		if y_files[idx].replace('Surfy_','').replace('.bin','').split('_')[1]==timestep_label:
			target_y_files.append(y_files[idx])
	target_y_files.sort()

	if len(target_x_files)!=len(target_y_files):
		print("Mismatch reading coordinate files")

	x_coords=[read_surface_coordinates(os.path.join(data_dir,fname)) for fname in target_x_files]
	x_unique=np.unique(np.concatenate(x_coords))
	y_coords=[read_surface_coordinates(os.path.join(data_dir,fname)) for fname in target_y_files]
	y_unique=np.unique(np.concatenate(y_coords))

	processors=[]
	x_indices=[]
	y_indices=[]
	for idx in range(len(target_x_files)):
		processors.append(target_x_files[idx].replace('Surfx_','').replace('.bin','').split('_')[0])
		x_indices.append([np.where(x_coords[idx][0]==x_unique)[0][0],np.where(x_coords[idx][-1]==x_unique)[0][0]])
		y_indices.append([np.where(y_coords[idx][0]==y_unique)[0][0],np.where(y_coords[idx][-1]==y_unique)[0][0]])

	return x_unique,y_unique,processors,np.array(x_indices,dtype=np.int32),np.array(y_indices,dtype=np.int32)


#processors and timestep_label must be correct strings
#dimension must be "x" or "y"
def get_split_surface_velocity(data_dir,dimension,processors,timestep_label,full_dims,x_indices,y_indices,ignore_ghost_cells=True,ng=2,target_steps=None):
	data_type="d"
	data_size=8

	chunk_size=(y_indices[0,1]-y_indices[0,0]+1)*(x_indices[0,1]-x_indices[0,0]+1)*data_size
	filename=os.path.join(data_dir, 'SurfV'+dimension+'_'+processors[0]+"_"+timestep_label+".bin")
	filesize=os.path.getsize(filename)
	num_steps=filesize//chunk_size
	for idx in range(1,len(processors)):
		chunk_size=(y_indices[idx,1]-y_indices[idx,0]+1)*(x_indices[idx,1]-x_indices[idx,0]+1)*data_size
		filename=os.path.join(data_dir, 'SurfV'+dimension+'_'+processors[idx]+"_"+timestep_label+".bin")
		filesize=os.path.getsize(filename)
		num_steps=min(num_steps,filesize//chunk_size)

	if target_steps==None:
		num_outputs=num_steps
	elif hasattr(target_steps, '__iter__'):
		target_steps.sort()
		num_outputs=0
		for quant in target_steps:
			if quant<=num_steps: num_outputs+=1
	else:
		target_steps=list(target_steps)
		num_outputs=1

	if ignore_ghost_cells:
		velocity=np.zeros((full_dims[1]-2*ng,full_dims[0]-2*ng,num_outputs))
	else:
		velocity=np.zeros((full_dims[1],full_dims[0],num_outputs))

	for idx in range(len(processors)):
		chunk_size=(y_indices[idx,1]-y_indices[idx,0]+1)*(x_indices[idx,1]-x_indices[idx,0]+1)*data_size
		filename=os.path.join(data_dir, 'SurfV'+dimension+'_'+processors[idx]+"_"+timestep_label+".bin")
		infile=open(filename,"rb")
		if num_outputs==num_steps:
			for idx_s in range(num_steps):
				raw_data=infile.read(chunk_size)
				raw_array=np.array(struct.unpack(data_type*(chunk_size//data_size),raw_data)).reshape((y_indices[idx,1]-y_indices[idx,0]+1,x_indices[idx,1]-x_indices[idx,0]+1))
				if ignore_ghost_cells:
					velocity[y_indices[idx,0]:y_indices[idx,1]+1-2*ng,x_indices[idx,0]:x_indices[idx,1]+1-2*ng,idx_s]=raw_array[ng:-ng,ng:-ng]
				else:
					velocity[y_indices[idx,0]:y_indices[idx,1]+1,x_indices[idx,0]:x_indices[idx,1]+1,idx_s]=raw_array
		elif num_outputs==1:
			infile.seek(chunk_size*target_steps[0],0)
			raw_data=infile.read(chunk_size)
			raw_array=np.array(struct.unpack(data_type*(chunk_size//data_size),raw_data)).reshape((y_indices[idx,1]-y_indices[idx,0]+1,x_indices[idx,1]-x_indices[idx,0]+1))
			if ignore_ghost_cells:
				velocity[y_indices[idx,0]:y_indices[idx,1]+1-2*ng,x_indices[idx,0]:x_indices[idx,1]+1-2*ng,0]=raw_array[ng:-ng,ng:-ng]
			else:
				velocity[y_indices[idx,0]:y_indices[idx,1]+1,x_indices[idx,0]:x_indices[idx,1]+1,0]=raw_array
		else:
			for idx_s in range(len(target_steps)):
				infile.seek(chunk_size*target_steps[idx_s],0)
				raw_data=infile.read(chunk_size)
				raw_array=np.array(struct.unpack(data_type*(chunk_size//data_size),raw_data)).reshape((y_indices[idx,1]-y_indices[idx,0]+1,x_indices[idx,1]-x_indices[idx,0]+1))
				if ignore_ghost_cells:
					velocity[y_indices[idx,0]:y_indices[idx,1]+1-2*ng,x_indices[idx,0]:x_indices[idx,1]+1-2*ng,idx_s]=raw_array[ng:-ng,ng:-ng]
				else:
					velocity[y_indices[idx,0]:y_indices[idx,1]+1,x_indices[idx,0]:x_indices[idx,1]+1,idx_s]=raw_array
		infile.close()

	return velocity


def merge_split_surface_velocity_coords_to_file(data_dir,ignore_ghost_cells=True,ng=2):
	coordinate_timesteps,vx_timesteps,vy_timesteps,vz_timesteps=list_split_surface_velocity_files(data_dir)
	for idx in range(len(coordinate_timesteps)):
		x_grid,y_grid,processors,x_indices,y_indices=get_split_surface_velocity_coords(data_dir,coordinate_timesteps[idx])
		outfilename=os.path.join(data_dir, 'Surfx.bin')
		outfile=open(outfilename,"wb")
		if ignore_ghost_cells:
			outfile.write(bytes(x_grid[ng:-ng]))
		else:
			outfile.write(bytes(x_grid))
		outfile.close()
		outfilename=os.path.join(data_dir, 'Surfy.bin')
		outfile=open(outfilename,"wb")
		if ignore_ghost_cells:
			outfile.write(bytes(y_grid[ng:-ng]))
		else:
			outfile.write(bytes(y_grid))
		outfile.close()


#processors and timestep_label must be correct strings
#dimension must be "x" or "y"
def merge_split_surface_velocity_to_file(data_dir,dimension,processors,timestep_label,full_dims,x_indices,y_indices,make_new,ignore_ghost_cells=True,ng=2):
	data_type="d"
	data_size=8

	if ignore_ghost_cells:
		velocity=np.zeros((full_dims[1]-2*ng,full_dims[0]-2*ng))
	else:
		velocity=np.zeros((full_dims[1],full_dims[0]))

	if not os.path.exists(os.path.join(data_dir, 'Surft.bin')):
		make_new=True

	infilename=os.path.join(data_dir, 'Surft_'+timestep_label+".bin")
	filesize=os.path.getsize(infilename)
	num_steps=filesize//data_size
	for idx in range(len(processors)):
		chunk_size=(y_indices[idx,1]-y_indices[idx,0]+1)*(x_indices[idx,1]-x_indices[idx,0]+1)*data_size
		infilename=os.path.join(data_dir, 'SurfV'+dimension+'_'+processors[idx]+"_"+timestep_label+".bin")
		filesize=os.path.getsize(infilename)
		num_steps=min(num_steps,filesize//chunk_size)

	if make_new:
		idx_start=0
		outfilename=os.path.join(data_dir, 'SurfV'+dimension+".bin")
		outfile1=open(outfilename,"wb")
		outfile2=open(os.path.join(data_dir, 'Surft.bin'),"wb")
	else: #append to existing files
		outfile2=open(os.path.join(data_dir, 'Surft.bin'),"rb")
		outfile2.seek(-data_size,2)
		final_time=struct.unpack(data_type,outfile2.read(data_size))
		outfile2.close()
		
		infile2=open(os.path.join(data_dir, 'Surft_'+timestep_label+'.bin'),"rb")
		for idx_start in range(num_steps):
			if struct.unpack(data_type,infile2.read(data_size))>final_time:
				break
		infile2.close()
		outfilename=os.path.join(data_dir, 'SurfV'+dimension+".bin")
		outfile1=open(outfilename,"ab")
		outfile2=open(os.path.join(data_dir, 'Surft.bin'),"ab")

	infile2=open(os.path.join(data_dir, 'Surft_'+timestep_label+'.bin'),"rb")
	infile2.seek(data_size*idx_start,0)
	for idx_s in range(idx_start,num_steps):
		outfile2.write(infile2.read(data_size))
		for idx in range(len(processors)):
			chunk_size=(y_indices[idx,1]-y_indices[idx,0]+1)*(x_indices[idx,1]-x_indices[idx,0]+1)*data_size
			infilename=os.path.join(data_dir, 'SurfV'+dimension+'_'+processors[idx]+"_"+timestep_label+".bin")
			infile1=open(infilename,"rb")
			infile1.seek(chunk_size*idx_s,0)

			raw_data=infile1.read(chunk_size)
			raw_array=np.array(struct.unpack(data_type*(chunk_size//data_size),raw_data)).reshape((y_indices[idx,1]-y_indices[idx,0]+1,x_indices[idx,1]-x_indices[idx,0]+1))
			if ignore_ghost_cells:
				velocity[y_indices[idx,0]:y_indices[idx,1]+1-2*ng,x_indices[idx,0]:x_indices[idx,1]+1-2*ng]=raw_array[ng:-ng,ng:-ng]
			else:
				velocity[y_indices[idx,0]:y_indices[idx,1]+1,x_indices[idx,0]:x_indices[idx,1]+1]=raw_array
			infile1.close()

		for idx_y in range(len(velocity[:,0])):
			outfile1.write(bytes(velocity[idx_y,:]))
	infile2.close()
	outfile1.close()
	outfile2.close()


def list_split_surface_velocity_files(data_dir):
	x_files = [f for f in os.listdir(data_dir) if os.path.isfile(os.path.join(data_dir, f)) and 'Surfx_' in f and f.count('_')==2]
	y_files = [f for f in os.listdir(data_dir) if os.path.isfile(os.path.join(data_dir, f)) and 'Surfy_' in f and f.count('_')==2]

	coordinate_timesteps_x=[]
	for idx in range(len(x_files)):
		if x_files[idx].replace('Surfx_','').replace('.bin','').split('_')[1] not in coordinate_timesteps_x:
			coordinate_timesteps_x.append(x_files[idx].replace('Surfx_','').replace('.bin','').split('_')[1])

	y_filelabels=np.zeros((len(y_files),2),dtype=np.int32)
	coordinate_timesteps_y=[]
	for idx in range(len(y_files)):
		if y_files[idx].replace('Surfy_','').replace('.bin','').split('_')[1] not in coordinate_timesteps_y:
			coordinate_timesteps_y.append(y_files[idx].replace('Surfy_','').replace('.bin','').split('_')[1])
		y_filelabels[idx,:]=y_files[idx].replace('Surfy_','').replace('.bin','').split('_')
	coordinate_timesteps=[]
	for idx in range(len(coordinate_timesteps_x)):
		if coordinate_timesteps_x[idx] in coordinate_timesteps_y:
			coordinate_timesteps.append(coordinate_timesteps_x[idx])

	vx_files = [f for f in os.listdir(data_dir) if os.path.isfile(os.path.join(data_dir, f)) and 'SurfVx_' in f and f.count('_')==2]
	vx_timesteps=[]
	for idx in range(len(vx_files)):
		if vx_files[idx].replace('SurfVx_','').replace('.bin','').split('_')[1] not in vx_timesteps:
			vx_timesteps.append(vx_files[idx].replace('SurfVx_','').replace('.bin','').split('_')[1])


	vy_files = [f for f in os.listdir(data_dir) if os.path.isfile(os.path.join(data_dir, f)) and 'SurfVy_' in f and f.count('_')==2]
	vy_timesteps=[]
	for idx in range(len(vy_files)):
		if vy_files[idx].replace('SurfVy_','').replace('.bin','').split('_')[1] not in vy_timesteps:
			vy_timesteps.append(vy_files[idx].replace('SurfVy_','').replace('.bin','').split('_')[1])

	vz_files = [f for f in os.listdir(data_dir) if os.path.isfile(os.path.join(data_dir, f)) and 'SurfVz_' in f and f.count('_')==2]
	vz_timesteps=[]
	for idx in range(len(vz_files)):
		if vz_files[idx].replace('SurfVz_','').replace('.bin','').split('_')[1] not in vz_timesteps:
			vz_timesteps.append(vz_files[idx].replace('SurfVz_','').replace('.bin','').split('_')[1])

	coordinate_timesteps.sort()
	vx_timesteps.sort()
	vy_timesteps.sort()
	vz_timesteps.sort()

	return coordinate_timesteps,vx_timesteps,vy_timesteps,vz_timesteps


def write_surface_cage_file(x,y,filename):
	np_data_type="float64"

	x_out=np.array(x,dtype=np_data_type)
	y_out=np.array(y,dtype=np_data_type)
	outfile=open(filename,"wb")
	for idx in range(min(len(x_out),len(y_out))):
		outfile.write(x_out[idx])
		outfile.write(y_out[idx])
	outfile.close()


def read_surface_cage_file(filename):
	data_type="d"
	data_size=8

	filesize=os.path.getsize(filename)
	num_g=filesize//(data_size*2)
	x=np.zeros((num_g))
	y=np.zeros((num_g))
	infile=open(filename,"rb")
	for idx in range(num_g):
		x[idx]=struct.unpack(data_type, infile.read(data_size))[0]
		y[idx]=struct.unpack(data_type, infile.read(data_size))[0]
	infile.close()
	return x,y


def save_connectivity_map(y,x,connectivity_map,connectivity_filename):
	"""
	Store connectivity map
	"""
	len_x=len(connectivity_map[:,0])
	len_y=len(connectivity_map[0,:])
	outfile=open(connectivity_filename,'wb')
	outfile.write(struct.pack('i',len_x))
	outfile.write(struct.pack('i',len_y))
	for idx_x in range(len_x):
		outfile.write(struct.pack('f',x[idx_x,0]))
	for idx_y in range(len_y):
		outfile.write(struct.pack('f',y[0,idx_y]))
	for idx_x in range(len_x):
		outfile.write(bytes(connectivity_map[idx_x,:]))
	outfile.close()


def load_connectivity_map(connectivity_filename):
	infile=open(connectivity_filename,"rb")
	len_x=struct.unpack('i',infile.read(4))[0]
	len_y=struct.unpack('i',infile.read(4))[0]
	x_grid=np.zeros((len_x))
	y_grid=np.zeros((len_y))
	connectivity_map=np.zeros((len_x,len_y),dtype=np.int32)
	for idx_x in range(len_x):
		x_grid[idx_x]=struct.unpack('f',infile.read(4))[0]
	for idx_y in range(len_y):
		y_grid[idx_y]=struct.unpack('f',infile.read(4))[0]
	y_grid,x_grid=np.meshgrid(y_grid,x_grid)
	for idx_x in range(len_x):
		connectivity_map[idx_x,:]=np.array(struct.unpack('i'*len_y,infile.read(4*len_y)))
	infile.close()
	return y_grid,x_grid,connectivity_map


