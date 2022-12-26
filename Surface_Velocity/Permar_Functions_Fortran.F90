!Filenames, directories, lines in control files are capped at 200 characters by len=200 magic number
!!

module Permar_Functions_Fortran
      implicit none

      LOGICAL :: vx_active,vy_active,vz_active,x_periodic,x_mirrored,y_periodic,y_mirrored
      INTEGER, PARAMETER :: num = KIND(1.D0)
      INTEGER, PARAMETER :: num_sz = 8
      CHARACTER (LEN = 200) :: data_dir
      INTEGER :: num_x, num_y,num_cage, num_Qx, num_Qy
      REAL(num) :: start_time
      REAL(num), DIMENSION(:), ALLOCATABLE :: output_times
      INTEGER, DIMENSION(:), ALLOCATABLE :: infile_indices

      LOGICAL :: code_run, next_timestep
      LOGICAL :: coord_x_regular,coord_y_regular
      INTEGER :: idx_x,idx_y,idx_Qx,idx_Qy,idx_Qx2,idx_Qy2
      REAL(num) :: delta_x,delta_y,rk_4_step(2),delta_Qx,delta_Qy
      INTEGER :: idx_t,idx_st,idx_ot,idx_g,idx_f
      INTEGER :: num_t
      REAL(num) :: time1,time2,time_curr,dt,z_actual
      CHARACTER (LEN = 200) :: filename_full
      CHARACTER (LEN = 200) :: t_filename
      CHARACTER (LEN = 200) :: x_filename
      CHARACTER (LEN = 200) :: y_filename
      CHARACTER (LEN = 200) :: vx_filename
      CHARACTER (LEN = 200) :: vy_filename
      CHARACTER (LEN = 200) :: vz_filename
      CHARACTER (LEN = 200) :: cage_filename
      CHARACTER (LEN = 200), DIMENSION(:), ALLOCATABLE :: QSL_files
      CHARACTER (LEN = 200), DIMENSION(:), ALLOCATABLE :: connectivity_files
      REAL(num), DIMENSION(:), ALLOCATABLE :: cage_x
      REAL(num), DIMENSION(:), ALLOCATABLE :: cage_y
      REAL(num), DIMENSION(:), ALLOCATABLE :: coord_x
      REAL(num), DIMENSION(:), ALLOCATABLE :: coord_y
      REAL(num), DIMENSION(:,:), ALLOCATABLE :: vx1
      REAL(num), DIMENSION(:,:), ALLOCATABLE :: vx2
      REAL(num), DIMENSION(:,:), ALLOCATABLE :: vy1
      REAL(num), DIMENSION(:,:), ALLOCATABLE :: vy2
      REAL(num), DIMENSION(:,:), ALLOCATABLE :: vz1
      REAL(num), DIMENSION(:,:), ALLOCATABLE :: vz2
      REAL(num), DIMENSION(:,:), ALLOCATABLE :: Qx_start
      REAL(num), DIMENSION(:,:), ALLOCATABLE :: Qy_start
      REAL(num), DIMENSION(:,:), ALLOCATABLE :: Q_start
      REAL(num), DIMENSION(:,:), ALLOCATABLE :: Qx_end
      REAL(num), DIMENSION(:,:), ALLOCATABLE :: Qy_end
      REAL(num), DIMENSION(:,:), ALLOCATABLE :: Q_end
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: connectivity_map

      public :: initialize_variables, read_surface_coords_file,read_cage_file,write_cage_file,parse_QSL_2dbinfile

      contains


      subroutine initialize_variables
        start_time = -1.0E20_num
        t_filename = 'Surft.bin'
        x_filename = 'Surfx.bin'
        y_filename = 'Surfy.bin'
        vx_filename = 'SurfVx.bin'
        vy_filename = 'SurfVy.bin'
        vz_filename = 'SurfVz.bin'
      end subroutine initialize_variables


      subroutine read_surface_coords_file(filename,num_c,coordinates)
        CHARACTER (LEN = 200) :: filename
        INTEGER :: num_c
        REAL(num), DIMENSION(:), ALLOCATABLE :: coordinates

        INTEGER :: idx

        INQUIRE(FILE=filename, SIZE=num_c)
        num_c=num_c/num_sz
        allocate(coordinates(num_c))

        open (unit=802,file=filename,form='unformatted',access='direct',recl=num_sz)
        do idx=1,num_c
          read (802,rec=idx) coordinates(idx)
        end do
        close(802)
      end subroutine read_surface_coords_file


      subroutine read_cage_file(filename,num_cage,cage_x,cage_y)
        CHARACTER (LEN = 200) :: filename
        INTEGER :: num_cage
        REAL(num), DIMENSION(:), ALLOCATABLE :: cage_x
        REAL(num), DIMENSION(:), ALLOCATABLE :: cage_y

        INTEGER :: idx

        INQUIRE(FILE=filename, SIZE=num_cage)
        num_cage=num_cage/2/num_sz
        allocate(cage_x(num_cage))
        allocate(cage_y(num_cage))
        open (unit=806,file=filename,form='unformatted',access='direct',recl=num_sz)
        do idx=1,num_cage
          read (806,rec=idx*2-1) cage_x(idx)
          read (806,rec=idx*2) cage_y(idx)
        end do
        close(806)
      end subroutine read_cage_file


      subroutine write_cage_file(filename,num_cage,cage_x,cage_y)
        CHARACTER (LEN = 200) :: filename
        INTEGER :: num_cage
        REAL(num), DIMENSION(:), ALLOCATABLE :: cage_x
        REAL(num), DIMENSION(:), ALLOCATABLE :: cage_y

        INTEGER :: idx

        open (unit=806,file=filename,form='unformatted',status='REPLACE',recl=num_sz)
        close(806)
        open (unit=806,file=filename,form='unformatted',access='direct',recl=num_sz)
        do idx=1,num_cage
          write (806,rec=idx*2-1) cage_x(idx)
          write (806,rec=idx*2) cage_y(idx)
        end do
        close(806)
      end subroutine write_cage_file


      function point_in_bounds(point_in,coordinates,num_c,coord_periodic,coord_mirrored) result(point_out)
        REAL(num) :: point_in
        REAL(num), DIMENSION(:), ALLOCATABLE :: coordinates
        INTEGER :: num_c
        LOGICAL :: coord_periodic,coord_mirrored
        REAL(num) :: point_out

        if (coord_periodic) then
          point_out=MODULO(point_in-coordinates(1),coordinates(num_c)-coordinates(1))+coordinates(1)
        else if (coord_mirrored) then
          point_out=point_in !TODO: Change this
        else
          point_out=MIN(MAX(point_in,coordinates(1)),coordinates(num_c))
        end if
      end function point_in_bounds

      !f00=f(x0,y0),	f10=f(x1,y0), etc
      function bilinear_interpolation(x,y,x0,x1,y0,y1,f00,f01,f10,f11) result(f_interp)
        REAL(num) :: x,y,x0,x1,y0,y1,f00,f01,f10,f11
        REAL(num) :: f_interp

        REAL(num) :: shift_x0,shift_x1,shift_y0,shift_y1

        shift_x0=x-x0
        shift_x1=x1-x
        shift_y0=y-y0
        shift_y1=y1-y

        f_interp=(f00*shift_x1*shift_y1+f10*shift_x0*shift_y1+f01*shift_x1*shift_y0+f11*shift_x0*shift_y0) &
        &/(x1-x0)/(y1-y0)
      end function bilinear_interpolation


      function is_array_regular(coordinates,num_c) result(is_regular)
        REAL(num), DIMENSION(:), ALLOCATABLE :: coordinates
        INTEGER :: num_c
        LOGICAL :: is_regular

        INTEGER :: idx
        REAL(num) :: delta
        REAL(num) :: accuracy

        accuracy=1E-10_num
        delta=coordinates(2)-coordinates(1)
        is_regular=.true.
        do idx=2,num_c-1
          if (abs(coordinates(idx+1)-coordinates(idx)-delta) .gt. accuracy) then
            is_regular=.false.
          end if
        end do
      end function is_array_regular


      function find_index_regular(target,start,num_c,delta) result(idx)
        REAL(num) :: target,start,delta
        INTEGER :: num_c
        INTEGER :: idx

        idx=MIN(FLOOR((target-start)/delta)+1,num_c-1)
      end function find_index_regular


      function find_index_irregular(target,coordinates,num_c,idx_prev) result(idx)
        REAL(num) :: target
        REAL(num), DIMENSION(:), ALLOCATABLE :: coordinates
        INTEGER :: num_c,idx_prev
        INTEGER :: idx

        INTEGER :: idx2

        if ((coordinates(idx_prev) .le. target) .and. (coordinates(idx_prev+1) .ge. target)) then
          idx=idx_prev
        else
          idx2=1
          do while ((idx2 .lt. num_c-1) .and. (coordinates(idx2) .ge. target))
            idx2=idx2+1
          end do
        end if
      end function find_index_irregular


      subroutine get_velocity(x,y,t,t1,t2,vx_loc_curr,vy_loc_curr)
        REAL(num) :: x,y,t,t1,t2

        REAL(num) :: vx_loc1,vx_loc2,vx_loc_curr,vy_loc1,vy_loc2,vy_loc_curr

        if (coord_x_regular) then
          idx_x=find_index_regular(x,coord_x(1),num_x,delta_x)
        else
          idx_x=find_index_irregular(x,coord_x,num_x,idx_x)
        end if
        if (coord_y_regular) then
          idx_y=find_index_regular(y,coord_y(1),num_y,delta_y)
        else
          idx_y=find_index_irregular(y,coord_y,num_y,idx_y)
        end if
            
        if (vx_active) then
          vx_loc1=bilinear_interpolation(x,y,coord_x(idx_x),coord_x(idx_x+1),coord_y(idx_y),coord_y(idx_y+1),&
          &vx1(idx_x,idx_y),vx1(idx_x,idx_y+1),vx1(idx_x+1,idx_y),vx1(idx_x+1,idx_y+1))
          vx_loc2=bilinear_interpolation(x,y,coord_x(idx_x),coord_x(idx_x+1),coord_y(idx_y),coord_y(idx_y+1),&
          &vx2(idx_x,idx_y),vx2(idx_x,idx_y+1),vx2(idx_x+1,idx_y),vx2(idx_x+1,idx_y+1))
          vx_loc_curr=((t2-t)*vx_loc1+(t-t1)*vx_loc2)/(t2-t1)
        else
          vx_loc_curr=0.0_num
        end if
        if (vy_active) then
          vy_loc1=bilinear_interpolation(x,y,coord_x(idx_x),coord_x(idx_x+1),coord_y(idx_y),coord_y(idx_y+1),&
          &vy1(idx_x,idx_y),vy1(idx_x,idx_y+1),vy1(idx_x+1,idx_y),vy1(idx_x+1,idx_y+1))
          vy_loc2=bilinear_interpolation(x,y,coord_x(idx_x),coord_x(idx_x+1),coord_y(idx_y),coord_y(idx_y+1),&
          &vy2(idx_x,idx_y),vy2(idx_x,idx_y+1),vy2(idx_x+1,idx_y),vy2(idx_x+1,idx_y+1))
          vy_loc_curr=((t2-t)*vy_loc1+(t-t1)*vy_loc2)/(t2-t1)
        else
          vy_loc_curr=0.0_num
        end if
      end subroutine get_velocity


      function surface_move_rk_4(x_in,y_in,t,delta_t,t1,t2) result(arr_out)
        REAL(num) :: x_in,y_in,t,delta_t,t1,t2
        REAL(num) :: arr_out(2)

        REAL(num) :: x_curr,y_curr,vx_loc,vy_loc
        REAL(num) :: k_1_x,k_1_y,k_2_x,k_2_y,k_3_x,k_3_y,k_4_x,k_4_y,delta_t_half,delta_t_sixth

        delta_t_half=delta_t*0.5_num
        delta_t_sixth=delta_t/6.0_num

        x_curr=point_in_bounds(x_in,coord_x,num_x,x_periodic,x_mirrored)
        y_curr=point_in_bounds(y_in,coord_y,num_y,y_periodic,y_mirrored)
        call get_velocity(x_curr,y_curr,t,t1,t2,vx_loc,vy_loc)
        k_1_x=vx_loc
        k_1_y=vy_loc

        x_curr=point_in_bounds(x_curr+k_1_x*delta_t_half,coord_x,num_x,x_periodic,x_mirrored)
        y_curr=point_in_bounds(y_curr+k_1_y*delta_t_half,coord_y,num_y,y_periodic,y_mirrored)
        call get_velocity(x_curr,y_curr,t+delta_t_half,t1,t2,vx_loc,vy_loc)
        k_2_x=vx_loc
        k_2_y=vy_loc

        x_curr=point_in_bounds(x_curr+k_2_x*delta_t_half,coord_x,num_x,x_periodic,x_mirrored)
        y_curr=point_in_bounds(y_curr+k_2_y*delta_t_half,coord_y,num_y,y_periodic,y_mirrored)
        call get_velocity(x_curr,y_curr,t+delta_t_half,t1,t2,vx_loc,vy_loc)
        k_3_x=vx_loc
        k_3_y=vy_loc

        x_curr=point_in_bounds(x_curr+k_3_x*delta_t,coord_x,num_x,x_periodic,x_mirrored)
        y_curr=point_in_bounds(y_curr+k_3_y*delta_t,coord_y,num_y,y_periodic,y_mirrored)
        call get_velocity(x_curr,y_curr,t+delta_t,t1,t2,vx_loc,vy_loc)
        k_4_x=vx_loc
        k_4_y=vy_loc

        arr_out(1)=point_in_bounds(x_in+(k_1_x+2.0_num*k_2_x+2.0_num*k_3_x+k_4_x)*delta_t_sixth,coord_x,&
        &num_x,x_periodic,x_mirrored)
        arr_out(2)=point_in_bounds(y_in+(k_1_y+2.0_num*k_2_y+2.0_num*k_3_y+k_4_y)*delta_t_sixth,coord_y,&
        &num_y,y_periodic,y_mirrored)
      end function surface_move_rk_4


      subroutine parse_QSL_2dbinfile(filename,z_actual,num_y,num_x,y,x,Q)
        character (len = 200) :: filename
        REAL(num) :: z_actual
        INTEGER :: num_y
        INTEGER :: num_x
        REAL(num), DIMENSION(:,:), ALLOCATABLE :: y
        REAL(num), DIMENSION(:,:), ALLOCATABLE :: x
        REAL(num), DIMENSION(:,:), ALLOCATABLE :: Q

        INTEGER :: num_points,filesize,idx_y,idx_x,idx
        REAL(8) :: x_1,x_temp,y_temp,z_temp
        REAL(4) :: Q_temp

        inquire(file=filename,size=filesize)
        num_points=filesize/28

        open (unit=807,file=filename,form='unformatted',access='direct',recl=28)
        read (807,rec=1) x_1,y_temp,z_temp,Q_temp
        num_y=1
        idx_y=2
        do while (idx_y<=num_points)
          read (807,rec=idx_y) x_temp,y_temp,z_temp,Q_temp
          if (x_temp==x_1) then
            num_y=num_y+1
            idx_y=idx_y+1
          else
            idx_y=num_points+1
          end if
        end do
        close(807)
	
	num_x=num_points/num_y
	allocate(y(num_x,num_y))
	allocate(x(num_x,num_y))
	allocate(Q(num_x,num_y))

        idx=1	
        open (unit=807,file=filename,form='unformatted',access='direct',recl=28)
        do idx_x=1,num_x
          do idx_y=1,num_y
            read (807,rec=idx) x_temp,y_temp,z_temp,Q_temp
            y(idx_x,idx_y)=real(y_temp,num_sz)
            x(idx_x,idx_y)=real(x_temp,num_sz)
            Q(idx_x,idx_y)=real(Q_temp,num_sz)
            idx=idx+1
          end do
        end do
        close(807)
        z_actual=real(z_temp,num_sz)
      end subroutine parse_QSL_2dbinfile


      subroutine save_connectivity_map(num_y,num_x,y,x,connectivity_map,connectivity_filename)
        INTEGER :: num_y
        INTEGER :: num_x
        REAL(num), DIMENSION(:,:), ALLOCATABLE :: y
        REAL(num), DIMENSION(:,:), ALLOCATABLE :: x
        INTEGER, DIMENSION(:,:), ALLOCATABLE :: connectivity_map
        character (len = 200) :: connectivity_filename

        INTEGER :: idx_y,idx_x,idx

        open (unit=808,file=connectivity_filename,form='unformatted',status='REPLACE',recl=4)
        close(808)
        open (unit=808,file=connectivity_filename,form='unformatted',access='direct',recl=4)
        write(808,rec=1) num_x
        write(808,rec=2) num_y
        idx=3
        do idx_x=1,num_x
          write(808,rec=idx) REAL(x(idx_x,1),4)
          idx=idx+1
        end do
        do idx_y=1,num_y
          write(808,rec=idx) REAL(y(1,idx_y),4)
          idx=idx+1
        end do
        do idx_x=1,num_x
          do idx_y=1,num_y
            write(808,rec=idx) connectivity_map(idx_x,idx_y)
            idx=idx+1
          end do
        end do
        close(808)
      end subroutine save_connectivity_map


end module Permar_Functions_Fortran

