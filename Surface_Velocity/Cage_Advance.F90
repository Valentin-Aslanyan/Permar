!To compile, run the following line (without leading exclamation mark):
!gfortran Permar_Functions_Fortran.F90 Cage_Advance.F90 -o Cage_Advance

!Filenames, directories, lines in control files are capped at 200 characters by len=200 magic number
!For specialized outputs from Lare3d, move points in a plane of constant z (as defined when outputting)
!Use Convert_SurfaceV_Files.py beforehand, ensure ghost cells are trimmed from all files

program Cage_Advance
      use Permar_Functions_Fortran
      implicit none
      call initialize_variables

      !Define variables below, then recompile to run
      vx_active = .TRUE.
      vy_active = .FALSE.
      vz_active = .FALSE.    !Does nothing for now
      x_periodic = .TRUE.
      x_mirrored = .FALSE.   !Does not work yet (will probably crash)
      y_periodic = .TRUE.
      y_mirrored = .FALSE.   !Does not work yet (will probably crash)

      !Use normalized units; comment out start_time to use start of output file
      !start_time = 0.0_num
      !Must be strictly in increasing order, use continuation character & after 132 characters
      output_times = (/ 0.8406499516990934_num,1.6814491781989176_num,2.5218007579086446_num,3.3623561378794125_num, &
      & 4.202985795644654_num,5.043727963985752_num,5.884264607019615_num,6.724924063948593_num,7.565429760289294_num, &
      & 8.40587774195161_num,9.246699209663571_num,10.086990799758821_num,10.927711682215376_num,11.768288245114924_num, &
      & 12.60887271441676_num,13.449333893229891_num,14.290139326760476_num,15.130739807798575_num,15.97127156844316_num, &
      & 16.811719178527138_num /)


      data_dir = './Data/'


      !!!!!!!!!!!!!!!!!!!!!!!
      ! Main program begins !
      !!!!!!!!!!!!!!!!!!!!!!!
      filename_full=TRIM(data_dir)//TRIM(t_filename)
      INQUIRE(FILE=filename_full, SIZE=num_t)
      num_t=num_t/num_sz
      open (unit=801,file=filename_full,form='unformatted',access='direct',recl=num_sz)

      if (num_t .lt. 2) then
        print *, 'Error, not enough timesteps in inputs'
        code_run=.false.
      else if (start_time .eq. -1.0E20_num) then !default start time at beginning of file
        code_run=.true.

        idx_t=1
        read (801,rec=idx_t) time1
        read (801,rec=idx_t+1) time2
        start_time=time1
      else !start_time has been set
        idx_t=1
        read (801,rec=idx_t) time1
        read (801,rec=idx_t+1) time2
        if (start_time .lt. time1) then
          print *, 'Error, start time earlier than input'
          code_run=.false.
        else
          do while ((start_time .ge. time2) .and. (idx_t .lt. num_t))
            idx_t=idx_t+1
            time1=time2
            read (801,rec=idx_t+1) time2
          end do
          if (idx_t .ge. num_t) then
            print *, 'Error, start time later than input'
            code_run=.false.
          else
            code_run=.true.
          end if
        end if
      end if

      if (code_run) then
        time_curr=start_time
        idx_ot=1
        do while (idx_ot .lt. SIZE(output_times) .and. output_times(idx_ot) .lt. time_curr)
          idx_ot=idx_ot+1
        end do
        if (output_times(idx_ot) .lt. time_curr) then
          idx_ot=idx_ot+1
        end if

        !Load cage, coords and velocity
        WRITE(cage_filename, '("SurfaceCage0.bin")')
        filename_full=TRIM(data_dir)//TRIM(cage_filename)
        call read_cage_file(filename_full,num_cage,cage_x,cage_y)
        filename_full=TRIM(data_dir)//TRIM(x_filename)
        call read_surface_coords_file(filename_full,num_x,coord_x)
        delta_x=coord_x(2)-coord_x(1)
        coord_x_regular=is_array_regular(coord_x,num_x)
        idx_x=1
        filename_full=TRIM(data_dir)//TRIM(y_filename)
        call read_surface_coords_file(filename_full,num_y,coord_y)
        delta_y=coord_y(2)-coord_y(1)
        coord_y_regular=is_array_regular(coord_y,num_y)
        idx_y=1
        if (vx_active) then
          allocate(vx1(num_x,num_y))
          allocate(vx2(num_x,num_y))
          filename_full=TRIM(data_dir)//TRIM(vx_filename)
          open (unit=803,file=filename_full,form='unformatted',access='direct',recl=num_sz*num_x*num_y)
          read (803,rec=idx_t) vx1
          read (803,rec=idx_t+1) vx2
        end if
        if (vy_active) then
          allocate(vy1(num_x,num_y))
          allocate(vy2(num_x,num_y))
          filename_full=TRIM(data_dir)//TRIM(vy_filename)
          open (unit=804,file=filename_full,form='unformatted',access='direct',recl=num_sz*num_x*num_y)
          read (804,rec=idx_t) vy1
          read (804,rec=idx_t+1) vy2
        end if
        if (vz_active) then
          allocate(vz1(num_x,num_y))
          allocate(vz2(num_x,num_y))
          filename_full=TRIM(data_dir)//TRIM(vz_filename)
          open (unit=805,file=filename_full,form='unformatted',access='direct',recl=num_sz*num_x*num_y)
          read (805,rec=idx_t) vz1
          read (805,rec=idx_t+1) vz2
        end if

        !Main loop
        do while ((idx_t .lt. num_t) .and. (idx_ot .le. SIZE(output_times)))
          if (output_times(idx_ot) .gt. time2) then
            dt=time2-time_curr
            next_timestep=.true.
          else
            dt=output_times(idx_ot)-time_curr
            next_timestep=.false.
          end if

          !Advance cage by one timestep
          do idx_g=1,num_cage
            rk_4_step=surface_move_rk_4(cage_x(idx_g),cage_y(idx_g),time_curr,dt,time1,time2)
            cage_x(idx_g)=rk_4_step(1)
            cage_y(idx_g)=rk_4_step(2)
          end do

          if (next_timestep) then
            idx_t=idx_t+1
            if (idx_t .lt. num_t) then
              time_curr=time2
              time1=time2
              read (801,rec=idx_t+1) time2
              if (vx_active) then
                vx1=vx2
                read (803,rec=idx_t+1) vx2
              end if
              if (vy_active) then
                vy1=vy2
                read (804,rec=idx_t+1) vy2
              end if
              if (vz_active) then
                vz1=vz2
                read (805,rec=idx_t+1) vz2
              end if
            end if
          else
            print *, "Outputting slice ",idx_ot," at time ",output_times(idx_ot)
            WRITE(cage_filename, '("SurfaceCage"I0".bin")') idx_ot
            filename_full=TRIM(data_dir)//TRIM(cage_filename)
            call write_cage_file(filename_full,num_cage,cage_x,cage_y)
            time_curr=output_times(idx_ot)
            idx_ot=idx_ot+1
          end if
        end do

        if (vx_active) then
          close(803)
        end if
        if (vy_active) then
          close(804)
        end if
        if (vz_active) then
          close(805)
        end if
      end if

      close(801)


end program Cage_Advance


