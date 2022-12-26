!To compile, run the following line (without leading exclamation mark):
!gfortran Permar_Functions_Fortran.F90 Connectivity_Map_Create.F90 -o Connectivity_Map_Create

!Filenames, directories, lines in control files are capped at 200 characters by len=200 magic number
!For specialized outputs from Lare3d, move points in a plane of constant z (as defined when outputting)
!Use Convert_SurfaceV_Files.py beforehand, ensure ghost cells are trimmed from all files

program Connectivity_Map_Generate
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

      !Use normalized units; connectivity map will be made relative to first listed time
      !Must be strictly in increasing order, use continuation character & after 132 characters
      output_times = (/ 0.0_num,0.8406499516990934_num,1.6814491781989176_num,2.5218007579086446_num,3.3623561378794125_num, &
      & 4.202985795644654_num,5.043727963985752_num,5.884264607019615_num,6.724924063948593_num,7.565429760289294_num, &
      & 8.40587774195161_num,9.246699209663571_num,10.086990799758821_num,10.927711682215376_num,11.768288245114924_num, &
      & 12.60887271441676_num,13.449333893229891_num,14.290139326760476_num /)

      !Full path to each QSL run corresponding to output_times
      QSL_files = [character(len=100) :: './QSL/0000/qslLare.qsl','./QSL/0001/qslLare.qsl','./QSL/0002/qslLare.qsl', &
      & './QSL/0003/qslLare.qsl','./QSL/0004/qslLare.qsl','./QSL/0005/qslLare.qsl','./QSL/0006/qslLare.qsl', &
      & './QSL/0007/qslLare.qsl','./QSL/0008/qslLare.qsl','./QSL/0009/qslLare.qsl','./QSL/0010/qslLare.qsl', &
      & './QSL/0011/qslLare.qsl','./QSL/0012/qslLare.qsl','./QSL/0013/qslLare.qsl','./QSL/0014/qslLare.qsl', &
      & './QSL/0015/qslLare.qsl','./QSL/0016/qslLare.qsl', './QSL/0017/qslLare.qsl' ]

      !Full path to file to be written
      connectivity_files = [character(len=100) :: './Data/Connectivity_0_0.bin','./Data/Connectivity_0_1.bin', &
      & './Data/Connectivity_0_2.bin','./Data/Connectivity_0_3.bin','./Data/Connectivity_0_4.bin', &
      & './Data/Connectivity_0_6.bin','./Data/Connectivity_0_7.bin','./Data/Connectivity_0_8.bin', &
      & './Data/Connectivity_0_9.bin','./Data/Connectivity_0_10.bin','./Data/Connectivity_0_11.bin', &
      & './Data/Connectivity_0_12.bin','./Data/Connectivity_0_13.bin','./Data/Connectivity_0_14.bin', &
      & './Data/Connectivity_0_15.bin','./Data/Connectivity_0_16.bin','./Data/Connectivity_0_17.bin', &
      & './Data/Connectivity_0_2.bin' ]


      data_dir = './Data/'


      !!!!!!!!!!!!!!!!!!!!!!!
      ! Main program begins !
      !!!!!!!!!!!!!!!!!!!!!!!
      filename_full=TRIM(data_dir)//TRIM(t_filename)
      INQUIRE(FILE=filename_full, SIZE=num_t)
      num_t=num_t/num_sz
      open (unit=801,file=filename_full,form='unformatted',access='direct',recl=num_sz)

      if ((SIZE(output_times) .ne. SIZE(QSL_files)) .and. (SIZE(output_times) .ne. SIZE(connectivity_files))) then
        print *, 'Different sizes for output_times, QSL_files and connectivity_files'
        code_run=.false.
      else if (num_t .lt. 2) then
        print *, 'Error, not enough timesteps in inputs'
        code_run=.false.
      else
        idx_ot=1
        read (801,rec=idx_ot) time1
        if (output_times(1) .lt. time1) then
          print *, 'Error, initial time earlier than input'
          code_run=.false.
        else
          do while ((output_times(1) .lt. time1) .and. (idx_ot .lt. num_t))
            idx_ot=idx_ot+1
            read (801,rec=idx_ot) time1
          end do
          if (idx_ot .ge. num_t) then
            print *, 'Error, initial time later than input'
            code_run=.false.
          else
            code_run=.true.
          end if
        end if
      end if

      if (code_run) then
        !Read initial Q, output first step
        call parse_QSL_2dbinfile(QSL_files(1),z_actual,num_Qy,num_Qx,Qy_start,Qx_start,Q_start)
        delta_Qx=Qx_start(2,1)-Qx_start(1,1)
        delta_Qy=Qy_start(1,2)-Qy_start(1,1)
	allocate(connectivity_map(num_Qx,num_Qy))
        do idx_Qx=1,num_Qx
          do idx_Qy=1,num_Qy
            if (Q_start(idx_Qx,idx_Qy)>0.0) then
              connectivity_map(idx_Qx,idx_Qy)=3
            else
              connectivity_map(idx_Qx,idx_Qy)=0
            end if
          end do
        end do
        print *, "Outputting ", TRIM(connectivity_files(1))
        call save_connectivity_map(num_Qy,num_Qx,Qy_start,Qx_start,connectivity_map,connectivity_files(1))

        !Load coords and velocity
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
        end if
        if (vy_active) then
          allocate(vy1(num_x,num_y))
          allocate(vy2(num_x,num_y))
          filename_full=TRIM(data_dir)//TRIM(vy_filename)
          open (unit=804,file=filename_full,form='unformatted',access='direct',recl=num_sz*num_x*num_y)
        end if
        if (vz_active) then
          allocate(vz1(num_x,num_y))
          allocate(vz2(num_x,num_y))
          filename_full=TRIM(data_dir)//TRIM(vz_filename)
          open (unit=805,file=filename_full,form='unformatted',access='direct',recl=num_sz*num_x*num_y)
        end if

        !Main loop
        do idx_f=2,SIZE(output_times)
          time_curr=output_times(idx_f)
          read (801,rec=1) time2
          read (801,rec=2) time1
          idx_st=2
          do while ((time_curr .gt. time1) .and. (idx_st .lt. num_t))
            idx_st=idx_st+1
            time2=time1
            read (801,rec=idx_st) time1
          end do
          if (time_curr .gt. time1) then
            print *, "Time ",time_curr," is outside available range"
            print *, "Latest available time is ",time1
          else
            call parse_QSL_2dbinfile(QSL_files(idx_f),z_actual,num_Qy,num_Qx,Qy_end,Qx_end,Q_end)
            if (vx_active) then
              read (803,rec=idx_st) vx1
              read (803,rec=idx_st-1) vx2
            end if
            if (vy_active) then
              read (804,rec=idx_st) vy1
              read (804,rec=idx_st-1) vy2
            end if
            if (vz_active) then
              read (805,rec=idx_st) vz1
              read (805,rec=idx_st-1) vz2
            end if

            do idx_t=idx_st-1,idx_ot+1,-1
              if (output_times(1) .gt. time2) then
                dt=output_times(1)-time_curr
              else
                dt=time2-time_curr
              end if

              !Advance cage by one timestep
              do idx_Qx=1,num_Qx
                do idx_Qy=1,num_Qy
                  rk_4_step=surface_move_rk_4(Qx_end(idx_Qx,idx_Qy),Qy_end(idx_Qx,idx_Qy),time_curr,dt,time1,time2)
                  Qx_end(idx_Qx,idx_Qy)=rk_4_step(1)
                  Qy_end(idx_Qx,idx_Qy)=rk_4_step(2)
                end do
              end do

              time_curr=time2
              time1=time2
              read (801,rec=idx_t-1) time2
              if (vx_active) then
                vx1=vx2
                read (803,rec=idx_t-1) vx2
              end if
              if (vy_active) then
                vy1=vy2
                read (804,rec=idx_t-1) vy2
              end if
              if (vz_active) then
                vz1=vz2
                read (805,rec=idx_t-1) vz2
              end if
            end do

            !Final timestep
            if (output_times(1) .gt. time2) then
              dt=output_times(1)-time_curr
            else
              dt=time2-time_curr
            end if
            do idx_Qx=1,num_Qx
              do idx_Qy=1,num_Qy
                rk_4_step=surface_move_rk_4(Qx_end(idx_Qx,idx_Qy),Qy_end(idx_Qx,idx_Qy),time_curr,dt,time1,time2)
                Qx_end(idx_Qx,idx_Qy)=rk_4_step(1)
                Qy_end(idx_Qx,idx_Qy)=rk_4_step(2)
              end do
            end do

            !Classify connectivity
            do idx_Qx=1,num_Qx
              do idx_Qy=1,num_Qy
                connectivity_map(idx_Qx,idx_Qy)=0
                if (Q_start(idx_Qx,idx_Qy)>0.0) then
                  connectivity_map(idx_Qx,idx_Qy)=connectivity_map(idx_Qx,idx_Qy)+1
                end if
                if (Q_end(idx_Qx,idx_Qy)>0.0) then
                  connectivity_map(idx_Qx,idx_Qy)=connectivity_map(idx_Qx,idx_Qy)+2
                end if
                idx_Qx2=MAX(MIN(NINT((Qx_end(idx_Qx,idx_Qy)-Qx_start(1,1))/delta_Qx)+1,num_Qx),1)
                idx_Qy2=MAX(MIN(NINT((Qy_end(idx_Qx,idx_Qy)-Qy_start(1,1))/delta_Qy)+1,num_Qy),1)
                if (Q_start(idx_Qx2,idx_Qy2)>0.0 .and. Q_end(idx_Qx,idx_Qy)<0.0) then
                  connectivity_map(idx_Qx,idx_Qy)=connectivity_map(idx_Qx,idx_Qy)+4
                else if (Q_start(idx_Qx2,idx_Qy2)<0.0 .and. Q_end(idx_Qx,idx_Qy)>0.0) then
                  connectivity_map(idx_Qx,idx_Qy)=connectivity_map(idx_Qx,idx_Qy)+4
                end if
              end do
            end do
            print *, "Outputting ", TRIM(connectivity_files(idx_f))
            call save_connectivity_map(num_Qy,num_Qx,Qy_start,Qx_start,connectivity_map,connectivity_files(idx_f))
            deallocate(Qy_end,Qx_end,Q_end)
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


end program Connectivity_Map_Generate


