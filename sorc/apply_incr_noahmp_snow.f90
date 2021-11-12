 program apply_incr_noahmp_snow

 use netcdf

 implicit none

 include 'mpif.h'

! defined types
 type noahmp_type
    double precision, allocatable :: swe                (:)
    double precision, allocatable :: snow_depth         (:)
    double precision, allocatable :: active_snow_layers (:)
    double precision, allocatable :: swe_previous       (:)
    double precision, allocatable :: snow_soil_interface(:,:)
    double precision, allocatable :: temperature_snow   (:,:)
    double precision, allocatable :: snow_ice_layer     (:,:)
    double precision, allocatable :: snow_liq_layer     (:,:)
    double precision, allocatable :: temperature_soil   (:)
    double precision, allocatable :: sealand_mask       (:) ! should really be an integer, but is double in the file
 end type noahmp_type

 type observation_type
    double precision, allocatable :: snow_depth (:)
 end type observation_type

 type(noahmp_type)      :: noahmp_state
 type(observation_type) :: obs

 integer :: res, len_rst
 character(len=8) :: date_str 
 character(len=2) :: hour_str
 character(len=1) :: rankch

 ! restart arrays 
 type(noahmp_type) :: sfc_rst

 integer :: ierr, nprocs, myrank, lunit, ncid
 logical :: file_exists

 namelist /noahmp_snow/ date_str,hour_str, res 
!
    call mpi_init(ierr)
    call mpi_comm_size(mpi_comm_world, nprocs, ierr)
    call mpi_comm_rank(mpi_comm_world, myrank, ierr)

    print*
    print*,"starting apply_incr_noahmp_snow program on rank ", myrank

    ! READ NAMELIST 

    inquire (file='apply_incr_nml', exist=file_exists) 

    if (.not. file_exists) then
        write (6, *) 'ERROR: apply_incr_nml does not exist'
        call mpi_abort(mpi_comm_world, 10)
    end if

    open (action='read', file='apply_incr_nml', iostat=ierr, newunit=lunit)
    read (nml=noahmp_snow, iostat=ierr, unit=lunit)

    ! SET-UP THE NOAH-MP STATE 

    len_rst = res*res 

    allocate(noahmp_state%swe                (len_rst))
    allocate(noahmp_state%snow_depth         (len_rst))
    allocate(noahmp_state%active_snow_layers (len_rst))
    allocate(noahmp_state%swe_previous       (len_rst))
    allocate(noahmp_state%snow_soil_interface(len_rst,7))
    allocate(noahmp_state%temperature_snow   (len_rst,3))
    allocate(noahmp_state%snow_ice_layer     (len_rst,3))
    allocate(noahmp_state%snow_liq_layer     (len_rst,3))
    allocate(noahmp_state%temperature_soil   (len_rst))
    allocate(noahmp_state%sealand_mask       (len_rst))

    ! READ RESTART FILE 
    call   read_fv3_restart(myrank, date_str, hour_str, res, ncid, & 
                noahmp_state )

    call   write_fv3_restart(ncid, noahmp_state, res)


    ! CLOSE RESTART FILE 
    print*
    print*,"closing restart, apply_incr_noahmp_snow program on rank ", myrank
    ierr = nf90_close(ncid)

    call mpi_finalize(ierr)

 contains 

!--------------------------------------------------------------
! if at netcdf call returns an error, print out a message
! and stop processing.
!--------------------------------------------------------------
 subroutine netcdf_err( err, string )


        use netcdf

        implicit none

        include 'mpif.h'

        integer, intent(in) :: err
        character(len=*), intent(in) :: string
        character(len=80) :: errmsg

        if( err == nf90_noerr )return
        errmsg = nf90_strerror(err)
        print*,''
        print*,'fatal error: ', trim(string), ': ', trim(errmsg)
        print*,'stop.'
        call mpi_abort(mpi_comm_world, 999)

        return
 end subroutine netcdf_err


!--------------------------------------------------------------
! open fv3 restart, and read in required variables
! file is opened as read/write and remains open
!--------------------------------------------------------------
 subroutine read_fv3_restart(myrank, date_str, hour_str, res, ncid, & 
                noahmp_state)

 use netcdf 

 implicit none 

 include 'mpif.h'

 integer, intent(in) :: myrank, res
 character(len=8), intent(in) :: date_str 
 character(len=2), intent(in) :: hour_str 

 integer, intent(inout) :: ncid
 type(noahmp_type), intent(out)  :: noahmp_state

 character(len=100) :: restart_file
 character(len=1) :: rankch
 logical :: file_exists
 integer :: ierr
 integer :: id_dim, id_var, fres
 double precision :: dummy2D(res, res) 

    ! OPEN FILE
    write(rankch, '(i1.1)') (myrank+1)
    restart_file = date_str//"."//hour_str//"0000.sfc_data.tile"//rankch//".nc"

    inquire(file=trim(restart_file), exist=file_exists)

    if (.not. file_exists) then
            print *, 'restart_file does not exist, ', &
                    trim(restart_file) , ' exiting'
            call mpi_abort(mpi_comm_world, 10) 
    endif

    write (6, *) 'opening ', trim(restart_file)

    ierr=nf90_open(trim(restart_file),nf90_write,ncid)
    call netcdf_err(ierr, 'opening file: '//trim(restart_file) )

    ! CHECK DIMENSIONS
    ierr=nf90_inq_dimid(ncid, 'xaxis_1', id_dim)
    call netcdf_err(ierr, 'reading xaxis_1' )
    ierr=nf90_inquire_dimension(ncid,id_dim,len=fres)
    call netcdf_err(ierr, 'reading xaxis_1' )

    if ( fres /= res) then
       print*,'fatal error: dimensions wrong.'
       call mpi_abort(mpi_comm_world, ierr)
    endif

    ! READ VARS
    ierr=nf90_inq_varid(ncid, "slmsk", id_var)
    call netcdf_err(ierr, 'reading slmsk id' )
    ierr=nf90_get_var(ncid, id_var, dummy2D)
    call netcdf_err(ierr, 'reading slmsk' )
 
    noahmp_state%sealand_mask = reshape(dummy2D, (/res*res/))


!slmsk, sheleg, snwdph, sncovr, stc


end subroutine read_fv3_restart


!--------------------------------------------------------------
! write updated fields tofv3_restarts  open on ncid
!--------------------------------------------------------------
 subroutine write_fv3_restart(ncid, noahmp_state, res) 

 use netcdf 

 implicit none 

 integer, intent(in) :: ncid, res
 type(noahmp_type), intent(in) :: noahmp_state

 !double precision :: dummy2D(res, res) 

 integer :: ierr, id_var 

     ierr=nf90_inq_varid(ncid, "slmsk", id_var)
     call netcdf_err(ierr, 'reading slmsk id' )
     
     !dummy2d = reshape(noahmp_state%slmsk, (/res,res/))
     !ierr = nf90_put_var( ncid, id_var, slmsk)
     !call netcdf_err(ierr, 'writing slmsk' )

 end subroutine write_fv3_restart

!--------------------------------------------------------------
! add increment to all snow layers, maintaining previous distribution
! between layers, and maintaining previous density. 
!--------------------------------------------------------------
  subroutine UpdateAllLayers(vector_length, noahmp, obs, increment)
  
  type(noahmp_type)      :: noahmp
  type(observation_type) :: obs
  integer                :: vector_length
  double precision       :: increment(vector_length), to_remove
  double precision       :: layer_density, swe_increment, liq_ratio, remove_ratio
  integer                :: iloc, ilayer, iinter, active_layers, vector_loc, pathway, removed
  double precision       :: soil_interfaces(7) = (/0.0,0.0,0.0,0.1,0.4,1.0,2.0/)
  double precision       :: partition_ratio, layer_depths(3)
  
  associate( &
      obs_snow_depth => obs%snow_depth            ,&
                 swe => noahmp%swe                ,&
          snow_depth => noahmp%snow_depth         ,&
  active_snow_layers => noahmp%active_snow_layers ,&
        swe_previous => noahmp%swe_previous       ,&
 snow_soil_interface => noahmp%snow_soil_interface,&
    temperature_snow => noahmp%temperature_snow   ,&
      snow_ice_layer => noahmp%snow_ice_layer     ,&
      snow_liq_layer => noahmp%snow_liq_layer     ,&
    temperature_soil => noahmp%temperature_soil )

  
!   increment = obs_snow_depth - snow_depth  ! snow to add or remove [mm]
  
  do iloc = 1, vector_length
    
    pathway = 0
    
    if(obs_snow_depth(iloc) == 0.0) then

      swe                (iloc)   = 0.0
      snow_depth         (iloc)   = 0.0
      active_snow_layers (iloc)   = 0.0
      swe_previous       (iloc)   = 0.0
      snow_soil_interface(iloc,:) = (/0.0,0.0,0.0,-0.1,-0.4,-1.0,-2.0/)
      temperature_snow   (iloc,:) = 0.0
      snow_ice_layer     (iloc,:) = 0.0
      snow_liq_layer     (iloc,:) = 0.0

    else
    
      active_layers = nint(active_snow_layers(iloc))  ! number of active layers (0,-1,-2,-3)

      if(active_layers < 0) then  ! in multi-layer mode
      
        layer_depths(1) = snow_soil_interface(iloc,1)
        layer_depths(2) = snow_soil_interface(iloc,2)-snow_soil_interface(iloc,1)
        layer_depths(3) = snow_soil_interface(iloc,3)-snow_soil_interface(iloc,2)

!if(iloc ==  10962)then
!print*, 'depths',  iloc,obs_snow_depth(iloc)  ,snow_depth(iloc)
!print*, 'increment',  increment(iloc)  ,layer_depths
!print*, 'interfaces',    snow_soil_interface(iloc,:)
!print*, 'ice',    snow_ice_layer(iloc,:), swe(iloc)
!print*, 'liq',    snow_liq_layer(iloc,:)
!end if
        if(increment(iloc) > 0.0) then  ! add snow in multi-layer mode

          pathway = 1
    
          vector_loc = 4 + active_layers  ! location in vector of top layer
          
          layerloop: do ilayer = vector_loc, 3
          
            partition_ratio = -layer_depths(ilayer)/snow_depth(iloc)*1000.d0
            layer_density = (snow_ice_layer(iloc,ilayer)+snow_liq_layer(iloc,ilayer)) / &
                              (-snow_soil_interface(iloc,ilayer))
            swe_increment = partition_ratio * increment(iloc) * layer_density / 1000.d0
            liq_ratio = snow_liq_layer(iloc,ilayer) / &
                          ( snow_ice_layer(iloc,ilayer) + snow_liq_layer(iloc,ilayer) )
            snow_ice_layer(iloc,ilayer) = snow_ice_layer(iloc,ilayer) + &
                                              (1.0 - liq_ratio) * swe_increment
            snow_liq_layer(iloc,ilayer) = snow_liq_layer(iloc,ilayer) + &
                                              liq_ratio * swe_increment
            do iinter = ilayer, 3  ! remove snow from each snow layer
              snow_soil_interface(iloc,iinter) = snow_soil_interface(iloc,iinter) - &
                                                   partition_ratio * increment(iloc)/1000.d0
            end do

          end do layerloop
            
        elseif(increment(iloc) < 0.0) then  ! remove snow in multi-layer mode

          pathway = 2
          
          vector_loc = 4 + active_layers  ! location in vector of top layer
          
          layerloop: do ilayer = vector_loc, 3
          
            partition_ratio = -layer_depths(ilayer)/snow_depth(iloc)*1000.d0
            layer_density = (snow_ice_layer(iloc,ilayer)+snow_liq_layer(iloc,ilayer)) / &
                              (-layer_depths(ilayer))
            swe_increment = partition_ratio * increment(iloc) * layer_density / 1000.d0
            liq_ratio = snow_liq_layer(iloc,ilayer) / &
                          ( snow_ice_layer(iloc,ilayer) + snow_liq_layer(iloc,ilayer) )
            snow_ice_layer(iloc,ilayer) = snow_ice_layer(iloc,ilayer) + &
                                              (1.0 - liq_ratio) * swe_increment
            snow_liq_layer(iloc,ilayer) = snow_liq_layer(iloc,ilayer) + &
                                              liq_ratio * swe_increment
            do iinter = ilayer, 3  ! remove snow from each snow layer
              snow_soil_interface(iloc,iinter) = snow_soil_interface(iloc,iinter) - &
                                                   partition_ratio * increment(iloc)/1000.d0
            end do

          end do layerloop
          
        end if  ! increment
        
        ! For multi-layer mode, recalculate interfaces and sum depth/swe

        do ilayer = 4, 7
          snow_soil_interface(iloc,ilayer) = snow_soil_interface(iloc,3) - soil_interfaces(ilayer)
        end do

        snow_depth(iloc) = -snow_soil_interface(iloc,3) * 1000.d0

        swe(iloc) = 0.0

        do ilayer = 1, 3
          swe(iloc) = swe(iloc) + snow_ice_layer(iloc,ilayer) + snow_liq_layer(iloc,ilayer)
        end do

        swe_previous(iloc) = swe(iloc)

        if(snow_depth(iloc) < 25.d0) then  ! go out of multi-layer mode
          active_snow_layers (iloc) = 0.d0
          snow_soil_interface(iloc,:) = (/0.0,0.0,0.0,-0.1,-0.4,-1.0,-2.0/)
          temperature_snow   (iloc,:) = 0.0
          snow_ice_layer     (iloc,:) = 0.0
          snow_liq_layer     (iloc,:) = 0.0
        end if

      elseif(active_layers == 0) then  ! snow starts in zero-layer mode

        if(increment(iloc) > 0.0) then  ! add snow in zero-layer mode
    
          if(snow_depth(iloc) == 0) then   ! no snow present, so assume density based on soil temperature
            pathway = 3
            layer_density = max(80.0,min(120.,67.92+51.25*exp((temperature_soil(iloc)-273.15)/2.59)))
          else   ! use existing density
            pathway = 4
            layer_density = swe(iloc) / snow_depth(iloc) * 1000.d0
          end if
          snow_depth(iloc) = snow_depth(iloc) + increment(iloc)
          swe(iloc) = swe(iloc) + increment(iloc) * layer_density / 1000.d0
          swe_previous(iloc) = swe(iloc)

          active_snow_layers(iloc)      = 0.0
          snow_ice_layer(iloc,:)        = 0.0
          snow_liq_layer(iloc,:)        = 0.0
          temperature_snow(iloc,:)      = 0.0
          snow_soil_interface(iloc,1:3) = 0.0

          if(snow_depth(iloc) > 25.0) then  ! snow depth is > 25mm so put in a layer
            pathway = 5
            active_snow_layers(iloc) = -1.0
            snow_ice_layer(iloc,3)   = swe(iloc)
            temperature_snow(iloc,3) = temperature_soil(iloc)
            do ilayer = 3, 7
              snow_soil_interface(iloc,ilayer) = snow_soil_interface(iloc,ilayer) - snow_depth(iloc)/1000.d0
            end do
          end if
          
        elseif(increment(iloc) < 0.0) then  ! remove snow in zero-layer mode

          pathway = 6
    
          if(snow_depth(iloc) <= 0.0) stop "inconsistency in snow_depth and increment"
          layer_density = swe(iloc) / snow_depth(iloc) * 1000.d0
          snow_depth(iloc) = snow_depth(iloc) + increment(iloc)
          swe(iloc) = swe(iloc) + increment(iloc) * layer_density / 1000.d0
          swe_previous(iloc) = swe(iloc)

          active_snow_layers(iloc)      = 0.0
          snow_ice_layer(iloc,:)        = 0.0
          snow_liq_layer(iloc,:)        = 0.0
          temperature_snow(iloc,:)      = 0.0
          snow_soil_interface(iloc,1:3) = 0.0

        end if  ! increment

      end if  ! active_layers
        
    end if  ! obs_snow_depth == 0
    
    ! do some gross checks

    if(abs(snow_soil_interface(iloc,7) - snow_soil_interface(iloc,3) + 2.d0) > 0.0000001) then
      print*, "Depth of soil not 2m"
      print*, pathway
      print*, snow_soil_interface(iloc,7), snow_soil_interface(iloc,3)
!      stop
    end if

    if(active_snow_layers(iloc) < 0.0 .and. abs(snow_depth(iloc) + 1000.d0*snow_soil_interface(iloc,3)) > 0.0000001) then
      print*, "snow_depth and snow_soil_interface inconsistent"
      print*, pathway
      print*, active_snow_layers(iloc), snow_depth(iloc), snow_soil_interface(iloc,3)
!      stop
    end if

    if(abs(obs_snow_depth(iloc) - snow_depth(iloc)) > 0.0000001) then
      print*, "observed snow and updated model snow inconsistent"
      print*, pathway
      print*, obs_snow_depth(iloc), snow_depth(iloc)
!      stop
    end if

    if(snow_depth(iloc) < 0.0 .or. snow_soil_interface(iloc,3) > 0.0 ) then
      print*, "observed snow and updated model snow inconsistent"
      print*, pathway
      print*, snow_depth(iloc), snow_soil_interface(iloc,3)
!      stop
    end if

  end do
  
  end associate
   
  end subroutine UpdateAllLayers
 
 end program apply_incr_noahmp_snow
