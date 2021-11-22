 program apply_incr_noahmp_snow

 use netcdf

 use NoahMPdisag_module, only : noahmp_type 

 implicit none

 include 'mpif.h'

! defined types

 type(noahmp_type)      :: noahmp_state

 integer :: res, len_land_vec
 character(len=8) :: date_str 
 character(len=2) :: hour_str

 ! restart arrays 
 type(noahmp_type) :: sfc_rst
 ! index to map between tile and vector space 
 integer, allocatable :: tile2vector(:,:)

 integer :: ierr, nprocs, myrank, lunit, ncid
 logical :: file_exists

 namelist /noahmp_snow/ date_str, hour_str, res 
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

    ! GET MAPPING INDEX (see subroutine comments re: source of land/sea mask)

    call get_fv3_mapping(myrank, date_str, hour_str, res, len_land_vec, tile2vector)
  
    ! SET-UP THE NOAH-MP STATE 

    allocate(noahmp_state%swe                (len_land_vec))
    allocate(noahmp_state%snow_depth         (len_land_vec))
    allocate(noahmp_state%active_snow_layers (len_land_vec))
    allocate(noahmp_state%swe_previous       (len_land_vec))
    allocate(noahmp_state%snow_soil_interface(len_land_vec,7))
    allocate(noahmp_state%temperature_snow   (len_land_vec,3))
    allocate(noahmp_state%snow_ice_layer     (len_land_vec,3))
    allocate(noahmp_state%snow_liq_layer     (len_land_vec,3))
    allocate(noahmp_state%temperature_soil   (len_land_vec))

    ! READ RESTART FILE 
    call   read_fv3_restart(myrank, date_str, hour_str, res, ncid, & 
                len_land_vec, tile2vector, noahmp_state )

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
! Get land sea mask from fv3 restart, and use to create 
! index for mapping from tiles (FV3 UFS restart) to vector
!  of land locations (offline Noah-MP restart)
! NOTE: slmsk in the restarts counts grid cells as land if 
!       they have a non-zero land fraction. Excludes grid 
!       cells that are surrounded by sea (islands). The slmsk 
!       in the oro_grid files (used by JEDI for screening out 
!       obs is different, and counts grid cells as land if they 
!       are more than 50% land (same exclusion of islands). If 
!       we want to change these definitations, may need to use 
!       land_frac field from the oro_grid files.
!--------------------------------------------------------------

 subroutine get_fv3_mapping(myrank, date_str, hour_str, res, & 
                len_land_vec, tile2vector)

 implicit none 

 include 'mpif.h'

 integer, intent(in) :: myrank, res
 character(len=8), intent(in) :: date_str 
 character(len=2), intent(in) :: hour_str 
 integer, allocatable, intent(out) :: tile2vector(:,:)
 integer :: len_land_vec

 character(len=100) :: restart_file
 character(len=1) :: rankch
 logical :: file_exists
 integer :: ierr,  ncid
 integer :: id_dim, id_var, fres
 integer :: slmsk(res,res) ! saved as double in the file, but i think this is OK
 integer :: i, j, nn

    ! OPEN FILE
    write(rankch, '(i1.1)') (myrank+1)
    restart_file = date_str//"."//hour_str//"0000.sfc_data.tile"//rankch//".nc"

    inquire(file=trim(restart_file), exist=file_exists)

    if (.not. file_exists) then
            print *, 'restart_file does not exist, ', &
                    trim(restart_file) , ' exiting'
            call mpi_abort(mpi_comm_world, 10) 
    endif

    write (6, *) 'calculate mapping from land mask in ', trim(restart_file)

    ierr=nf90_open(trim(restart_file),nf90_write,ncid)
    call netcdf_err(ierr, 'opening file: '//trim(restart_file) )

    ! READ MASK and GET MAPPING
    ierr=nf90_inq_varid(ncid, "slmsk", id_var)
    call netcdf_err(ierr, 'reading slmsk id' )
    ierr=nf90_get_var(ncid, id_var, slmsk)
    call netcdf_err(ierr, 'reading slmsk' )
 
    ! get number of land points  (note: slmsk is double)
    len_land_vec = 0
    do i = 1, res 
        do j = 1, res 
             if ( slmsk(i,j) == 1)  len_land_vec = len_land_vec+ 1  
        enddo 
    enddo
    
    write(6,*) 'Number of land points on rank ', myrank, ' :',  len_land_vec

    allocate(tile2vector(len_land_vec,2)) 
  
    nn=0
    do i = 1, res 
        do j = 1, res 
             if ( slmsk(i,j) == 1)   then 
                nn=nn+1
                tile2vector(nn,1) = i 
                tile2vector(nn,2) = j 
             endif
        enddo 
    enddo

end subroutine get_fv3_mapping


!--------------------------------------------------------------
! open fv3 restart, and read in required variables
! file is opened as read/write and remains open
!--------------------------------------------------------------
 subroutine read_fv3_restart(myrank, date_str, hour_str, res, ncid, & 
                len_land_vec,tile2vector, noahmp_state)

 implicit none 

 include 'mpif.h'

 integer, intent(in) :: myrank, res, len_land_vec
 character(len=8), intent(in) :: date_str 
 character(len=2), intent(in) :: hour_str 
 integer, intent(in) :: tile2vector(len_land_vec,2)

 integer, intent(out) :: ncid
 type(noahmp_type), intent(inout)  :: noahmp_state

 character(len=100) :: restart_file
 character(len=1) :: rankch
 logical :: file_exists
 integer :: ierr 
 integer :: id_dim, id_var, fres
 double precision :: dummy2D(res, res) 
 integer :: nn

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

    ! read swe (file name: sheleg) 
    ierr=nf90_inq_varid(ncid, "sheleg", id_var)
    call netcdf_err(ierr, 'reading sheleg id' )
    ierr=nf90_get_var(ncid, id_var, dummy2D)
    call netcdf_err(ierr, 'reading sheleg' )

    do nn=1,len_land_vec 
        noahmp_state%swe(nn) = dummy2D(tile2vector(nn,1), tile2vector(nn,2))
    enddo
 
    ! read swe (file name: sheleg, vert dim 1) 
    ! read snow_depth (file name: snwdph, vert dim 1)
    ! read active_snow_layers (file name: snowxy, vert dim: 1) 
    ! read swe_previous (file name: sneqvoxy, vert dim: 1) 
    ! read snow_soil_interface (file name: zsnsoxy, vert dim: 7) 
    ! read temperature_snow (file name: tsnoxy, vert dim: 3) 
    ! read snow_ice_layer (file name:  snicexy, vert dim: 3) 
    ! read snow_liq_layer (file name: snliqxy, vert dim: 3) 
    ! read temperature_soil (file name: stc, use layer 1 only, vert dim: 1) 


end subroutine read_fv3_restart
!--------------------------------------------------------------
! write updated fields tofv3_restarts  open on ncid
!--------------------------------------------------------------
 subroutine write_fv3_restart(ncid, noahmp_state, res) 

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




 end program apply_incr_noahmp_snow
