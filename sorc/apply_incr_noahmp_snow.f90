 program apply_incr_noahmp_snow

 use netcdf

 implicit none

 include 'mpif.h'

 integer :: res 
 character(len=8) :: date_str 
 character(len=2) :: hour_str
 character(len=100) :: restart_file
 character(len=1) :: rankch

 integer :: ierr, nprocs, myrank, lunit
 logical :: file_exists
 integer :: ncid, id_dim, fres

    

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

    ! READ RESTART FILE 

    write(rankch, '(i1.1)') (myrank+1)

    restart_file = date_str//"."//hour_str//"0000.sfc_data.tile"//rankch//".nc"

    inquire(file=trim(restart_file), exist=file_exists)

    if (.not. file_exists) then
            print *, 'restart_file does not exist, ', &
                    trim(restart_file) , ' exiting'
            call mpi_abort(mpi_comm_world, 10) 
    endif

    write (6, *) 'opening ', trim(restart_file)

    ierr=nf90_open(trim(restart_file),nf90_nowrite,ncid)
!    call netcdf_err(ierr, 'opening file: '//trim(restart_file) )

    ierr=nf90_inq_dimid(ncid, 'xaxis_1', id_dim)
!    call netcdf_err(ierr, 'reading xaxis_1' )
    ierr=nf90_inquire_dimension(ncid,id_dim,len=fres)
!    call netcdf_err(ierr, 'reading xaxis_1' )

    if ( fres /= res) then
       print*,'fatal error: dimensions wrong.'
       call mpi_abort(mpi_comm_world, ierr)
    endif

    call mpi_finalize(ierr)

 end program apply_incr_noahmp_snow

 subroutine netcdf_err( err, string )

    !--------------------------------------------------------------
    ! if at netcdf call returns an error, print out a message
    ! and stop processing.
    !--------------------------------------------------------------

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


