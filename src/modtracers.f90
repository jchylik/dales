!> \file modtracers.f90
!! Definitions and functions for passive and reactive tracers

!>
!!  \author Ruud Janssen, TNO
!  This file is part of DALES.
!
! DALES is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 3 of the License, or
! (at your option) any later version.
!
! DALES is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.
!
!  Copyright 1993-2023 Delft University of Technology, Wageningen
!  University, Utrecht University, KNMI, TNO
!

module modtracers

  use modglobal,      only: nsv, i1, ih, j1, jh, k1, kmax, cexpnr
  use modtracer_type, only: T_tracer
  use modprecision,   only: field_r
  use modfields,      only: svm, sv0, svp, sv0av, svprof
  use go,             only: goSplitString_s
  use modstat_nc
  use utils

  implicit none

  private

  save

  public :: inittracers
  public :: add_tracer
  public :: allocate_tracers
  public :: exittracers
  public :: tracer_profs_from_netcdf

  integer :: iname
  integer :: nsv_user !< Number of user-provided tracers

  type(T_tracer), allocatable, public, protected :: tracer_prop(:) !< List of tracers
  logical,                     protected         :: ltracers = .false.
  character(6),                protected         :: &
    tracernames(200) = (/ ('', iname=1, 200)/)            !< For compatibility

  logical :: file_exists

contains

  !> Initialize tracer definition.
  !!
  !! Read the namelist NAMTRACERS from the namoptions file, distribute
  !! the parameters to all processes and allocate the tracer (SV) arrays.
  subroutine inittracers
    ! read namelist
    ! init tracer type

    use modglobal,        only : ifnamopt, fname_options, checknamelisterror
    use modmpi,           only : myid, comm3d, d_mpi_bcast

    implicit none

    ! Auxiliary variables
    integer :: ierr

    ! Namelist definition
    namelist /NAMTRACERS/ &
        ltracers, tracernames

    ! Read namelist
    if (myid == 0) then
        open(ifnamopt, file=fname_options, status='old', iostat=ierr)
        read(ifnamopt, NAMTRACERS, iostat=ierr)
        call checknamelisterror(ierr, ifnamopt, 'NAMTRACERS')
        write(6, NAMTRACERS)
        close(ifnamopt)
    end if


    ! Broadcast namelist values to all MPI tasks
    call d_mpi_bcast(ltracers,             1, 0, comm3d, ierr)
    call d_mpi_bcast(tracernames(1:200), 200, 0, comm3d, ierr)

    if (ltracers) then
      call tracer_props_from_ascii('scalar.inp.'//cexpnr, 'tracerdata.inp')
    else
      call tracer_props_from_netcdf('tracers.'//cexpnr//'.nc')
    end if ! ltracers


  end subroutine inittracers

  !> \brief Define a new tracer
  !!
  !! \param name Short name of the tracer.
  !! \param long_name Full name of the tracer.
  !! \param unit Unit.
  !! \param molar_mass Molar mass (g mol^-1)
  !! \param lemis Tracer is emitted.
  !! \param lreact Tracer is reactive.
  !! \param ldep Tracer is deposited.
  !! \param lags Tracer is photosynthesized.
  !! \param lmicro Tracer is involved in cloud microphysics.
  !! \param laero Tracer is involved in aerosol microphyiscs.
  !! \note All tracers should be added before readinitfiles is called!
  subroutine add_tracer(name, long_name, unit, molar_mass, lemis, lreact, &
                        ldep, lags, lmicro, isv)
    character(*),  intent(in)            :: name
    character(*),  intent(in),  optional :: long_name
    character(*),  intent(in),  optional :: unit
    real(field_r), intent(in),  optional :: molar_mass
    logical,       intent(in),  optional :: lemis
    logical,       intent(in),  optional :: lreact
    logical,       intent(in),  optional :: ldep
    logical,       intent(in),  optional :: lags
    logical,       intent(in),  optional :: lmicro
    integer,       intent(out), optional :: isv

    integer                     :: s
    type(T_tracer), allocatable :: tmp(:)

    ! Check if the tracer already exists. If so, don't add a new one.
    if (nsv > 0) then
      do s = 1, nsv
        if (trim(to_lower(name)) == trim(to_lower(tracer_prop(s) % tracname))) then
          write(*,*) "Tracer", name, "already exists!"
          if(present(isv)) isv = s
          return
        end if
      end do
    end if

    nsv = nsv + 1

    if (.not. allocated(tracer_prop)) then
      allocate(tracer_prop(nsv))
    else
      allocate(tmp(nsv))
      tmp(1:nsv-1) = tracer_prop(1:nsv-1)
      call move_alloc(tmp, tracer_prop)
    end if

    tracer_prop(nsv) % tracname = name
    tracer_prop(nsv) % trac_idx = nsv
    if (present(long_name)) tracer_prop(nsv) % traclong = trim(long_name)
    if (present(unit)) tracer_prop(nsv) % unit = unit
    if (present(molar_mass)) tracer_prop(nsv) % molar_mass = molar_mass
    if (present(lemis)) tracer_prop(nsv) % lemis = lemis
    if (present(lreact)) tracer_prop(nsv) % lreact = lreact
    if (present(ldep)) tracer_prop(nsv) % ldep = ldep
    if (present(lags)) tracer_prop(nsv) % lags = lags
    if (present(lmicro)) tracer_prop(nsv) % lmicro = lmicro

    if (present(isv)) isv = nsv

  end subroutine add_tracer

  !> Allocates all tracer fields
  subroutine allocate_tracers

    integer :: s

    do s = 1, nsv
      call tracer_prop(s) % print_properties()
    end do

    allocate(svm(2-ih:i1+ih,2-jh:j1+jh,k1,nsv), &
             sv0(2-ih:i1+ih,2-jh:j1+jh,k1,nsv), &
             svp(2-ih:i1+ih,2-jh:j1+jh,k1,nsv), &
             sv0av(k1,nsv), svprof(k1,nsv))

    svm(:,:,:,:) = 0
    sv0(:,:,:,:) = 0
    svp(:,:,:,:) = 0
    sv0av(:,:) = 0
    svprof(:,:) = 0

    !$acc enter data copyin(svm(2-ih:i1+ih,2-jh:j1+jh,1:k1,1:nsv), &
    !$acc&                  sv0(2-ih:i1+ih,2-jh:j1+jh,1:k1,1:nsv), &
    !$acc&                  svp(2-ih:i1+ih,2-jh:j1+jh,1:k1,1:nsv), &
    !$acc&                  sv0av(1:k1,1:nsv), svprof(1:k1,1:nsv))

  end subroutine allocate_tracers

  !> Deallocates all tracers fields
  subroutine exittracers

    !$acc exit data delete(svm(2-ih:i1+ih,2-jh:j1+jh,1:k1,1:nsv), &
    !$acc&                 sv0(2-ih:i1+ih,2-jh:j1+jh,1:k1,1:nsv), &
    !$acc&                 svp(2-ih:i1+ih,2-jh:j1+jh,1:k1,1:nsv), &
    !$acc&                 sv0av(1:k1,1:nsv), svprof(1:k1,1:nsv))

    if (nsv > 0) then
      deallocate(tracer_prop)
      deallocate(svm, sv0, svp, sv0av, svprof)
    end if

  end subroutine exittracers

  !> \brief Setup tracers from ASCII input files
  !!
  !! \param file_profiles Name of file containing profiles of tracers. From 
  !! this file, the number of tracers is determined.
  !! \param file_properties Name of file containing tracer attributes.
  subroutine tracer_props_from_ascii(file_profiles, file_properties)
    character(len=*), intent(in) :: file_profiles
    character(len=*), intent(in) :: file_properties

    integer, parameter :: max_tracs  =  100 !<  Max. number of tracers that can be defined

    character(len=512) :: line
    character(len=7)   :: headers(max_tracs)
    integer            :: ierr
    integer            :: nheader
    integer            :: isv, n

    ! Buffers for tracer properties
    character(len=10) :: tracname_short(max_tracs) = "NA" ! Short tracer name
    character(len=32) :: tracname_long(max_tracs)  = "NA" ! Long tracer name
    character(len=10) :: tracer_unit(max_tracs) = "NA" ! Unit of tracer
    real(field_r)     :: molar_mass(max_tracs) = -1.0 ! Molar mass of tracer (g mol-1)
    logical           :: tracer_is_emitted(max_tracs) = .false. ! Tracer is emitted (T/F)
    logical           :: tracer_is_reactive(max_tracs) = .false. ! Tracer is reactive (T/F)
    logical           :: tracer_is_deposited(max_tracs) = .false. ! Tracer is deposited (T/F)
    logical           :: tracer_is_photosynth(max_tracs) = .false. ! Tracer is photosynthesized (T/F)
    logical           :: tracer_is_microphys(max_tracs) = .false. ! Tracer is involved in cloud microphysics (T/F)

    open(1, file=file_profiles, status='old', iostat=ierr)

    if (ierr /= 0) then
      if (myid == 0) write(6,*) "Error opening ", trim(file_profiles) ! DALESERROR
      error stop
    end if

    read(1, '(a512)') line
    read(1, '(a512)') line

    ! Determine the number of tracers from the header
    call goSplitString_s(line, nsv_user, headers, ierr, sep=' ')

    close(1)

    ! Read table of tracer properties, and put into the respective buffers
    open(1, file=file_properties, status='old', iostat=ierr)

    if (ierr /= 0) then
      if (myid == 0) write(6,*) "Error opening ", trim(file_properties) ! DALESERROR
      error stop
    end if

    ierr = 0
    isv = 0
    do while (ierr == 0)
      read(1, '(A)', iostat=ierr) line

      if (ierr == 0) then ! So no end of file is encountered
        if (line(1:1)=='#') then
          if (myid == 0)   print *, trim(line)
        else
          if (myid == 0)   print *, trim(line)
          isv = isv + 1
          read(line, *, iostat=ierr) tracname_short(isv), &
                                     tracname_long(isv), &
                                     tracer_unit(isv), &
                                     molar_mass(isv), &
                                     tracer_is_emitted(isv), &
                                     tracer_is_reactive(isv), &
                                     tracer_is_deposited(isv), &
                                     tracer_is_photosynth(isv), &
                                     tracer_is_microphys(isv)

          if ( (molar_mass(isv) < 0) .and. (molar_mass(isv) > -1.1) ) then
            if (myid == 0) then
              stop "MODTRACERS: a molar mass value is not set in the tracer &
                &input file"
            end if
          end if
        end if
      end if
    end do

    close(1)

    ! For every tracer, find the properties
    do n = 2, nsv_user ! Skip the first column, containing the heights
      call add_tracer( &
        name=trim(headers(n)), &
        long_name=trim(findval(headers(n), tracname_short, tracname_long, &
                               defltvalue='dummy longname')), & ! Default is 'dummy '
        unit=trim(findval(headers(n), tracname_short, &
                    tracer_unit, defltvalue='dummy unit')), & ! Default is 'dummy unit'
        molar_mass=findval(headers(n), tracname_short, &
                     molar_mass, defltvalue=-999._field_r), & ! Default is -999.
        lemis=findval(headers(n), tracname_short, &
                tracer_is_emitted, defltvalue=.false.), & ! Default is False
        lreact=findval(headers(n), tracname_short, &
                 tracer_is_reactive, defltvalue=.false.), & ! Default is False
        ldep=findval(headers(n), tracname_short, &
               tracer_is_deposited, defltvalue=.false.), & ! Default is False
        lags=findval(headers(n), tracname_short, &
               tracer_is_photosynth, defltvalue=.false.), & ! Default is False
        lmicro=findval(headers(n), tracname_short, &
                 tracer_is_microphys, defltvalue=.false.) & ! Default is False
      )
    end do
    
  end subroutine tracer_props_from_ascii

  !> \brief Read tracer properties from tracers.XXX.nc
  !!
  !! \param filename Name of the input file.
  subroutine tracer_props_from_netcdf(filename)
    character(*),  intent(in)  :: filename

    integer              :: ncid, nvars
    integer              :: ivar
    integer, allocatable :: varids(:)

    ! Tracer attributes
    character(NF90_MAX_NAME) :: name
    character(32) :: long_name
    character(16) :: unit
    real(field_r) :: molar_mass
    logical       :: lemis, lreact, ldep, lags

    inquire(file=filename, exist=file_exists)

    if (.not. file_exists) then
      write(6,*) "Warning: ", filename, " not found."
      return
    end if

    call nchandle_error(nf90_open(filename, NF90_NOWRITE, ncid))
    call nchandle_error(nf90_inquire(ncid, nVariables=nvars))

    allocate(varids(nvars))

    call nchandle_error(nf90_inq_varids(ncid, nvars, varids))

    do ivar = 1, nvars
      call nchandle_error(nf90_inquire_variable(ncid, varids(ivar), name=name))

      ! Read attributes
      call read_nc_attribute(ncid, varids(ivar), "long_name", long_name, default="dummy")
      call read_nc_attribute(ncid, varids(ivar), "unit", unit, default="---")
      call read_nc_attribute(ncid, varids(ivar), "molar_mass", molar_mass, default=-999._field_r)
      call read_nc_attribute(ncid, varids(ivar), "lemis", lemis, default=.false.)
      call read_nc_attribute(ncid, varids(ivar), "lreact", lreact, default=.false.)
      call read_nc_attribute(ncid, varids(ivar), "ldep", ldep, default=.false.)
      call read_nc_attribute(ncid, varids(ivar), "lags", lags, default=.false.)

      ! Setup tracer
      call add_tracer(trim(name), long_name=trim(long_name), unit=unit, &
                      molar_mass=molar_mass, lemis=lemis, lreact=lreact, &
                      ldep=ldep, lags=lags, lmicro=.false.)
    end do

    deallocate(varids)

  end subroutine tracer_props_from_netcdf

  !> \brief Read initial profiles for tracers
  !!
  !! \param filename Name of the input file.
  !! \param tracers List of tracers to read input for.
  !! \param svprof 2D array (z,s) to place initial profiles in.
  subroutine tracer_profs_from_netcdf(filename, tracers, nsv, svprof)
    character(*),   intent(in)  :: filename
    type(T_tracer), intent(in)  :: tracers(:)
    real(field_r),  intent(out) :: svprof(:,:)

    integer :: ncid
    integer :: ivar, nsv

    if (.not. file_exists) return

    call nchandle_error(nf90_open(filename, NF90_NOWRITE, ncid))

    nsv = size(tracers)

    do ivar = 1, nsv
      call read_nc_field(ncid, trim(tracers(ivar) % tracname), svprof(:,ivar), &
                         start=1, count=kmax, fillvalue=0._field_r)
    end do

  end subroutine tracer_profs_from_netcdf

end module modtracers
