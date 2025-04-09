module modstat_profiles

  use modglobal,      only: kmax, fname_options, ifnamopt, checknamelisterror, &
                            i1, j1, k1, kmax, ijtot, btime, ih, dt_lim, timee, &
                            tres, rtimee, imax
  use modmpi,         only: d_mpi_bcast, comm3d, mpierr, mpi_iallreduce, &
                            mpi_request, mpi_double, mpi_in_place, mpi_sum
  use modstat_nc
  use modprecision,   only: field_r, longint
  use modslabaverage, only: slabavg

  implicit none

  private

  character(len=*), parameter :: modname = 'modstat_profiles'

  public :: add_profile
  public :: init_profiles
  public :: sample_profiles
  public :: sample_profile
  public :: write_profiles

  interface sample_profile
    module procedure sample_profile
    module procedure sample_profile_masked
  end interface sample_profile

  ! Namelist options
  logical :: lstat
  real    :: dtav, timeav

  ! NetCDF file variables
  character(len=*) :: fname = 'new-profiles.xxx.nc'
  integer          :: ncid, nrec

  ! NetCDF data
  character(len=80), allocatable :: ncname(:,:)
  character(len=80)              :: tncname(1,4)
  real(field_r),     allocatable :: profiles(:,:)

  integer          :: nvar = 0
  integer(longint) :: idtav, itimeav, tnext, tnextwrite
  integer          :: nsamples
  logical          :: do_stats
  logical          :: write_stats

  real(field_r), allocatable :: slab_average(:)

contains

  subroutine add_profile(name, long_name, unit, dim)

    character(len=*), intent(in) :: name
    character(len=*), intent(in) :: long_name
    character(len=*), intent(in) :: unit
    character(len=*), intent(in) :: dim

    character(len=80), allocatable :: tmp_ncname(:,:)

    ! TODO: check for duplicate entries

    ! Allocate array for metadata
    if (.not. allocated(ncname)) then
      allocate(ncname(1,4))
    else
      ! If already allocated, grow in size by 1
      allocate(tmp_ncname(size(ncname, dim=1) + 1, 4))
      tmp_ncname(1:nvar,:) = ncname(1:nvar,:)
      call move_alloc(tmp_ncname, ncname)
    end if

    nvar = nvar + 1

    call ncinfo(ncname(nvar,:), name, long_name, unit, dim)

  end subroutine add_profile

  ! Read namelist and allocate memory
  subroutine init_profiles

    integer :: ierr

    namelist /NAMOUT1D/ lstat, dtav, timeav

    ! Read namelist
    if (myid == 0) then
      open(ifnamopt, file=fname_options, status='old', iostat=ierr)
      read(ifnamopt, NAMOUT1D, iostat=ierr)
      call checknamelisterror(ierr, ifnamopt, 'NAMOUT1D')
      close(ifnamopt)
    end if

    call d_mpi_bcast(lstat, 1, 0, comm3d, mpierr)
    call d_mpi_bcast(dtav, 1, 0, comm3d, mpierr)
    call d_mpi_bcast(timeav, 1, 0, comm3d, mpierr)
    call d_mpi_bcast(nvar, 1, 0, comm3d, mpierr)

    idtav = int(dtav / tres, kind=longint)
    itimeav = int(timeav / tres, kind=longint)

    tnext = idtav + btime
    tnextwrite = itimeav + btime
    nsamples = int(itimeav / idtav)

    if (.not. lstat) return

    allocate(profiles(kmax,nvar))
    allocate(slab_average(1:k1))

    profiles = 0

    ! Initialize the NetCDF file
    if (myid == 0) then
      call nctiminfo(tncname(1,:))
      call open_nc(fname, ncid, nrec, n3=kmax)

      if (nrec == 0) then
        call define_nc(ncid, 1, tncname)
        call writestat_dims_nc(ncid)
      end if

      call define_nc(ncid, nvar, ncname)
    end if

  end subroutine init_profiles

  !> Bookkeeping
  subroutine sample_profiles

    ! Do we need to sample profiles?
    if (timee >= tnext) then
      tnext = tnext + idtav
      do_stats = .true.
    else
      do_stats = .false.
    end if

    ! Do we need to write to file?
    if (timee >= tnextwrite) then
      tnextwrite = tnextwrite + itimeav
      write_stats = .true.
    else
      write_stats = .false.
    end if

    ! Limit time step if we need to sample soon
    ! Note: changing a variable from another module like this is VERY ugly!!
    dt_lim = minval([dt_lim, tnext - timee, tnextwrite - timee])

  end subroutine sample_profiles

  subroutine sample_profile(name, field)

    character(len=*), intent(in) :: name
    real(field_r),    intent(in) :: field(:,:,:)

    character(len=*), parameter :: routine = modname//':sample_profile'

    integer :: k, idx
    integer :: nh
    logical :: found = .false.

    if (.not. do_stats) return

    ! TODO: a hash is probably more efficient here
    do idx = 1, nvar
      if (trim(name) == trim(ncname(idx,1))) then
        found = .true.
        exit
      end if
    end do

    ! A bit hacky maybe: figure out if the given field has ghost cells
    nh = (size(field, dim=1) - imax) / 2

    call slabavg(field, nh, slab_average)

    do k = 1, kmax
      profiles(k,idx) = profiles(k,idx) + slab_average(k)
    end do

  end subroutine sample_profile

  subroutine sample_profile_masked(name, field, mask)

    character(len=*), intent(in) :: name
    real(field_r),    intent(in) :: field(:,:,:)
    logical,          intent(in) :: mask(:,:,:)

    character(len=*), parameter :: routine = modname//':sample_profile_masked'

    integer :: k, idx
    integer :: nh
    logical :: found = .false.

    do idx = 1, nvar
      if (trim(name) == trim(ncname(idx,1))) then
        found = .true.
        exit
      end if
    end do

    nh = (size(field, dim=1) - imax) / 2

    call slabavg(field, mask, nh, slab_average)

    do k = 1, kmax
      profiles(k,idx) = profiles(k,idx) + slab_average(k)
    end do

  end subroutine sample_profile_masked

  subroutine write_profiles

    integer :: k, n
    type(mpi_request) :: requests(nvar)

    if (.not. write_stats) return

    do n = 1, nvar
      do k = 1, kmax
        profiles(k,n) = profiles(k,n) / nsamples
      end do
    end do

    if (myid == 0) then
      call writestat_nc(ncid, 1, tncname, [rtimee], nrec, .true.)
      call writestat_nc(ncid, nvar, ncname, profiles(1:kmax,:), nrec, kmax)
    end if

    ! Reset averages
    do n = 1, nvar
      do k = 1, kmax
        profiles(k,n) = 0
      end do
    end do

  end subroutine write_profiles

end module modstat_profiles