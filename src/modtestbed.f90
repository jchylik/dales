!> \file modtestbed.f90
!!  Testbed continuous forcing & nudging
!>

!>
!!  Testbed continuous forcing & nudging
!>
!!  \author Roel Neggers, IGMK
!!  \par Revision list
!!  \todo Documentation
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
!  Copyright 1993-2009 Delft University of Technology, Wageningen University, Utrecht University, KNMI
!



module modtestbed

use netcdf

implicit none
PRIVATE
PUBLIC :: inittestbed, testbednudge, exittestbed, ltestbed,  & 
          testbed_getinttime, ntnudge, nknudge,              &
          ltb_nudge, ltb_u, ltb_v, ltb_w, ltb_thl, ltb_qt,   &
          tb_time, tb_lat, tb_lon,                           &
          tb_ps,tb_qts,tb_thls,tb_wqs,tb_wts,                &
          tb_z0h, tb_z0m, tb_alb, tb_Qnet,                   &
          tb_seaicefrct, tb_thls_ocean, tb_thls_seaice,      &
          ltb_soildata,                                      &
          ltb_setsurficefrac,                                &
          tb_u,tb_v,tb_w,tb_thl,tb_qt,tb_nccn,tb_ug,tb_vg,   &
          tb_dqtdxls,tb_dqtdyls,                             &
          tb_qtadv,tb_thladv,tb_uadv,tb_vadv,                &
          tb_ql,tb_qi,tb_qladv,tb_qiadv,                     &
          tb_dcl,tb_dci,                                     &  
          tb_sv,tb_svadv,tb_svs,                             &  !#sb3
          ltb_sv,                                            &  !#sb3
          ltb_setccn, ltb_setinp,ltb_clouds,                 &  !#sb3
          ltb_iceinit, ltb_iadv,                             &
          tb_xday,                                           &
          tb_n_inuc, tb_a_inuc, tb_b_inuc,                   &  !#sb3
          tb_n_inuc_r, tb_b_inuc_r,                          &  !#sb3
          ltb_spinnudge, tb_t_spinnudge, tb_tauhigh,         &  !#tb cloud fix
          tb_zmin_spinnudge, tb_zmid_spinnudge,              &
          tb_tsoilav,tb_phiwav,                              &
          tb_land_use,                                       &  !heterogeneous surface conditions
          tbrad_p, tbrad_ql, tbrad_qi, tbrad_qv, tbrad_t, tbrad_o3
SAVE
  real, dimension(:,:,:), allocatable :: tb_sv, tb_svadv  !#sb3
  real, dimension(:,:), allocatable :: tnudge,tb_u,tb_v,tb_w,tb_thl,tb_qt, &
                                       tb_nccn,tb_ql,tb_qi,                &
                                       tb_dcl, tb_dci,                     & 
                                       tb_ug,tb_vg,                        &
                                       tb_dqtdxls,tb_dqtdyls,              &
                                       tb_qtadv,tb_thladv,tb_uadv,tb_vadv, &
                                       tb_qladv,tb_qiadv,                  &
                                       tb_tsoilav,tb_phiwav,               &
                                       tbrad_p, tbrad_t, tbrad_qv,         &
                                       tb_land_use,                        &
                                       tbrad_ql, tbrad_qi, tbrad_o3
  real, dimension(:,:), allocatable :: tb_svs  !#sb3                                  
  real, dimension(:)  , allocatable :: tb_time, tb_lat, tb_lon,            &
                                       tb_ps, tb_qts, tb_thls,             &
                                       tb_wqs, tb_wts,                     &
                                       tb_z0h, tb_z0m, tb_alb, tb_Qnet,    &
                                       tb_seaicefrct,                      &
                                       tb_thls_ocean, tb_thls_seaice
  real, dimension(:)  , allocatable :: upnu,vpnu,wpnu,thlpnu               & ! nudging update profiles
                                      ,qtpnu 

  real :: tb_taunudge = 10800.                &  ! Nudging timescale
         ,tb_zminnudge = 0.                   &  ! Altitude above which is blended nudging applied  #tb
         ,tb_zmidnudge = 0.                   &  ! Altitude above which is full nudging applied     #tb
         ,tb_tauhigh   = 600.                 &  ! strong spinup nudging timescale
         ,tb_t_spinnudge    = 7200.           &  ! strong spinup nudging end 
         ,tb_zmin_spinnudge =    0.           &  ! Altitude above which is strong spinup blended nudging applied
         ,tb_zmid_spinnudge =  200.           &  ! Altitude above which is strong spinup full nudging applied
         ,tb_n_inuc    = 0.                   &  ! inp number parameter #sb3
         ,tb_a_inuc    = 0.                   &  ! inp coefficient a     #sb3
         ,tb_b_inuc    = 0.                   &  ! inp coefficient b #sb3
         ,tb_n_inuc_r  = 0.                   &  ! inp number parameter in Reisner correction  #sb3
         ,tb_b_inuc_r  = 0.                   &  ! inp parameter b in Reisner correction  #sb3
         ,tb_minzinv   = 100.0                &  ! lower bound for inversion search
         ,tb_maxzinv   = 5000.0               &  ! upper bound for inversion search
         ,tb_xday      = -1.0
  logical :: ltestbed   = .false.,            &
             ltb_nudge  = .false.,            &  ! whether to use nudging
             ltb_u      = .true. ,            &  ! nudge in u   (if ltb_nudge=.true.)
             ltb_v      = .true. ,            &  ! nudge in v   (if ltb_nudge=.true.)
             ltb_w      = .true. ,            &  ! nudge in subsidence  (if ltb_nudge=.true.)
             ltb_thl    = .true. ,            &  ! nudge in thl (if ltb_nudge=.true.)
             ltb_qt     = .true. ,            &  ! nudge in qt  (if ltb_nudge=.true.)
             ltb_latlon = .true. ,            &  ! whether use lat, lon from scm_in
             ltb_loadsv = .false.,            &  ! load scalar profiles from testbed #sb3
             ltb_svadv  = .false.,            &  ! load profiles of scalar advection from testbed #sb3
             ltb_svsurf = .false.,            &  ! load scalar surface properties from testbed #sb3
             ltb_sv     = .false.,            &  ! #sb3
             ltb_setccn    = .false.,           &  ! whether n_CCN set #sb3
             ltb_setinp = .true. ,            &  ! set inp properties
             ltb_clouds = .false.,             &  ! liquid cloud advection treated as advection clouds, not just humidity
             ltb_iceinit= .true.,            &  ! ice clouds initialised 
             ltb_iadv   = .true.,            &  ! ice cloud advection 
             ltb_calcsv = .false.,            &  ! additional calculations for scalars #sb3
             ltb_latestart      = .false.,    &  ! whether to start later given by XTIME
             ltb_usedate        = .false.,    &  ! date set based on attribute date
             ltb_spinnudge      =.false.,     &  ! to turn on nudging during spinup
             ltb_soildata       =.true.,      &  ! whether soil data are set 
             ltb_setsurficefrac =.false.         ! whether surfice frac variables are set
  integer :: nknudge,ntnudge

contains
  subroutine inittestbed

    !use modmpi,   only :myid  ,mpierr,comm3d    
    use modmpi,   only :myid,mpierr,comm3d,mpi_logical,mpi_integer &   ! #442fredrik
                      , D_MPI_BCAST                                    ! #442fredrik
    use modglobal,only :ifnamopt,fname_options,k1,&
                        nsv,                      & !#sb3
                        runtime,xlat,xlon,        & ! runtime, longitude and latitude
                        ifinput,cexpnr,           &
                        xtime,                    &
                        xday,                     & 
                        grav,rd,cp,pref0,rlv,zf,checknamelisterror
    use modsurfdata,only : ksoilmax,mpatch
    use modmicrodata3,only: rlvi
    use modforces, only : lforce_user
    ! use modnudge,  only : lnudge         ! #tb #443fredrik
    
    implicit none
    
    real, dimension(:,:,:), allocatable ::  dumsv, dumsv_adv  !#sb3
    real, dimension(:,:), allocatable :: dumomega,dumqv,dumql,dumqi,dumt,dumpf, dumo3,&
                                         dumheight,dumqt,dumnccn,dumthl,dumu,dumv,dumw, &
                                         dumug,dumvg,dumqtadv,dumthladv,dumuadv,dumvadv, &
                                         dumqadv,dumladv,dumiadv,dumtadv, &
                                         dumdcl, dumdci,                  &
                                         dumtsoilav,dumphiwav,dumswi,&
                                         dumlwnet,dumswnet
    real, dimension(:,:), allocatable :: dumsvs   ! #sb3
    real, dimension(:), allocatable :: dumheights

    real :: dumphifc,dumphiwp
    REAL, PARAMETER :: defdcl   = 2.5e-6    &   !< default cloud droplet size
                      ,defdci   = 6.0e-5    &   !< default ice crystal size 
                      ,defnccn  = 1.0e6         !< default initial ccn number
    INTEGER, PARAMETER ::  defnknudges = 1      !< default number of soil levels
    INTEGER  NCID, STATUS, VARID, timID
    INTEGER start2(2), count2(2)
    INTEGER start_sv(3), count_sv(3),start_svs(2), count_svs(2) 
    character(len = nf90_max_name) :: RecordDimName

    integer :: ierr,i,k,ik,ik0,ik1,ikdir,nknudgep1,nknudges
    integer :: i0tb
    real tv,rho,iexner,fac, dumdz, dummax
    real    :: tbtimeoffset, tfac

    character(len=15000) :: readbuffer

    namelist /NAMTESTBED/                   &
        ltestbed, ltb_nudge, tb_taunudge    &
       ,ltb_u, ltb_v, ltb_w                 &
       ,ltb_thl, ltb_qt                     &        
       ,tb_zminnudge                        & !#tb
       ,tb_zmidnudge                        & !#tb
       ,ltb_latestart, ltb_usedate          &
       ,ltb_spinnudge                       &
       ,ltb_clouds, ltb_iceinit, ltb_iadv   & !#tb
       ,tb_zmin_spinnudge,tb_zmid_spinnudge & !#tb
       ,tb_t_spinnudge, tb_tauhigh          &
       ,tb_minzinv, tb_maxzinv              & !#tb
       ,ltb_loadsv, ltb_svsurf, ltb_calcsv    !#sb3
       ! ltb_sv                               !#sb3
       ! ,ltb_latlon                          &

    if(myid==0)then

      open(ifnamopt,file=fname_options,status='old',iostat=ierr)
      read (ifnamopt,NAMTESTBED,iostat=ierr)
      call checknamelisterror(ierr, ifnamopt, 'NAMTESTBED')
      
      write(6 ,NAMTESTBED)
      close(ifnamopt)
      
      ! if (lnudge) then   ! #tb #442fredrik START   
      !   ! lnudge=.false.  ! overriding nudging
      !   write(6,*) "WARNING: set to use both nudging and testbed" !t
      !   write(6,*) "Big PROBLEM ! It can cause double nudging and other issues. "
      !   ! write(6,*) " OVERRIDING, setting lnudge=.false. "
      !endif ! #tb #442fredrik  END

      ! check for
      if (tb_zmidnudge.lt.tb_zminnudge) then
         tb_zmidnudge = tb_zminnudge
      endif
     
      if (nsv>0 .and. ltb_loadsv) then
         ltb_sv = .true.
      endif
      
      if (ltb_calcsv) then   ! #tb START
          write(6,*) "WARNING testbed: ltb_calcsv currently not used" !t
          ltb_calcsv=.false.
      endif ! #tb END
      if (ltb_svsurf) then   ! #tb START
          write(6,*) "WARNING testbed: ltb_svsurf not yet impelemented" !t
          ! ltb_svsurf=.false.
      endif ! #tb END
      if(ltb_svadv) then 
            write(6,*) 'WARNING testbed:','advective tendencies for scalars not yet included'
      endif
      if(ltb_sv.and. ltb_nudge) then 
            write(6,*) 'WARNING testbednudge:','nudging for scalars not yet included'
      endif
     
    end if
 
    call D_MPI_BCAST(ltestbed     , 1  ,0,comm3d,mpierr)
    call D_MPI_BCAST(ltb_nudge    , 1  ,0,comm3d,mpierr)
    call D_MPI_BCAST(ltb_u        , 1  ,0,comm3d,mpierr)
    call D_MPI_BCAST(ltb_v        , 1  ,0,comm3d,mpierr)
    call D_MPI_BCAST(ltb_w        , 1  ,0,comm3d,mpierr)
    call D_MPI_BCAST(ltb_thl      , 1  ,0,comm3d,mpierr)
    call D_MPI_BCAST(ltb_qt       , 1  ,0,comm3d,mpierr)
    ! call D_MPI_BCAST(ltb_latlon   , 1  ,0,comm3d,mpierr)
    call D_MPI_BCAST(tb_taunudge  , 1      ,0,comm3d,mpierr)
    call D_MPI_BCAST(tb_zminnudge , 1      ,0,comm3d,mpierr) !#tb
    call D_MPI_BCAST(tb_zmidnudge , 1      ,0,comm3d,mpierr) !#tb
    call D_MPI_BCAST(ltb_sv       , 1  ,0,comm3d,mpierr)
    call D_MPI_BCAST(ltb_svadv    , 1  ,0,comm3d,mpierr)
    call D_MPI_BCAST(ltb_svsurf   , 1  ,0,comm3d,mpierr)
    call D_MPI_BCAST(ltb_calcsv   , 1  ,0,comm3d,mpierr)
    call D_MPI_BCAST(ltb_spinnudge, 1  ,0,comm3d,mpierr)
    call D_MPI_BCAST(ltb_clouds   , 1  ,0,comm3d,mpierr)
    call D_MPI_BCAST(ltb_iadv     , 1  ,0,comm3d,mpierr)
    call D_MPI_BCAST(ltb_iceinit  , 1  ,0,comm3d,mpierr)
    call D_MPI_BCAST(ltb_usedate  , 1  ,0,comm3d,mpierr)
    call D_MPI_BCAST(tb_tauhigh       , 1      ,0,comm3d,mpierr)
    call D_MPI_BCAST(tb_t_spinnudge   , 1      ,0,comm3d,mpierr)     
    call D_MPI_BCAST(tb_zmin_spinnudge , 1      ,0,comm3d,mpierr) !#tb
    call D_MPI_BCAST(tb_zmid_spinnudge , 1      ,0,comm3d,mpierr) !#tb
    call D_MPI_BCAST(tb_minzinv   , 1      ,0,comm3d,mpierr) !#tb
    call D_MPI_BCAST(tb_maxzinv   , 1      ,0,comm3d,mpierr) !#tb
    
    if (.not. ltestbed) return
    
    ! lforce_user = .true.

    if(myid==0) then

        !--- open nc file ---
        STATUS = NF90_OPEN('scm_in.nc', nf90_nowrite, NCID)
        if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
          
          
        !--- get time & height dimensions ---
        status = nf90_inq_dimid(ncid, "time", timID)
        if (status /= nf90_noerr) call handle_err(status)
        status = nf90_inquire_dimension(NCID, timID, len=ntnudge, name=RecordDimName)
        if (STATUS .ne. nf90_noerr) call handle_err(STATUS)

!        write(6,'(a15,i10," ",a10)') 'scm_in time:',ntnudge,RecordDimName

        status = nf90_inq_dimid(ncid, "nlev", timID)
        if (status /= nf90_noerr) call handle_err(status)
        status = nf90_inquire_dimension(NCID, timID, len=nknudge, name=RecordDimName)
        if (STATUS .ne. nf90_noerr) call handle_err(STATUS)

        status = nf90_inq_dimid(ncid, "nlevp1", timID)
        if (status /= nf90_noerr) call handle_err(status)
        status = nf90_inquire_dimension(NCID, timID, len=nknudgep1, name=RecordDimName)
        if (STATUS .ne. nf90_noerr) call handle_err(STATUS)

        status = nf90_inq_dimid(ncid, "nlevs", timID)
        if (status .ne. nf90_noerr) then !call handle_err(status)
          write(6,*) '  inittestbed:  scm_in missing "nlevs", dimension of soil data, setting to ' &
                     ,defnknudges 
          nknudges = defnknudges
        else
          status = nf90_inquire_dimension(NCID, timID, len=nknudges, name=RecordDimName)
          if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
        endif

        write(6,*) '  inittestbed: scm_in: dimensions:',ntnudge, nknudge, nknudgep1, nknudges

    end if

    call D_MPI_BCAST(ntnudge    , 1  ,0,comm3d,mpierr)
    call D_MPI_BCAST(nknudge    , 1  ,0,comm3d,mpierr)
    call D_MPI_BCAST(nknudgep1  , 1  ,0,comm3d,mpierr)
    call D_MPI_BCAST(nknudges   , 1  ,0,comm3d,mpierr)

    !--- allocate space for input variables & reset---
    allocate(    tnudge    (ntnudge,k1), &
                 tb_u      (ntnudge,k1), &
                 tb_v      (ntnudge,k1), &
                 tb_w      (ntnudge,k1), &
                 tb_thl    (ntnudge,k1), &
                 tb_qt     (ntnudge,k1), &
                 tb_nccn   (ntnudge,k1), &
                 tb_ql     (ntnudge,k1), &
                 tb_qi     (ntnudge,k1), &
                 tb_dcl    (ntnudge,k1), &
                 tb_dci    (ntnudge,k1), &
                 tb_ug     (ntnudge,k1), &
                 tb_vg     (ntnudge,k1), &
                 tb_dqtdxls(ntnudge,k1), &
                 tb_dqtdyls(ntnudge,k1), &
                 tb_qtadv  (ntnudge,k1), &
                 tb_qladv  (ntnudge,k1), &
                 tb_qiadv  (ntnudge,k1), &
                 tb_thladv (ntnudge,k1), &
                 tb_uadv   (ntnudge,k1), &
                 tb_vadv   (ntnudge,k1) )

    allocate(    tb_time   (ntnudge), &
                 tb_lat    (ntnudge), &
                 tb_lon    (ntnudge), &
                 tb_ps     (ntnudge), &
                 tb_qts    (ntnudge), &
                 tb_thls   (ntnudge), &
                 tb_wts    (ntnudge), &
                 tb_wqs    (ntnudge), &
                 tb_z0m    (ntnudge), &
                 tb_z0h    (ntnudge), &
                 tb_alb    (ntnudge), &
                 tb_Qnet   (ntnudge), &
                 tb_seaicefrct (ntnudge), &
                 tb_thls_ocean (ntnudge), &
                 tb_thls_seaice(ntnudge) )

    allocate(    tb_tsoilav(ntnudge,ksoilmax), &
                 tb_phiwav (ntnudge,ksoilmax) )

    allocate(    tb_land_use(mpatch,mpatch))

    allocate(    tbrad_p    (ntnudge, nknudge), &
                 tbrad_t    (ntnudge, nknudge), &
                 tbrad_qv   (ntnudge, nknudge), &
                 tbrad_ql   (ntnudge, nknudge), &
                 tbrad_qi   (ntnudge, nknudge), &
                 tbrad_o3   (ntnudge, nknudge) )

    if (ltb_nudge) then   ! allocate nudging update arrays
      allocate( upnu   (k1)                     &
               ,vpnu   (k1)                     &
               ,wpnu   (k1)                     &
               ,thlpnu (k1)                     &
               ,qtpnu  (k1)                     &
      )
      upnu=0;vpnu=0;wpnu=0;thlpnu=0;qtpnu=0
    endif ! ltb_nudge
    if (nsv>0) then               
     allocate(   tb_sv      (ntnudge,k1,nsv),   &  ! #sb3 
                 tb_svadv   (ntnudge,k1,nsv),   &  ! #sb3
                 tb_svs     (ntnudge,nsv)       &  ! #sb3
             )
        tb_sv = 0     !#sb3
        tb_svadv = 0  !#sb3
        tb_svs = 0    !#sb3
    endif  

     tnudge = tb_taunudge     !nudging timescale

     tb_time=0;tb_lat=0;tb_lon=0;tb_ps=0;tb_qts=0;tb_thls=0;tb_wts=0;tb_wqs=0
     tb_z0m=0;tb_z0h=0;tb_alb=0;tb_Qnet=0; tb_seaicefrct=0;tb_thls_ocean=0;tb_thls_seaice=0
     tb_u=0;tb_v=0;tb_w=0;tb_thl=0;tb_qt=0;tb_nccn=0;tb_ql=0;tb_qi=0;tb_ug=0;tb_vg=0
     tb_dcl = 1.0e-6 ; tb_dci = 1.0e-6
     tb_dqtdxls=0;tb_dqtdyls=0;tb_qtadv=0;tb_thladv=0;tb_uadv=0;tb_vadv=0
     tb_qladv=0;tb_qiadv=0;
     tb_tsoilav=0;tb_phiwav=0
     tbrad_t=0;tbrad_qv =0;tbrad_ql =0;tbrad_qi =0;tbrad_o3 =0
     tb_land_use = 0
     tb_n_inuc = 0; tb_b_inuc_r =0; 

    
    if(myid==0) then

        allocate(dumomega (nknudge,ntnudge), &
                 dumheight(nknudge,ntnudge), &
                 dumpf    (nknudge,ntnudge), &
                 dumqv    (nknudge,ntnudge), &
                 dumql    (nknudge,ntnudge), &
                 dumqi    (nknudge,ntnudge), &
                 dumdcl   (nknudge,ntnudge), &
                 dumdci   (nknudge,ntnudge), &
                 dumo3    (nknudge,ntnudge), &
                 dumt     (nknudge,ntnudge), &
                 dumqt    (nknudge,ntnudge), &
                 dumnccn  (nknudge,ntnudge), &
                 dumthl   (nknudge,ntnudge), &
                 dumu     (nknudge,ntnudge), &
                 dumv     (nknudge,ntnudge), &
                 dumw     (nknudge,ntnudge), &
                 dumug    (nknudge,ntnudge), &
                 dumvg    (nknudge,ntnudge), &
                 dumqtadv (nknudge,ntnudge), &
                 dumthladv(nknudge,ntnudge), &
                 dumuadv  (nknudge,ntnudge), &
                 dumvadv  (nknudge,ntnudge), &
                 dumqadv  (nknudge,ntnudge), &
                 dumladv  (nknudge,ntnudge), &
                 dumiadv  (nknudge,ntnudge), & 
                 dumtadv  (nknudge,ntnudge)  & 
                )
        allocate(dumswnet  (nknudgep1,ntnudge), & 
                 dumlwnet  (nknudgep1,ntnudge), & 
                 dumheights (nknudges), & 
                 dumtsoilav (nknudges,ntnudge), & 
                 dumphiwav  (nknudges,ntnudge), & 
                 dumswi     (nknudges,ntnudge) &
                 )
        !if (nsv>0 .and. ltb_loadsv) then
        !if (ltb_sv) then  
        !   allocate(dumsv    (nsv,nknudge,ntnudge) ) ! 3D fields
        !   allocate(dumsv_adv(nsv,nknudge,ntnudge) ) ! add dumsv_adv
        !   allocate(dumsvs   (nsv,ntnudge) )         ! surface values 
        !endif 


        !--- timeseries ---

        !  time
        STATUS = NF90_INQ_VARID(NCID, 'time', VARID)
        if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
        STATUS = NF90_GET_VAR (NCID, VARID, tb_time, start=(/1/), count=(/ntnudge/) )
        if (STATUS .ne. nf90_noerr) call handle_err(STATUS) 
        !        write(6,'(a30,5f10.2)') 'inittestbed: tb_time:',&
        !             tb_time(1),tb_time(2),tb_time(3),tb_time(ntnudge-1),tb_time(ntnudge)
        ! 
        ! later start (or earlier)
        if (ltb_latestart) then  
           tbtimeoffset = 3600.*xtime - tb_time(1) ! offset
           write(6,*) "message testbed: shifting start by xtime = ", xtime," hours"
        else 
           tbtimeoffset = 0.0
        endif
        ! check runtime length 
        if (maxval(tb_time)<(runtime-tbtimeoffset)) then
           write(6,*) 'WARNING: highest time in scm_in.nc is less than runtime !' 
           runtime = maxval(tb_time)-1
           write(6,*) ' FIX: shortening runtime to maxval(time)'
           write(6,*) ' RUNTIME=  ' ,runtime
        endif

       ! lat and lon from scm_in if possible: 
        STATUS = NF90_INQ_VARID(NCID, 'lat', VARID)
        if (STATUS .ne. nf90_noerr) then
            write(6,*) 'WARNING: lat not found in scm_in.nc'
            ltb_latlon = .false.
        endif     
        STATUS = NF90_INQ_VARID(NCID, 'lon', VARID)
        if (STATUS .ne. nf90_noerr) then
            write(6,*) 'WARNING: lon not found in scm_in.nc'
            ltb_latlon = .false.
        endif       
        if (ltb_latlon) then
          !  latitude
          STATUS = NF90_INQ_VARID(NCID, 'lat', VARID)
          if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
          STATUS = NF90_GET_VAR (NCID, VARID, tb_lat, start=(/1/), count=(/ntnudge/) )
          if (STATUS .ne. nf90_noerr) call handle_err(STATUS) 
          !  longitude
          STATUS = NF90_INQ_VARID(NCID, 'lon', VARID)
          if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
          STATUS = NF90_GET_VAR (NCID, VARID, tb_lon, start=(/1/), count=(/ntnudge/) )
          if (STATUS .ne. nf90_noerr) call handle_err(STATUS) 
        
        else ! ltb_latlon = .false. 
          ! instead use modglobal one 
          tb_lat = xlat
          tb_lon = xlon
          write(6,*) 'WARNING: lat or lon not found in scm_in.nc'
          write(6,*) ' FIX: using xlat and xlon defined in namelist DOMAIN'
        endif

        !  surface pressure
        STATUS = NF90_INQ_VARID(NCID, 'ps', VARID)
        if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
        STATUS = NF90_GET_VAR (NCID, VARID, tb_ps, start=(/1/), count=(/ntnudge/) )
        if (STATUS .ne. nf90_noerr) call handle_err(STATUS) 

        !  surface temperature
        STATUS = NF90_INQ_VARID(NCID, 't_skin', VARID)
        if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
        STATUS = NF90_GET_VAR (NCID, VARID, tb_thls, start=(/1/), count=(/ntnudge/) )
        if (STATUS .ne. nf90_noerr) call handle_err(STATUS) 

        !  surface humidity
!        STATUS = NF90_INQ_VARID(NCID, '', VARID)
!        if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
!        STATUS = NF90_GET_VAR (NCID, VARID, tb_qts, start=(/1/), count=(/ntnudge/) )
!        if (STATUS .ne. nf90_noerr) call handle_err(STATUS) 

        !  surface T flux
        STATUS = NF90_INQ_VARID(NCID, 'sfc_sens_flx', VARID)
        if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
        STATUS = NF90_GET_VAR (NCID, VARID, tb_wts, start=(/1/), count=(/ntnudge/) )
        if (STATUS .ne. nf90_noerr) call handle_err(STATUS) 

        !  surface q flux
        STATUS = NF90_INQ_VARID(NCID, 'sfc_lat_flx', VARID)
        if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
        STATUS = NF90_GET_VAR (NCID, VARID, tb_wqs, start=(/1/), count=(/ntnudge/) )
        if (STATUS .ne. nf90_noerr) call handle_err(STATUS) 

        !  roughness length for momentum
        STATUS = NF90_INQ_VARID(NCID, 'mom_rough', VARID)
        if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
        STATUS = NF90_GET_VAR (NCID, VARID, tb_z0m, start=(/1/), count=(/ntnudge/) )
        if (STATUS .ne. nf90_noerr) call handle_err(STATUS) 
          
        !  roughness length for heat and moisture
        STATUS = NF90_INQ_VARID(NCID, 'heat_rough', VARID)
        if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
        STATUS = NF90_GET_VAR (NCID, VARID, tb_z0h, start=(/1/), count=(/ntnudge/) )
        if (STATUS .ne. nf90_noerr) call handle_err(STATUS) 
          
        !  surface albedo, for radiation
        STATUS = NF90_INQ_VARID(NCID, 'albedo', VARID)
        if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
        STATUS = NF90_GET_VAR (NCID, VARID, tb_alb, start=(/1/), count=(/ntnudge/) )
        if (STATUS .ne. nf90_noerr) call handle_err(STATUS) 
        
        ! Sea Ice fraction when defined
        STATUS = NF90_INQ_VARID(NCID, 'sea_ice_frct', VARID)
        if (STATUS .ne. nf90_noerr) then
          write(6,*) "NOTICE testbed: sea_ice_frct not in scm_in file"
          ltb_setsurficefrac = .false. ! surficefrac variables cannot be set from modtestbed
        else
          STATUS = NF90_GET_VAR (NCID, VARID, tb_seaicefrct, start=(/1/), count=(/ntnudge/) )
          if (STATUS .ne. nf90_noerr) call handle_err(STATUS) 
          write(6,*) "message testbed: sea_ice_frct loaded from scm_in file"
          ltb_setsurficefrac = .true. ! if following variables are also set
        endif

        ! skin temperatures over open ocean and sea ice when defined
        STATUS = NF90_INQ_VARID(NCID, 't_skin_ocean', VARID)
        if (STATUS .ne. nf90_noerr) then
          write(6,*) "NOTICE testbed: t_skin_ocean not in scm_in file"
          ltb_setsurficefrac = .false. ! surficefrac variables cannot be set from modtestbed
        else
          STATUS = NF90_GET_VAR (NCID, VARID, tb_thls_ocean, start=(/1/), count=(/ntnudge/) )
          if (STATUS .ne. nf90_noerr) call handle_err(STATUS) 
          write(6,*) "message testbed: t_skin_ocean loaded from scm_in file"
        endif
        STATUS = NF90_INQ_VARID(NCID, 't_skin_seaice', VARID)
        if (STATUS .ne. nf90_noerr) then
          write(6,*) "NOTICE testbed: t_skin_seaice not in scm_in file"
          ltb_setsurficefrac = .false. ! surficefrac variables cannot be set from modtestbed
        else
          STATUS = NF90_GET_VAR (NCID, VARID, tb_thls_seaice, start=(/1/), count=(/ntnudge/) )
          if (STATUS .ne. nf90_noerr) call handle_err(STATUS) 
          write(6,*) "message testbed: t_skin_seaice loaded from scm_in file"
        endif

        do i=1,ntnudge
        
          rho = tb_ps(i) / (rd * tb_thls(i))
          tb_wts(i) = -tb_wts(i) / (cp  * rho)        !Change sign: upward = positive in LES, but by convention upward = negative in most GCMs.
          tb_wqs(i) = -tb_wqs(i) / (rlv * rho)
          iexner = (tb_ps(i)/pref0)**(-rd/cp) 

          tb_thls       (i) = iexner * tb_thls       (i)
          tb_thls_ocean (i) = iexner * tb_thls_ocean (i)
          tb_thls_seaice(i) = iexner * tb_thls_seaice(i)

        end do


        !--- profiles full levels ---
        start2 = (/ 1      , 1       /)
        count2 = (/ nknudge, ntnudge /)

        !  height
        STATUS = NF90_INQ_VARID(NCID, 'height_f', VARID)
        if (STATUS .ne. nf90_noerr) then
          STATUS = NF90_INQ_VARID(NCID, 'zf', VARID)
          if (STATUS .ne. nf90_noerr) then
            write(6,*) 'height array not found in scm_in.nc'
            call handle_err(STATUS)
          else
            write(6,*) 'modetestbed: loading variable zf (KPT naming)' !#lagtraj
          endif
        endif
        STATUS = NF90_GET_VAR (NCID, VARID, dumheight, start=start2, count=count2)
        if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
         
        !  pressure
        STATUS = NF90_INQ_VARID(NCID, 'pressure_f', VARID)
        if (STATUS .ne. nf90_noerr) then
          STATUS = NF90_INQ_VARID(NCID, 'pres', VARID)
          if (STATUS .ne. nf90_noerr) then
            write(6,*) 'pressure array not found in scm_in.nc'
            call handle_err(STATUS)
          else
            write(6,*) 'modetestbed: loading variable pres (KPT naming)' !#lagtraj
          endif
        endif
        STATUS = NF90_GET_VAR (NCID, VARID, dumpf, start=start2, count=count2)
        if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
        
        !  u
        !STATUS = NF90_INQ_VARID(NCID, 'u', VARID)
        STATUS = NF90_INQ_VARID(NCID, 'u_local', VARID)
        if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
        STATUS = NF90_GET_VAR (NCID, VARID, dumu, start=start2, count=count2)
        if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
          
        !  v
        !STATUS = NF90_INQ_VARID(NCID, 'v', VARID)
        STATUS = NF90_INQ_VARID(NCID, 'v_local', VARID)
        if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
        STATUS = NF90_GET_VAR (NCID, VARID, dumv, start=start2, count=count2)
        if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
          
        !  qt
        !STATUS = NF90_INQ_VARID(NCID, 'q', VARID)
        STATUS = NF90_INQ_VARID(NCID, 'q_local', VARID)
        if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
        STATUS = NF90_GET_VAR (NCID, VARID, dumqv, start=start2, count=count2)
        if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
        
        !  ql
        !STATUS = NF90_INQ_VARID(NCID, 'ql', VARID)
        STATUS = NF90_INQ_VARID(NCID, 'ql_local', VARID)
        if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
        STATUS = NF90_GET_VAR (NCID, VARID, dumql, start=start2, count=count2)
        if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
        
        !  qi
        !STATUS = NF90_INQ_VARID(NCID, 'qi', VARID)
        STATUS = NF90_INQ_VARID(NCID, 'qi_local', VARID)
        if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
        STATUS = NF90_GET_VAR (NCID, VARID, dumqi, start=start2, count=count2)
        if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
        
        ! D_cl
        STATUS = NF90_INQ_VARID(NCID, 'd_cl', VARID)
        if (STATUS .ne. nf90_noerr) then
          write(6,*) 'WARNING: d_cl not in scm_in.nc, setting tb_dcl to default value'
          dumdcl = defdcl ! 1.0e-7
          ! 
        else
          STATUS = NF90_GET_VAR (NCID, VARID, dumdcl, start=start2, count=count2)
          if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
          ! 
        endif 
        
        ! D_ci
        STATUS = NF90_INQ_VARID(NCID, 'd_ci', VARID)
        if (STATUS .ne. nf90_noerr) then
          write(6,*) 'WARNING: d_ci not in scm_in.nc, setting tb_dci to default value'
          dumdci = defdci ! 1.0e-6
          ! 
        else
          STATUS = NF90_GET_VAR (NCID, VARID, dumdci, start=start2, count=count2)
          if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
          ! 
        endif   
        
        
        ! n_ccn
        STATUS = NF90_INQ_VARID(NCID, 'n_ccn', VARID)
        if (STATUS .ne. nf90_noerr) then
          write(6,*) 'WARNING: n_ccn not in scm_in.nc'
          dumnccn = defnccn ! 0.
          ltb_setccn = .false.
        else
          STATUS = NF90_GET_VAR (NCID, VARID, dumnccn, start=start2, count=count2)
          if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
          ltb_setccn = .true.
        endif
          
        !  thl
        !STATUS = NF90_INQ_VARID(NCID, 't', VARID)
        STATUS = NF90_INQ_VARID(NCID, 't_local', VARID)
        if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
        STATUS = NF90_GET_VAR (NCID, VARID, dumt, start=start2, count=count2)
        if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
        
        !if (ltb_iceinit) then
         do i=1,ntnudge
         do k=1,nknudge
          dumqt(k,i) = dumqv(k,i) + dumql(k,i) ! ice treated separately
          iexner = (dumpf(k,i)/pref0)**(-rd/cp)
          ! ice treated separately 
          ! theta_l corrected by cloud liquid water content 
          dumthl(k,i) = dumt(k,i) * iexner                   & 
                       - rlv * dumql(k,i)* iexner / cp
         enddo
         enddo        
        !else  ! ltb_iceinit = .false.
        ! do i=1,ntnudge
        ! do k=1,nknudge
        !  dumqt(k,i) = dumqv(k,i) + dumql(k,i) + dumqi(k,i)
        !  iexner = (dumpf(k,i)/pref0)**(-rd/cp)
        !  ! theta_l corrected by cloud liquid water content and cloud ice
        !  ! multiplied by rlvi to account for heat that could be release from freezing
        !  dumthl(k,i) = dumt(k,i) * iexner                    &
        !               - rlv * dumql(k,i)* iexner / cp        &
        !               - rlvi *dumqi(k,i)* iexner / cp  
        ! enddo
        ! enddo
        !endif
        !do i=1,ntnudge
        !do k=1,nknudge
        !  iexner = (dumpf(k,i)/pref0)**(-rd/cp)
        !  dumthl(k,i) = dumt(k,i) * iexner - rlv * dumql(k,i)* iexner / cp
        !enddo
        !enddo

        !  O3
        STATUS = NF90_INQ_VARID(NCID, 'o3', VARID)
        if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
        STATUS = NF90_GET_VAR (NCID, VARID, dumo3, start=start2, count=count2)
        if (STATUS .ne. nf90_noerr) call handle_err(STATUS)

        !  w
        STATUS = NF90_INQ_VARID(NCID, 'omega', VARID)
        if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
        STATUS = NF90_GET_VAR (NCID, VARID, dumomega, start=start2, count=count2)
        if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
                
        do i=1,ntnudge
        do k=1,nknudge
          tv  = dumt(k,i) * (1.+0.61*dumqv(k,i))
          rho = dumpf(k,i) / (rd*tv)
          dumw(k,i) = - dumomega(k,i) / ( rho * grav )     !convert from Pa/s to m/s
        enddo
        enddo

        !  ug
        STATUS = NF90_INQ_VARID(NCID, 'ug', VARID)
        if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
        STATUS = NF90_GET_VAR (NCID, VARID, dumug, start=start2, count=count2)
        if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
          
        !  vg
        STATUS = NF90_INQ_VARID(NCID, 'vg', VARID)
        if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
        STATUS = NF90_GET_VAR (NCID, VARID, dumvg, start=start2, count=count2)
        if (STATUS .ne. nf90_noerr) call handle_err(STATUS)          
         
        !  uadv
        STATUS = NF90_INQ_VARID(NCID, 'uadv', VARID)
        if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
        STATUS = NF90_GET_VAR (NCID, VARID, dumuadv, start=start2, count=count2)
        if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
          
        !  vadv
        STATUS = NF90_INQ_VARID(NCID, 'vadv', VARID)
        if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
        STATUS = NF90_GET_VAR (NCID, VARID, dumvadv, start=start2, count=count2)
        if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
          
        !  qtadv
        STATUS = NF90_INQ_VARID(NCID, 'qadv', VARID)
        if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
        STATUS = NF90_GET_VAR (NCID, VARID, dumqadv, start=start2, count=count2)
        if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
        
        !  qladv
        STATUS = NF90_INQ_VARID(NCID, 'ladv', VARID)
        if (STATUS .ne. nf90_noerr) then
          STATUS = NF90_INQ_VARID(NCID, 'qladv', VARID)
          if (STATUS .ne. nf90_noerr) then
            write(6,*) 'modtestbed: neither ldav nor qladv found in scm_in.nc'
            dumladv = 0.
            write(6,*) 'modtestbed: => setting advection of liquid clouds to 0. '
          else
            write(6,*) 'modetestbed: loading variable qladv (KPT naming)' !#lagtraj
            STATUS = NF90_GET_VAR (NCID, VARID, dumiadv, start=start2, count=count2)
            if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
          endif          
        else
          STATUS = NF90_GET_VAR (NCID, VARID, dumladv, start=start2, count=count2)
          if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
        endif
        
        !  qiadv
        STATUS = NF90_INQ_VARID(NCID, 'iadv', VARID)
        if (STATUS .ne. nf90_noerr) then
          STATUS = NF90_INQ_VARID(NCID, 'qiadv', VARID)
          if (STATUS .ne. nf90_noerr) then
            write(6,*) 'modtestbed: neither idav nor qiadv found in scm_in.nc'
            dumiadv = 0.
            write(6,*) 'modtestbed: => setting advection of ice clouds to 0. '
          else
            write(6,*) 'modetestbed: loading variable qiadv (KPT naming)' !#lagtraj
            STATUS = NF90_GET_VAR (NCID, VARID, dumiadv, start=start2, count=count2)
            if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
          endif
        else
          STATUS = NF90_GET_VAR (NCID, VARID, dumiadv, start=start2, count=count2)
          if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
        endif

          
        !  thladv
        STATUS = NF90_INQ_VARID(NCID, 'tadv', VARID)
        if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
        STATUS = NF90_GET_VAR (NCID, VARID, dumtadv, start=start2, count=count2)
        if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
        
        !if (ltb_iadv) then
         do i=1,ntnudge
         do k=1,nknudge
           dumqtadv(k,i) = dumqadv(k,i) + dumladv(k,i) ! ice treated separately
           iexner = (dumpf(k,i)/pref0)**(-rd/cp)
           ! decreasing heat advection by excess advected humidity
           ! that is likely to turn into liquid droplets
           dumthladv(k,i) = dumtadv(k,i) * iexner               &
                            - rlv * dumladv(k,i)*iexner/cp
         enddo
         enddo        
        !else ! ltb_iadv = .false. 
        ! do i=1,ntnudge
        ! do k=1,nknudge
        !   dumqtadv(k,i) = dumqadv(k,i) + dumladv(k,i) + dumiadv(k,i)
        !   iexner = (dumpf(k,i)/pref0)**(-rd/cp)
        !   ! decreasing heat advection by excess advected humidity
        !   ! that is likely to turn into ice or liquid droplets
        !   dumthladv(k,i) = dumtadv(k,i) * iexner               &
        !                    - rlv * dumladv(k,i)*iexner/cp      &
        !                    - rlvi* dumiadv(k,i)*iexner/cp 
        ! enddo
        ! enddo
        !endif
        !do i=1,ntnudge
        !do k=1,nknudge
        !  iexner = (dumpf(k,i)/pref0)**(-rd/cp)
        !  !dumthladv(k,i) = dumtadv(k,i) * iexner
        !  dumthladv(k,i) = dumtadv(k,i) * iexner - rlv*dumladv(k,i)*iexner/cp
        !enddo
        !enddo

        !read land_use from txt file

        open(ifinput, file='land_use.inp.'//cexpnr)
        ierr = 0
        do while (ierr == 0)
                read(ifinput,'(A)',iostat=ierr) readbuffer  
                if (ierr == -1) then
                        if (myid == 0) then
                                print *, "Modtestbed: No land_use file. Right behaviour if (lhetero .and. ltestbed) == false"
                        endif
                endif
                if (ierr == 0) then
                        if (readbuffer(1:1) == '#') then
                                if (myid == 0) then
                                        print * , trim(readbuffer)
                                endif
                        else
                                read(readbuffer, *, iostat=ierr) tb_land_use
                                exit !only read one line
                        endif
                endif
        enddo
        close(ifinput)
       

        !--- profiles half levels ---
        start2 = (/ 1        , 1       /)
        count2 = (/ nknudgep1, ntnudge /)
          
        !  net SW downward flux
        STATUS = NF90_INQ_VARID(NCID, 'fradSWnet', VARID)
        if (STATUS .ne. nf90_noerr) then
          write(6,*) 'fradSWnet not in scm_in.nc'
          dumswnet = 0.
        else
          STATUS = NF90_GET_VAR (NCID, VARID, dumswnet, start=start2, count=count2)
          if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
        endif
          
        !  net LW downward flux
        STATUS = NF90_INQ_VARID(NCID, 'fradLWnet', VARID)
        if (STATUS .ne. nf90_noerr) then
          write(6,*) 'fradLWnet not in scm_in.nc'
          dumswnet = 0.
        else
          STATUS = NF90_GET_VAR (NCID, VARID, dumlwnet, start=start2, count=count2)
          if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
        endif
          
        do i=1,ntnudge
          tb_Qnet(i) = dumswnet(nknudgep1,i) + dumlwnet(nknudgep1,i)      !flux at surface is stored in lowest half level of profile
!          write(6,*) "modtestbed: qnet:",i,tb_Qnet(i),nknudge,nknudgep1
        enddo


        !--- soil profiles ---
        
        STATUS = NF90_INQ_VARID(NCID, 'h_soil', VARID)
        if (STATUS .ne. nf90_noerr) then ! call handle_err(STATUS)
          write(6,*) 'h_soil not in scm_in.nc'
          write(6,*) 'land-surface data cannot be loaded due to absence of h_soil'
          ltb_soildata = .false.  ! lack of soil data does not exit sim 
        else
          STATUS = NF90_GET_VAR (NCID, VARID, dumheights, start=(/1/), count=(/nknudges/) )
          if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
          ltb_soildata = .true.
        endif 

        start2 = (/ 1       , 1       /)
        count2 = (/ nknudges, ntnudge /)

        if (ltb_soildata) then  ! i.e. to skip loading soil data if missing
        
          !  tsoilav
          STATUS = NF90_INQ_VARID(NCID, 't_soil', VARID)
          if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
          STATUS = NF90_GET_VAR (NCID, VARID, dumtsoilav, start=start2, count=count2)
          if (STATUS .ne. nf90_noerr) call handle_err(STATUS)

          !  phiwav
          STATUS = NF90_INQ_VARID(NCID, 'q_soil', VARID)
          if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
          STATUS = NF90_GET_VAR (NCID, VARID, dumphiwav, start=start2, count=count2)
          if (STATUS .ne. nf90_noerr) call handle_err(STATUS)


          !--- soil scalars ---
        
          !  field capacity
          status = nf90_inquire_attribute(ncid, nf90_global, "field_capacity")
          if (status /= nf90_noerr) call handle_err(status)
          status = nf90_get_att(ncid, nf90_global, "field_capacity", dumphifc)
          if (status /= nf90_noerr) call handle_err(status)
        
          !  wilting point
          status = nf90_inquire_attribute(ncid, nf90_global, "wilting_point")
          if (status /= nf90_noerr) call handle_err(status)
          status = nf90_get_att(ncid, nf90_global, "wilting_point", dumphiwp)
          if (status /= nf90_noerr) call handle_err(status)
          
          dumswi = ( dumphiwav - dumphiwp ) / ( dumphifc - dumphiwp )   !soil wetness index, using input values for wilting point and field capacity
          ! (dumswi calculation is avoided if data are not set)
        endif
        
        !--- get the date 
        if (ltb_usedate) then
          STATUS = nf90_inquire_attribute(ncid, nf90_global, "dayofyear")
          if (STATUS .ne. nf90_noerr) then
            write(6,*) 'ERROR modtestbed: dayofyear not in scm_in.nc'
            write(6,*) '  ltb_usedate=.TRUE. requires dayofyear in scm_in.nc'
            call handle_err(STATUS)  
          else
            STATUS = nf90_get_att(ncid, nf90_global, "dayofyear", tb_xday)
            if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
            write(6,*) 'inittestbed: dayofyear loaded from scm_in'
            write(6,*) 'inittestbed: setting xday=', tb_xday
          endif        
        endif

        
        !--- INP settings ---
        !  number parameter
        STATUS = nf90_inquire_attribute(ncid, nf90_global, "in_n_inuc")
        if (STATUS .ne. nf90_noerr) then
          write(6,*) 'NOTICE modtestbed: in_n_inuc not in scm_in.nc'
          ltb_setinp = .false.
        else
          STATUS = nf90_get_att(ncid, nf90_global, "in_n_inuc", tb_n_inuc)
          if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
        endif
        
        !  coefficient a
        STATUS = nf90_inquire_attribute(ncid, nf90_global, "in_a_inuc")
        if (STATUS .ne. nf90_noerr) then
          write(6,*) 'NOTICE modtestbed: in_b_inuc not in scm_in.nc'
          ltb_setinp  = .false.      
        else
          STATUS = nf90_get_att(ncid, nf90_global, "in_a_inuc", tb_a_inuc)
          if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
        endif
        ! coefficient b
        STATUS = nf90_inquire_attribute(ncid, nf90_global, "in_b_inuc")
        if (STATUS .ne. nf90_noerr) then
          write(6,*) 'NOTICE modtestbed: in_b_inuc not in scm_in.nc'
          ltb_setinp  = .false.      
        else
          STATUS = nf90_get_att(ncid, nf90_global, "in_b_inuc", tb_b_inuc)
          if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
        endif
        ! coefficient N_inuc_R
        STATUS = nf90_inquire_attribute(ncid, nf90_global, "in_n_inucr")
        if (STATUS .ne. nf90_noerr) then
          write(6,*) 'NOTICE modtestbed: in_b_inuc not in scm_in.nc'
          ltb_setinp  = .false.      
        else
          STATUS = nf90_get_att(ncid, nf90_global, "in_n_inucr", tb_n_inuc_r)
          if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
        endif
        ! coefficient b_inuc_R
        STATUS = nf90_inquire_attribute(ncid, nf90_global, "in_b_inucr")
        if (STATUS .ne. nf90_noerr) then
          write(6,*) 'NOTICE modtestbed: in_b_inuc not in scm_in.nc'
          ltb_setinp  = .false.      
        else
          STATUS = nf90_get_att(ncid, nf90_global, "in_b_inucr", tb_b_inuc_r)
          if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
        endif  
        ! add message about that 
        if (ltb_setinp) then
          write(6,*) 'INFO modtestbed: inp properties loaded from in scm_in.nc'
          ! write(6,*) 'INFO modtestbed: '
        endif

        ! loading scalars
          start_sv  = (/ 1       , 1       , 1       /)
          count_sv  = (/ nsv     , nknudge , ntnudge /) 
          !
          start_svs = (/ 1       , 1       /)
          count_svs = (/ nsv     , ntnudge /)
          ! ltb_sv = .true.  <- done before
          ! other options:
          !2)
          ! count_sv  = (/ nsv     , ntnudge , nknudge /) 
          !3)
          ! count_sv  = (/ nknudge , ntnudge ,  nsv     /)
          !     
        !if (ltb_loadsv) then ! if (nsv>0 .and. ltb_loadsv) then
        !
        !  ! surface timeseries
        !  write(6,*) 'Now loading sv variables' 
        !  write(6,*) count_sv
        !  !  sv profiles for nsv species
        !  STATUS = NF90_INQ_VARID(NCID, 'sv', VARID)
        !  if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
        !  STATUS = NF90_GET_VAR (NCID, VARID, dumsv, start=start_sv, count=count_sv)
        !  if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
        !
        ! endif
        
        !if (ltb_svadv) then  
        !  write(6,*) 'Now loading sv variables' 
        !  write(6,*) count_sv
        !  !  sv profiles for nsv species
        !  STATUS = NF90_INQ_VARID(NCID, 'svadv', VARID)
        !  if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
        !  STATUS = NF90_GET_VAR (NCID, VARID, dumsv_adv, start=start_sv, count=count_sv)
        !  if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
        !
        !endif
        
        !if (ltb_svsurf) then ! if (nsv>0 .and. ltb_svsurf) then
        !  !
        !  write(6,*) 'Now loading svs variables' 
        !  write(6,*) count_svs
        !  STATUS = NF90_INQ_VARID(NCID, 'sv_flx', VARID)
        !  if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
        !  STATUS = NF90_GET_VAR (NCID, VARID, dumsvs, start=start_svs, count=count_svs)
        !  if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
        !         
        !endif
         
        !--- close nc file ---
        STATUS = NF90_CLOSE(NCID)  
        if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
        


        !--- interpolate towards LES levels, reverse height-order if needed, switch dimensions ---
        do i=1,ntnudge


          !ik = nknudge
          do k=1,k1
            
            !  look for level closest to zf
            dummax = maxval(zf)
            ik0 = k1
            do ik=2,nknudge-1
              dumdz = abs(zf(k)-dumheight(ik,i))
              if ( dumdz .lt. dummax ) then
                dummax = dumdz
                ik0 = ik 
              endif 
            enddo
            !  determine in which way the height array increases 
            ikdir = int( sign( 1., dumheight(2,i) - dumheight(1,i) ) )
            if ( zf(k).gt.dumheight(ik0,i) ) then
              ik1 = ik0 + ikdir 
            else
              ik1 = ik0 - ikdir 
            endif
          
            fac = ( zf(k)-dumheight(ik0,i) ) / ( dumheight(ik1,i)-dumheight(ik0,i) )
        
            tb_thl    (i,k)   = dumthl    (ik0,i) + fac * ( dumthl    (ik1,i)    - dumthl    (ik0,i)    )
            tb_qt     (i,k)   = dumqt     (ik0,i) + fac * ( dumqt     (ik1,i)    - dumqt     (ik0,i)    )
            tb_nccn   (i,k)   = dumnccn   (ik0,i) + fac * ( dumnccn   (ik1,i)    - dumnccn   (ik0,i)    )
            tb_ql     (i,k)   = dumql     (ik0,i) + fac * ( dumql     (ik1,i)    - dumql     (ik0,i)    )
            tb_qi     (i,k)   = dumqi     (ik0,i) + fac * ( dumqi     (ik1,i)    - dumqi     (ik0,i)    )
            tb_dcl    (i,k)   = dumdcl    (ik0,i) + fac * ( dumdcl    (ik1,i)    - dumdcl    (ik0,i)    )
            tb_dci    (i,k)   = dumdci    (ik0,i) + fac * ( dumdci    (ik1,i)    - dumdci    (ik0,i)    )
            tb_u      (i,k)   = dumu      (ik0,i) + fac * ( dumu      (ik1,i)    - dumu      (ik0,i)    )
            tb_v      (i,k)   = dumv      (ik0,i) + fac * ( dumv      (ik1,i)    - dumv      (ik0,i)    )
            tb_w      (i,k)   = dumw      (ik0,i) + fac * ( dumw      (ik1,i)    - dumw      (ik0,i)    )
            tb_ug     (i,k)   = dumug     (ik0,i) + fac * ( dumug     (ik1,i)    - dumug     (ik0,i)    )
            tb_vg     (i,k)   = dumvg     (ik0,i) + fac * ( dumvg     (ik1,i)    - dumvg     (ik0,i)    )
            tb_uadv   (i,k)   = dumuadv   (ik0,i) + fac * ( dumuadv   (ik1,i)    - dumuadv   (ik0,i)    )
            tb_vadv   (i,k)   = dumvadv   (ik0,i) + fac * ( dumvadv   (ik1,i)    - dumvadv   (ik0,i)    )
            tb_qtadv  (i,k)   = dumqtadv  (ik0,i) + fac * ( dumqtadv  (ik1,i)    - dumqtadv  (ik0,i)    )
            tb_qladv  (i,k)   = dumladv   (ik0,i) + fac * ( dumladv   (ik1,i)    - dumladv   (ik0,i)    )
            tb_qiadv  (i,k)   = dumiadv   (ik0,i) + fac * ( dumiadv   (ik1,i)    - dumiadv   (ik0,i)    )
            tb_thladv (i,k)   = dumthladv (ik0,i) + fac * ( dumthladv (ik1,i)    - dumthladv (ik0,i)    )

            !if (i.eq.1) write(6,*)  k, zf(k), " : ", ik, dumheight(ik,i), ik-1, dumheight(ik-1,i)
            ! and for scalars
            
            ! if (ltb_sv .and. i.eq.1) then  ! #442fredrik commented out because it cannot be D_BCAST
            !if (ltb_sv) then
            !  !tb_sv(i,k,:) = dumsv(:,ik,i) + fac * ( dumsv(:,ik-1,i) - dumsv(:,ik,i) )
            !  tb_sv(i,k,:) = dumsv(:,ik0,i) + fac * ( dumsv(:,ik1,i) - dumsv(:,ik0,i) )
            !  !tb_svadv(i,k,:) = dumsv_adv(:,ik,i) + fac * ( dumsv_adv(:,ik-1,i) - dumsv_adv(:,ik,i) )
            !  tb_svadv(i,k,:) = dumsv_adv(:,ik0,i) + fac * ( dumsv_adv(:,ik1,i) - dumsv_adv(:,ik0,i) )
            !endif           
          enddo
          

          !--- soil & surface properties ---
          !
          !    note: The following surface properties from scm_in.nc are used in the simulation: 
          !            tb_lat, tb_lon, tb_thls, tb_ps, tb_alb, tb_wts, tb_wqs, tb_z0m, tb_z0h,
          !            tb_seaicefrct, tb_thls_ocean, tb_thls_seaice
          !          See for example modstartup.f90 and modtimedep.f90 (timedepsurf) for ltestbed=T.
          !
          !          The following properties are read from scm_in.nc but still need to be passed on in the code: 
          !            tb_tsoilav, tb_phiwav
          !          For these the values provided through namoptions are still used
          !   
          if (ltb_soildata) then
           tb_qts(i) = tb_qt(i,1)    !qts seems not really used anymore (see subr. timedepsurf in modtimedep.f90)

           tb_tsoilav(i,:) = dumtsoilav(:,i)

           tb_phiwav (i,:) = dumphiwav (:,i)      
           !  calculate soil moisture matching LES soil model, by scaling input soil wetness index with native field capacity & wilting point
           !tb_phiwav (i,:) = phiwp + dumswi(:,i) * (phifc - phiwp )
           write(6,*) 'inittestbed: surface properties loaded but currently not used' 
          endif

          !--- profiles for full radiation scheme (iradiation=1: see modradfull, d4stream_tb_setup) ---
          !
          do k = 1, nknudge
            tbrad_p(i,k)  = dumpf(k,i)
            tbrad_t(i,k)  = dumt(k,i)
            tbrad_qv(i,k) = dumqv(k,i)
            tbrad_ql(i,k) = dumql(k,i) + dumqi(k,i) 
            tbrad_qi(i,k) = dumqi(k,i) 
            tbrad_o3(i,k) = dumo3(k,i)
          end do

        enddo
        
        ! ------- late start, ----------------------------------
        ! purpose: interpolate first value and shift later ones
        if(ltb_latestart) then
          ! find the time index 
          i = 1    ! start loop
          do while (tb_time(i).le.tbtimeoffset) 
            i = i+1
          end do 
          i0tb = i-1
          write(6,*) "message testbed: original scm_in profiles used after index ",i0tb
         ! new first index 
         tfac = (tbtimeoffset - tb_time(i0tb))/(tb_time(i0tb+1)-tb_time(i0tb))
         i=1 
         ! shift time 
          tb_time (i)  =  0.0 ! tbtimeoffset
         ! following lines generated by a script
         ! surface 1d variables 
          tb_lat   (i) = tb_lat  (i0tb) +tfac*(tb_lat  (i0tb+1)-tb_lat  (i0tb)) 
          tb_lon   (i) = tb_lon  (i0tb) +tfac*(tb_lon  (i0tb+1)-tb_lon  (i0tb)) 
          tb_ps    (i) = tb_ps   (i0tb) +tfac*(tb_ps   (i0tb+1)-tb_ps   (i0tb)) 
          tb_thls  (i) = tb_thls (i0tb) +tfac*(tb_thls (i0tb+1)-tb_thls (i0tb)) 
          tb_thls_ocean (i) = tb_thls_ocean(i0tb) +tfac*(tb_thls_ocean(i0tb+1)-tb_thls_ocean(i0tb)) 
          tb_thls_seaice(i) = tb_thls_seaice(i0tb) +tfac*(tb_thls_seaice(i0tb+1)-tb_thls_seaice(i0tb)) 
          tb_wts   (i) = tb_wts  (i0tb) +tfac*(tb_wts  (i0tb+1)-tb_wts  (i0tb)) 
          tb_wqs   (i) = tb_wqs  (i0tb) +tfac*(tb_wqs  (i0tb+1)-tb_wqs  (i0tb)) 
          tb_seaicefrct(i)  = tb_seaicefrct(i0tb) +tfac*(tb_seaicefrct(i0tb+1)-tb_seaicefrct(i0tb)) 
          tb_alb   (i) = tb_alb  (i0tb) +tfac*(tb_alb  (i0tb+1)-tb_alb  (i0tb)) 
          tb_z0h   (i) = tb_z0h  (i0tb) +tfac*(tb_z0h  (i0tb+1)-tb_z0h  (i0tb)) 
          tb_z0m   (i) = tb_z0m  (i0tb) +tfac*(tb_z0m  (i0tb+1)-tb_z0m  (i0tb)) 
          tb_Qnet  (i) = tb_Qnet (i0tb) +tfac*(tb_Qnet (i0tb+1)-tb_Qnet (i0tb)) 
          tb_qts   (i) = tb_qts  (i0tb) +tfac*(tb_qts  (i0tb+1)-tb_qts  (i0tb)) 
         ! profiles variables 
         do k=1,k1 
          tb_thl   (i,k) = tb_thl  (i0tb,k) +tfac*(tb_thl  (i0tb+1,k)-tb_thl  (i0tb,k)) 
          tb_qt    (i,k) = tb_qt   (i0tb,k) +tfac*(tb_qt   (i0tb+1,k)-tb_qt   (i0tb,k)) 
          tb_nccn  (i,k) = tb_nccn (i0tb,k) +tfac*(tb_nccn (i0tb+1,k)-tb_nccn (i0tb,k)) 
          tb_ql    (i,k) = tb_ql   (i0tb,k) +tfac*(tb_ql   (i0tb+1,k)-tb_ql   (i0tb,k)) 
          tb_qi    (i,k) = tb_qi   (i0tb,k) +tfac*(tb_qi   (i0tb+1,k)-tb_qi   (i0tb,k)) 
          tb_dcl   (i,k) = tb_dcl  (i0tb,k) +tfac*(tb_dcl  (i0tb+1,k)-tb_dcl  (i0tb,k))
          tb_dci   (i,k) = tb_dci  (i0tb,k) +tfac*(tb_dci  (i0tb+1,k)-tb_dci  (i0tb,k))
          tb_u     (i,k) = tb_u    (i0tb,k) +tfac*(tb_u    (i0tb+1,k)-tb_u    (i0tb,k)) 
          tb_v     (i,k) = tb_v    (i0tb,k) +tfac*(tb_v    (i0tb+1,k)-tb_v    (i0tb,k)) 
          tb_w     (i,k) = tb_w    (i0tb,k) +tfac*(tb_w    (i0tb+1,k)-tb_w    (i0tb,k)) 
          tb_ug    (i,k) = tb_ug   (i0tb,k) +tfac*(tb_ug   (i0tb+1,k)-tb_ug   (i0tb,k)) 
          tb_vg    (i,k) = tb_vg   (i0tb,k) +tfac*(tb_vg   (i0tb+1,k)-tb_vg   (i0tb,k)) 
          tb_uadv  (i,k) = tb_uadv (i0tb,k) +tfac*(tb_uadv (i0tb+1,k)-tb_uadv (i0tb,k)) 
          tb_vadv  (i,k) = tb_vadv (i0tb,k) +tfac*(tb_vadv (i0tb+1,k)-tb_vadv (i0tb,k)) 
          tb_qtadv (i,k) = tb_qtadv(i0tb,k) +tfac*(tb_qtadv(i0tb+1,k)-tb_qtadv(i0tb,k)) 
          tb_qladv (i,k) = tb_qladv(i0tb,k) +tfac*(tb_qladv(i0tb+1,k)-tb_qladv(i0tb,k))
          tb_qiadv (i,k) = tb_qiadv(i0tb,k) +tfac*(tb_qiadv(i0tb+1,k)-tb_qiadv(i0tb,k))
          tb_thladv (i,k)  = tb_thladv(i0tb,k) +tfac*(tb_thladv(i0tb+1,k)-tb_thladv(i0tb,k)) 
         end do 
         ! soil profiles variables 
         do k=1,ksoilmax 
          tb_tsoilav(i,k) = tb_tsoilav(i0tb,k) +tfac*(tb_tsoilav(i0tb+1,k)-tb_tsoilav(i0tb,k)) 
          tb_phiwav (i,k) = tb_phiwav(i0tb,k) +tfac*(tb_phiwav(i0tb+1,k)-tb_phiwav(i0tb,k)) 
         end do 
         ! radiation profiles variables 
         do k = 1, nknudge 
          tbrad_p     (i,k) = tbrad_p    (i0tb,k) +tfac*(tbrad_p    (i0tb+1,k)-tbrad_p    (i0tb,k)) 
          tbrad_t     (i,k) = tbrad_t    (i0tb,k) +tfac*(tbrad_t    (i0tb+1,k)-tbrad_t    (i0tb,k)) 
          tbrad_qv    (i,k) = tbrad_qv   (i0tb,k) +tfac*(tbrad_qv   (i0tb+1,k)-tbrad_qv   (i0tb,k)) 
          tbrad_ql    (i,k) = tbrad_ql   (i0tb,k) +tfac*(tbrad_ql   (i0tb+1,k)-tbrad_ql   (i0tb,k)) 
          tbrad_qi    (i,k) = tbrad_qi   (i0tb,k) +tfac*(tbrad_qi   (i0tb+1,k)-tbrad_qi   (i0tb,k)) 
          tbrad_o3    (i,k) = tbrad_o3   (i0tb,k) +tfac*(tbrad_o3   (i0tb+1,k)-tbrad_o3   (i0tb,k))
         end do 
         ! shift later indices 
         do i=2,(ntnudge-i0tb+1) 
         ! time
          tb_time (i)  = tb_time (i0tb+i-1) - tbtimeoffset
         ! surface 1d variables
          tb_lat   (i) = tb_lat  (i0tb+i-1) 
          tb_lon   (i) = tb_lon  (i0tb+i-1) 
          tb_ps    (i) = tb_ps   (i0tb+i-1) 
          tb_thls  (i) = tb_thls (i0tb+i-1) 
          tb_thls_ocean (i)= tb_thls_ocean(i0tb+i-1) 
          tb_thls_seaice(i)= tb_thls_seaice(i0tb+i-1) 
          tb_wts   (i) = tb_wts  (i0tb+i-1) 
          tb_wqs   (i) = tb_wqs  (i0tb+i-1) 
          tb_seaicefrct(i) = tb_seaicefrct(i0tb+i-1) 
          tb_alb   (i) = tb_alb  (i0tb+i-1) 
          tb_z0h   (i) = tb_z0h  (i0tb+i-1) 
          tb_z0m   (i) = tb_z0m  (i0tb+i-1) 
          tb_Qnet  (i) = tb_Qnet (i0tb+i-1) 
          tb_qts   (i) = tb_qts  (i0tb+i-1) 
         ! profiles variables 
         do k=1,k1 
          tb_thl   (i,k) = tb_thl  (i0tb+i-1,k) 
          tb_qt    (i,k) = tb_qt   (i0tb+i-1,k) 
          tb_nccn  (i,k) = tb_nccn (i0tb+i-1,k) 
          tb_ql    (i,k) = tb_ql   (i0tb+i-1,k) 
          tb_qi    (i,k) = tb_qi   (i0tb+i-1,k) 
          tb_dcl   (i,k) = tb_dcl  (i0tb+i-1,k)
          tb_dci   (i,k) = tb_dci  (i0tb+i-1,k)
          tb_u     (i,k) = tb_u    (i0tb+i-1,k) 
          tb_v     (i,k) = tb_v    (i0tb+i-1,k) 
          tb_w     (i,k) = tb_w    (i0tb+i-1,k) 
          tb_ug    (i,k) = tb_ug   (i0tb+i-1,k) 
          tb_vg    (i,k) = tb_vg   (i0tb+i-1,k) 
          tb_uadv  (i,k) = tb_uadv (i0tb+i-1,k) 
          tb_vadv  (i,k) = tb_vadv (i0tb+i-1,k) 
          tb_qtadv (i,k) = tb_qtadv(i0tb+i-1,k) 
          tb_qladv (i,k) = tb_qladv(i0tb+i-1,k)
          tb_qiadv (i,k) = tb_qiadv(i0tb+i-1,k)
          tb_thladv(i,k) = tb_thladv(i0tb+i-1,k) 
         end do 
         ! soil profiles variables 
         do k=1,ksoilmax 
          tb_tsoilav (i,k)= tb_tsoilav(i0tb+i-1,k) 
          tb_phiwav (i,k) = tb_phiwav(i0tb+i-1,k) 
         end do 
          ! radiation profiles variables 
         do k = 1, nknudge 
          tbrad_p     (i,k) = tbrad_p    (i0tb+i-1,k) 
          tbrad_t     (i,k) = tbrad_t    (i0tb+i-1,k) 
          tbrad_qv    (i,k) = tbrad_qv   (i0tb+i-1,k) 
          tbrad_ql    (i,k) = tbrad_ql   (i0tb+i-1,k) 
          tbrad_qi    (i,k) = tbrad_qi   (i0tb+i-1,k) 
          tbrad_o3    (i,k) = tbrad_o3   (i0tb+i-1,k)
         end do 
        end do 
        ! #dbg output
         write(6,*) "message testbed: scm_in init profile interpolated at tbtimeoffset=", tbtimeoffset
         write(6,*) "height    thl      qt         u      v  "
         do k=(k1-1),1,-1
          write (6,'(f7.1,f8.1,e12.4,2f7.1)') &   ! #dbg
                zf     (k), &
                tb_thl (1,k), &
                tb_qt  (1,k), &
                tb_u   (1,k), &
                tb_v   (1,k)
         enddo
        endif

        
        !--- clean-up ---
        deallocate(dumomega,dumqv,dumql,dumqi,dumt,dumpf,dumo3)
        deallocate(dumheight,dumqt,dumnccn,dumthl,dumu,dumv,dumw,dumug,dumvg)
        deallocate(dumuadv,dumvadv,dumqtadv,dumthladv,dumqadv,dumladv,dumiadv,dumtadv,dumswnet,dumlwnet)
        deallocate(dumheights,dumtsoilav,dumphiwav,dumswi)
        
        !if (ltb_sv) then ! if (ltb_loadsv) then
        !  deallocate(dumsv)
        !  deallocate(dumsv_adv)
        !  deallocate(dumsvs)
        !endif

        !--- do some output to screen ---
!        do i=1,2
!        !do i=1,ntnudge
!
!        write(6,'(a20,f10.2,a15,3f10.2)') 'modtestbed: scm_in time:',tb_time(i),' sfc pressure:',tb_ps(i),tb_thls(i),tb_qts(i)
!
!        write(6,*) ' zf       tnudge    tb_u    tb_v    tb_w    tb_thl    tb_qt    tb_ug    tb_vg'
!        do k=kmax,1,-1
!          write (6,'(f7.1,8e12.4)') &
!                zf          (k), &
!                tnudge      (i,k), &
!                tb_u        (i,k), &
!                tb_v        (i,k), &
!                tb_w        (i,k), &
!                tb_thl      (i,k), &
!                tb_qt       (i,k), &
!                tb_ug       (i,k), &
!                tb_vg       (i,k)
!        end do
!
!        end do

    end if

    call D_MPI_BCAST(ntnudge    , 1  ,0,comm3d,mpierr)

    call D_MPI_BCAST(tb_time    ,ntnudge         ,0,comm3d,mpierr)
    call D_MPI_BCAST(tb_lat     ,ntnudge         ,0,comm3d,mpierr)
    call D_MPI_BCAST(tb_lon     ,ntnudge         ,0,comm3d,mpierr)
    call D_MPI_BCAST(tb_ps      ,ntnudge         ,0,comm3d,mpierr)
    call D_MPI_BCAST(tb_qts     ,ntnudge         ,0,comm3d,mpierr)
    call D_MPI_BCAST(tb_thls    ,ntnudge         ,0,comm3d,mpierr)
    call D_MPI_BCAST(tb_wts     ,ntnudge         ,0,comm3d,mpierr)
    call D_MPI_BCAST(tb_wqs     ,ntnudge         ,0,comm3d,mpierr)
    call D_MPI_BCAST(tb_z0h     ,ntnudge         ,0,comm3d,mpierr)! not needed, see: inittimedep
    call D_MPI_BCAST(tb_z0m     ,ntnudge         ,0,comm3d,mpierr)! not needed, see: inittimedep
    call D_MPI_BCAST(tb_alb     ,ntnudge         ,0,comm3d,mpierr)! not needed, see: inittimedep
    call D_MPI_BCAST(tb_Qnet    ,ntnudge         ,0,comm3d,mpierr)
    call D_MPI_BCAST(tb_seaicefrct,ntnudge       ,0,comm3d,mpierr)! not needed, see: inittimedep
    call D_MPI_BCAST(tb_thls_ocean  ,ntnudge         ,0,comm3d,mpierr)! not needed, see: inittimedep
    call D_MPI_BCAST(tb_thls_seaice ,ntnudge         ,0,comm3d,mpierr)! not needed, see: inittimedep


    call D_MPI_BCAST(tnudge     ,ntnudge*k1      ,0,comm3d,mpierr)
    call D_MPI_BCAST(tb_u       ,ntnudge*k1      ,0,comm3d,mpierr)
    call D_MPI_BCAST(tb_v       ,ntnudge*k1      ,0,comm3d,mpierr)
    call D_MPI_BCAST(tb_w       ,ntnudge*k1      ,0,comm3d,mpierr)
    call D_MPI_BCAST(tb_thl     ,ntnudge*k1      ,0,comm3d,mpierr)
    call D_MPI_BCAST(tb_qt      ,ntnudge*k1      ,0,comm3d,mpierr)
    call D_MPI_BCAST(tb_ql      ,ntnudge*k1      ,0,comm3d,mpierr)
    call D_MPI_BCAST(tb_qi      ,ntnudge*k1      ,0,comm3d,mpierr)
    call D_MPI_BCAST(tb_dcl     ,ntnudge*k1      ,0,comm3d,mpierr)
    call D_MPI_BCAST(tb_dci     ,ntnudge*k1      ,0,comm3d,mpierr)    
    call D_MPI_BCAST(tb_nccn    ,ntnudge*k1      ,0,comm3d,mpierr)    
    call D_MPI_BCAST(tb_ug      ,ntnudge*k1      ,0,comm3d,mpierr)
    call D_MPI_BCAST(tb_vg      ,ntnudge*k1      ,0,comm3d,mpierr)
    call D_MPI_BCAST(tb_uadv    ,ntnudge*k1      ,0,comm3d,mpierr)
    call D_MPI_BCAST(tb_vadv    ,ntnudge*k1      ,0,comm3d,mpierr)
    call D_MPI_BCAST(tb_qtadv   ,ntnudge*k1      ,0,comm3d,mpierr)
    call D_MPI_BCAST(tb_qladv   ,ntnudge*k1      ,0,comm3d,mpierr)
    call D_MPI_BCAST(tb_qiadv   ,ntnudge*k1      ,0,comm3d,mpierr)
    call D_MPI_BCAST(tb_thladv  ,ntnudge*k1      ,0,comm3d,mpierr)
    call D_MPI_BCAST(tb_uadv    ,ntnudge*k1      ,0,comm3d,mpierr)
    call D_MPI_BCAST(tb_vadv    ,ntnudge*k1      ,0,comm3d,mpierr)
    call D_MPI_BCAST(tb_dqtdxls ,ntnudge*k1      ,0,comm3d,mpierr)
    call D_MPI_BCAST(tb_dqtdyls ,ntnudge*k1      ,0,comm3d,mpierr)

    call D_MPI_BCAST(tb_tsoilav ,ntnudge*ksoilmax      ,0,comm3d,mpierr)
    call D_MPI_BCAST(tb_phiwav  ,ntnudge*ksoilmax      ,0,comm3d,mpierr)

    call D_MPI_BCAST(tb_land_use, mpatch*mpatch         ,0, comm3d, mpierr)

    call D_MPI_BCAST(tbrad_p      ,ntnudge*nknudge      ,0,comm3d,mpierr)
    call D_MPI_BCAST(tbrad_qv     ,ntnudge*nknudge      ,0,comm3d,mpierr)
    call D_MPI_BCAST(tbrad_ql     ,ntnudge*nknudge      ,0,comm3d,mpierr)
    call D_MPI_BCAST(tbrad_qi     ,ntnudge*nknudge      ,0,comm3d,mpierr)
    call D_MPI_BCAST(tbrad_t      ,ntnudge*nknudge      ,0,comm3d,mpierr)
    call D_MPI_BCAST(tbrad_o3     ,ntnudge*nknudge      ,0,comm3d,mpierr)
    
    !if (nsv>0) then  ! #442fredrik : block commented out because it cannot be D_BCAST
    !  call D_MPI_BCAST(tb_sv        ,ntnudge*k1*nsv       ,0,comm3d,mpierr)! not needed, see: 
    !  !call D_MPI_BCAST(tb_svadv     ,ntnudge*k1*nsv       ,0,comm3d,mpierr)
    !endif 
    ! casting load flags
    call D_MPI_BCAST(ltb_setccn      , 1  ,0,comm3d,mpierr)
    call D_MPI_BCAST(ltb_setinp      , 1  ,0,comm3d,mpierr)
    call D_MPI_BCAST(tb_n_inuc       , 1      ,0,comm3d,mpierr)
    call D_MPI_BCAST(tb_a_inuc       , 1      ,0,comm3d,mpierr)
    call D_MPI_BCAST(tb_b_inuc       , 1      ,0,comm3d,mpierr)
    call D_MPI_BCAST(tb_n_inuc_r     , 1      ,0,comm3d,mpierr)
    call D_MPI_BCAST(tb_b_inuc_r     , 1      ,0,comm3d,mpierr)
    call D_MPI_BCAST(tb_xday         , 1      ,0,comm3d,mpierr)
    
    ! update of the time
    if (ltb_usedate.and.(tb_xday.gt.0.0)) then
      xday = tb_xday   ! modglobal time xday updated 
    endif
    
    ! outputs of tbrad arrays:
!     write(6,*) "inittestbed: scm inputs tbrad on myid", myid
!     write(6,*) "  tbrad_p:"
!     write(6,*) tbrad_p
!     write(6,*) "  tbrad_t:"
!     write(6,*) tbrad_t
!     write(6,*) "  tbrad_qv:"
!     write(6,*) tbrad_qv
!     write(6,*) "  tbrad_ql:"
!     write(6,*) tbrad_ql
!     write(6,*) "  tbrad_qi:"
!     write(6,*) tbrad_qi
!     write(6,*) "  tbrad_o3:"
!     write(6,*) tbrad_o3
!     write(6,*) "-----"

  end subroutine inittestbed

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine testbednudge
    use modglobal, only : timee,rtimee,i1,j1,kmax,rdt, zf ! #sb3
    use modfields, only : up,vp,wp,thlp, qtp,u0av,v0av,qt0av,thl0av
    implicit none

    integer k,t,kinv
    real :: dtm,dtp,currtnudge, qttnudge,qtthres, zfac, zfacu, zmin, zmid, dthldz_max1, dthldz(kmax)
    
    if (.not.(ltestbed .and. ltb_nudge)) return

    if (timee==0) return

    t=1
    do while(rtimee>tb_time(t))
      t=t+1
    end do
    if (rtimee>tb_time(1)) then
      t=t-1
    end if

    dtm = ( rtimee-tb_time(t) ) / ( tb_time(t+1)-tb_time(t) )
    dtp = ( tb_time(t+1)-rtimee)/ ( tb_time(t+1)-tb_time(t) )

    !--- nudging above BL top ---
    ! as defined by the level with the strongest thl gradient (thermal inversion), which can be time dependent
    ! to use this set tb_zminnudge=-1 in namoptions
    ! the value of tb_zmidnudge then sets the depth of the nudging transition layer between zmin (the inversion height) and zmid
    !
   if (ltb_spinnudge.and.(rtimee<tb_t_spinnudge)) then
      !-- use specified nudging height range --
      zmin =   tb_zmin_spinnudge ! tb_zminnudge
      zmid =   tb_zmid_spinnudge ! tb_zmidnudge
    
    qtthres = 1.e-6
    do k=1,kmax
    
      zfac = 0. 
      if ( zf(k).gt. zmin) then
        ! set factor in nudging intensity: linear increase from 0 to 1 in heightrange  tb_zmidnudge > z > tb_zminnudge
        zfac = zf(k) - zmin
        if (zmid.gt.zmin) then
          zfac = max( 0., min( 1., zfac / (zmid - zmin) ) )
        else
          zfac = 0.5 + sign( 0.5, zfac )
        endif
        !write(0,*) 'modtestbed: ', k, zf(k), zmin, zmid, zfac
      endif

      zfacu = zfac
      !zfacu = 1.
 
      ! time step protection - stronger nudging here
      currtnudge = max(rdt,tb_tauhigh) ! max(rdt,tnudge(t,k)*dtp+tnudge(t+1,k)*dtm)

      if (ltb_u)   upnu  (k) =  - zfacu * &
          ( u0av(k)   - (tb_u(t,k)  *dtp + tb_u(t+1,k)  *dtm) ) / currtnudge

      if (ltb_v)   vpnu  (k) =  - zfacu * &
          ( v0av(k)   - (tb_v(t,k)  *dtp + tb_v(t+1,k)  *dtm) ) / currtnudge

      if (ltb_w)   wpnu  (k) =  - zfac * &
          (           - (tb_w(t,k)  *dtp + tb_w(t+1,k)  *dtm) ) / currtnudge

      if (ltb_thl) thlpnu(k) =  - zfac * &
          ( thl0av(k) - (tb_thl(t,k)*dtp + tb_thl(t+1,k)*dtm) ) / currtnudge

      if (ltb_qt)  then
        if (qt0av(k)< qtthres) then
          qttnudge = rdt
        else
          qttnudge = currtnudge
        endif
        qtpnu (k) =    - zfac * &
          ( qt0av(k)  - (tb_qt(t,k) *dtp + tb_qt(t+1,k) *dtm) ) / qttnudge
      endif
      
      ! do updates 
      up  (2:i1,2:j1,k) = up  (2:i1,2:j1,k) + upnu  (k)
      vp  (2:i1,2:j1,k) = vp  (2:i1,2:j1,k) + vpnu  (k)
      wp  (2:i1,2:j1,k) = wp  (2:i1,2:j1,k) + wpnu  (k)
      thlp(2:i1,2:j1,k) = thlp(2:i1,2:j1,k) + thlpnu(k)
      qtp (2:i1,2:j1,k) = qtp (2:i1,2:j1,k) + qtpnu (k)

    end do    
   else ! i.e. ltb_spinnudge==.false. .or.(rtimee> tb_t_spinnudge)
    if (tb_zminnudge .le. -0.9 ) then
      dthldz_max1 = -100000.
      kinv = -1
      k = 2  ! initializing
      do while ( zf(k).le. tb_minzinv ) ! to avoid near ground inversion
        k = k+1
      end do
      do while ((k.lt.(kmax-1)).and.( zf(k).lt. tb_maxzinv) ) !do k=2,kmax-1
        dthldz(k) = ( thl0av(k+1) - thl0av(k-1) ) / ( zf(k+1) - zf(k-1) )
        !write(0,'(a,i4,a,f15.7,f8.2)') '  modtestbed:   k=',k,' dthl/dz=',dthldz(k),thl0av(k)
        if ( dthldz(k).gt.dthldz_max1 ) then
          kinv = k
          dthldz_max1 = dthldz(k)
        endif 
        k =  k+1
      end do
      if (kinv.gt.0) then
        ! write(6,'(a,i4,a,f8.2,a,f15.7)') '  modtestbed: inversion: kinv=',kinv,'  zinv=',zf(kinv),' dthl/dz=',dthldz(kinv)
        zmin = zf(kinv)
      else
        write(0,'(a)') '  modtestbed: inversion: not found'
        zmin = 0.
      endif
      zmid = zmin + tb_zmidnudge  
    else
      !-- use specified nudging height range --
      zmin = tb_zminnudge
      zmid = tb_zmidnudge
    endif
    
    qtthres = 1e-6
    do k=1,kmax
    
      zfac = 0. 
      if ( zf(k).gt. zmin) then
        ! set factor in nudging intensity: linear increase from 0 to 1 in heightrange  tb_zmidnudge > z > tb_zminnudge
        zfac = zf(k) - zmin
        if (zmid.gt.zmin) then
          zfac = max( 0., min( 1., zfac / (zmid - zmin) ) )
        else
          zfac = 0.5 + sign( 0.5, zfac )
        endif
        !write(0,*) 'modtestbed: ', k, zf(k), zmin, zmid, zfac
      endif

      zfacu = zfac
      !zfacu = 1.
 
      ! time step protection
      currtnudge = max(rdt,tnudge(t,k)*dtp+tnudge(t+1,k)*dtm)

      if (ltb_u)   upnu  (k) =  - zfacu * &
          ( u0av(k)   - (tb_u(t,k)  *dtp + tb_u(t+1,k)  *dtm) ) / currtnudge

      if (ltb_v)   vpnu  (k) =  - zfacu * &
          ( v0av(k)   - (tb_v(t,k)  *dtp + tb_v(t+1,k)  *dtm) ) / currtnudge

      if (ltb_w)   wpnu  (k) =  - zfac * &
          (           - (tb_w(t,k)  *dtp + tb_w(t+1,k)  *dtm) ) / currtnudge

      if (ltb_thl) thlpnu(k) =  - zfac * &
          ( thl0av(k) - (tb_thl(t,k)*dtp + tb_thl(t+1,k)*dtm) ) / currtnudge

      if (ltb_qt)  then
        if (qt0av(k)< qtthres) then
          qttnudge = rdt
        else
          qttnudge = currtnudge
        endif
        qtpnu (k)            =    - zfac * &
          ( qt0av(k)  - (tb_qt(t,k) *dtp + tb_qt(t+1,k) *dtm) ) / qttnudge
      endif
      
      ! do updates 
      up  (2:i1,2:j1,k) = up  (2:i1,2:j1,k) + upnu  (k)
      vp  (2:i1,2:j1,k) = vp  (2:i1,2:j1,k) + vpnu  (k)
      wp  (2:i1,2:j1,k) = wp  (2:i1,2:j1,k) + wpnu  (k)
      thlp(2:i1,2:j1,k) = thlp(2:i1,2:j1,k) + thlpnu(k)
      qtp (2:i1,2:j1,k) = qtp (2:i1,2:j1,k) + qtpnu (k)

    end do
   end if
   
    
    ! if(ltb_sv.and. ltb_nudge) then 
    !    write(6,*) 'testbednudge:','nudging for scalars not yet included'
    ! endif

    !write(6,*) 'testbednudge:', rtimee, t, tb_time(t), tb_time(t+1), currtnudge, dtm, dtp, qt0av (1),tb_qt (t,1),tb_qt (t+1,1)
    !write(6,*) 'testbednudge:', rtimee, t, tb_time(t), tb_time(t+1), currtnudge, dtm, dtp, qt0av (kmax),tb_qt (t,kmax),tb_qt (t+1,kmax)
    
    !write(6,*) 'testbednudge:', rtimee, t, tb_time(t), tb_time(t+1), currtnudge, dtm, dtp, thl0av(1),tb_thl(t,1),tb_thl(t+1,1)
    !write(6,*) 'testbednudge:', rtimee, t, tb_time(t), tb_time(t+1), currtnudge, dtm, dtp, thl0av(kmax),tb_thl(t,kmax),tb_thl(t+1,kmax)

  end subroutine testbednudge
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine testbed_getinttime(t, dtm, dtp)
     use modglobal, only : rtimee
!     use modfields, only : up,vp,wp,thlp, qtp,u0av,v0av,qt0av,thl0av
!     use modmpi,    only : myid
    implicit none
    integer, intent(out) :: t
    real, intent(out)    :: dtm, dtp


    t=1
    do while(rtimee>tb_time(t))
      t=t+1
    end do
    if (rtimee>tb_time(1)) then
      t=t-1
    end if

    dtm = ( rtimee-tb_time(t) ) / ( tb_time(t+1)-tb_time(t) )
    dtp = ( tb_time(t+1)-rtimee)/ ( tb_time(t+1)-tb_time(t) )


  end subroutine testbed_getinttime

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine exittestbed
   use modglobal, only: nsv
   if (allocated(tb_time)) then
    deallocate(  tnudge,tb_u,tb_v,tb_w,tb_thl,tb_qt,tb_nccn,tb_ug,tb_vg,tb_dqtdxls,tb_dqtdyls, &
                 tb_qtadv,tb_thladv,tb_uadv,tb_vadv,                                           &
                 tb_ql,tb_qi,tb_dcl,tb_dci,                                                    &
                 tb_qladv,tb_qiadv,                                                            &
                 tb_time,tb_lat,tb_lon,tb_ps,tb_qts,tb_thls,tb_wts,tb_wqs,                     &
                 tb_z0m,tb_z0h,tb_alb,tb_Qnet,tb_seaicefrct,tb_thls_ocean,tb_thls_seaice,      &
                 tb_tsoilav,tb_phiwav,tbrad_p,tbrad_t,tbrad_qv,tbrad_ql,tbrad_qi,tbrad_o3      &
               )
    if(ltb_nudge) then  ! deallocate nudging update profile
      deallocate(upnu,vpnu,wpnu,thlpnu)
      deallocate(qtpnu)
    endif
    if(nsv>0) then         
      deallocate(tb_sv        &
                ,tb_svadv     &
                ,tb_svs       &
               )
    end if
   end if
  end subroutine exittestbed

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
  subroutine handle_err(errcode)
      
  implicit none

  integer errcode
     
  write(6,*) 'Error: ', nf90_strerror(errcode)
  stop 2
      
  end subroutine handle_err


end module modtestbed
