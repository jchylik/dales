!>\file modthermodynamics.f90
!! Do the thermodynamics

!>
!! Do the thermodynamics
!>
!! Timeseries of the most relevant parameters. Written to tmser1.expnr and tmsurf.expnr
!! If netcdf is true, this module leads the tmser.expnr.nc output
!!  \author Pier Siebesma, K.N.M.I.
!!  \author Stephan de Roode,TU Delft
!!  \author Thijs Heus,MPI-M
!!  \par Revision list
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

module modthermodynamics
  use modprecision, only : field_r
  use modtimer
  implicit none
!   private
  public :: thermodynamics,calc_halflev
  public :: lqlnr
  logical :: lqlnr    = .true. !< switch for ql calc. with Newton-Raphson (on/off)
  real, allocatable :: th0av(:)
  real(field_r), allocatable :: thv0(:,:,:)
  real :: chi_half=0.5  !< set wet, dry or intermediate (default) mixing over the cloud edge
  real, allocatable :: thetah(:), qth(:), qlh(:)

contains

!> Allocate and initialize arrays
  subroutine initthermodynamics
    use modglobal, only : ih,i1,jh,j1,k1,tdn,tup,esatltab,esatitab,esatmtab,ttab
    use modmicrodata, only: imicro,imicro_bulk3
    implicit none
    real :: ilratio
    integer :: m

    allocate(th0av(k1))
    allocate(thv0(2-ih:i1+ih,2-jh:j1+jh,k1))
    allocate(thetah(k1), qth(k1), qlh(k1))

    th0av(:) = 0.

    !$acc enter data copyin(th0av, thv0, thetah, qth, qlh)

    ! esatltab(m) gives the saturation vapor pressure over water at T corresponding to m
    ! esatitab(m) is the same over ice
    ! esatmtab(m) is interpolated between the ice and liquid values with ilratio
    ! http://www.radiativetransfer.org/misc/atmlabdoc/atmlab/h2o/thermodynamics/e_eq_water_mk.html
    ! Murphy and Koop 2005 parameterization formula.
    do m=1,2000
       ttab(m)=150.+0.2*m
       esatltab(m)=exp(54.842763-6763.22/ttab(m)-4.21*log(ttab(m))+0.000367*ttab(m)+&
            tanh(0.0415*(ttab(m)-218.8))*(53.878-1331.22/ttab(m)-9.44523*log(ttab(m))+ 0.014025*ttab(m)))

       esatitab(m)=exp(9.550426-5723.265/ttab(m)+3.53068*log(ttab(m))-0.00728332*ttab(m))
       ilratio = max(0.,min(1.,(ttab(m)-tdn)/(tup-tdn)))
       if(imicro.eq.imicro_bulk3) then
          ! bulkmicro3 thermodynamics is for liquid only, ice is explicitely accounted for separately.
          esatmtab(m) = esatltab(m)
       else
          ! for all other microphysics, saturation is w.r.t. liquid and ice
          esatmtab(m) = ilratio*esatltab(m) + (1-ilratio)*esatitab(m)
       end if
    end do
  end subroutine initthermodynamics

!> Do moist thermodynamics.
!! Calculate the liquid water content, do the microphysics, calculate the mean hydrostatic pressure,
!! calculate the fields at the half levels, and finally calculate the virtual potential temperature.
  subroutine thermodynamics
    use modglobal,  only : lmoist,timee,k1,i1,j1,ih,jh,rd,rv,ijtot,cp,rlv,lnoclouds,lfast_thermo
    use modfields,  only : thl0, qt0, ql0, presf, exnf, thvh, thv0h, qt0av, ql0av, thvf, rhof
    use modmpi,     only : slabsum
    use modibm,     only : fluid_mask
    use modibmdata, only : lapply_ibm
    use modslabaverage, only : slabavg
    implicit none
    integer:: i, j, k

    call timer_tic('modthermodynamics/thermodynamics', 0)

    if (timee < 0.01) then
      call diagfld
    end if
    if (lmoist .and. (.not. lnoclouds)) then
       if (lfast_thermo) then
#if defined (DALES_GPU)
          call icethermo0_fast_gpu
#else
          call icethermo0_fast
#endif
       else
          call icethermo0
       end if
    else
       call calc_dry_tmp ! tmp0 is used in statistics
                         ! can consider calculating it only when needed
    end if
    call diagfld

    call calc_halflev !calculate halflevel values of qt0 and thl0

    if (lmoist .and. (.not. lnoclouds)) then
       if (lfast_thermo) then
#if defined (DALES_GPU)
          call icethermoh_fast_gpu
#else
          call icethermoh_fast
#endif
       else
          call icethermoh
       end if
    end if

    ! recalculate thv and rho on the basis of results
    call calthv

    !$acc parallel loop collapse(3) default(present) async(2)
    do k = 1, k1
      do j = 2, j1
        do i = 2, i1
          thv0(i,j,k) = (thl0(i,j,k)+rlv*ql0(i,j,k)/(cp*exnf(k))) &
                      * (1+(rv/rd-1)*qt0(i,j,k)-rv/rd*ql0(i,j,k))
        end do
      end do
    end do

    !$acc kernels default(present)
    thvh(:) = 0.0
    thvf(:) = 0.0
    !$acc end kernels

    !$acc host_data use_device(thvh, thv0h, thvf, thv0)
    if (.not. lapply_ibm) then
      call slabsum(thvh,1,k1,thv0h,2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1) ! redefine halflevel thv using calculated thv
      call slabsum(thvf,1,k1,thv0,2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1)
    else
      call slabavg(thv0h,fluid_mask,ih,thvh)
      call slabavg(thv0,fluid_mask,ih,thvf)
    end if
    !$acc end host_data

    !$acc kernels default(present) async(1)
    if (.not. lapply_ibm) then
      thvh(:) = thvh(:)/ijtot
      thvf(:) = thvf(:)/ijtot
    end if
    thvh(1) = th0av(1)*(1+(rv/rd-1)*qt0av(1)-rv/rd*ql0av(1)) ! override first level
    !$acc end kernels

    !$acc parallel loop default(present) async(1)
    do k = 1, k1
      rhof(k) = presf(k)/(rd*thvf(k)*exnf(k))
    end do

    !$acc wait
    call timer_toc('modthermodynamics/thermodynamics')
  end subroutine thermodynamics

!> Cleans up after the run
  subroutine exitthermodynamics
    implicit none
    !$acc exit data delete(th0av, thv0, thetah, qth, qlh)
    deallocate(th0av, thv0, thetah, qth, qlh)
  end subroutine exitthermodynamics

  !> Calculate real temperature tmp0 from thl0, for the dry case i.e. ql=0
  subroutine calc_dry_tmp
    use modglobal, only : i1,j1,k1
    use modfields, only : thl0,exnf
    use modfields, only : tmp0

    implicit none
    integer :: i, j, k

    !$acc parallel loop collapse(3) default(present) async
    do k = 1,k1
       do j = 2,j1
          do i = 2,i1
             tmp0(i,j,k) = exnf(k)*thl0(i,j,k)
          end do
       end do
    end do

  end subroutine calc_dry_tmp

!> Calculate thetav and dthvdz
  subroutine calthv
    use modglobal, only : lmoist,i1,j1,k1,kmax,zf,dzh,rlv,rd,rv,cp,eps1
    use modfields, only : thl0,thl0h,ql0,ql0h,qt0,qt0h,exnf,exnh,thv0h,dthvdz
    use modsurfdata,only : dthldz,dqtdz
    implicit none

    integer i, j, k
    real(field_r)    qs
    real(field_r)    a_surf,b_surf,dq,dth,dthv,temp
    real(field_r)    a_dry, b_dry, a_moist, b_moist, c_liquid, epsilon, eps_I, chi_sat, chi
    real(field_r)    del_thv_sat, del_thv_dry

    call timer_tic('modthermodynamics/calthv', 1)

    dthvdz = 0

    if (lmoist) then
      !$acc parallel loop collapse(3) default(present) async(1)
      do k = 2, k1
        do j = 2, j1
          do i = 2, i1
            thv0h(i,j,k) = (thl0h(i,j,k)+rlv*ql0h(i,j,k)/(cp*exnh(k))) &
                          *(1+(rv/rd-1)*qt0h(i,j,k)-rv/rd*ql0h(i,j,k))
          end do
        end do
      end do

      !TODO: fix the branching in this loop
      !$acc parallel loop collapse(3) default(present) &
      !$acc& private(a_dry, b_dry, a_moist, b_moist, c_liquid, epsilon, eps_I, chi_sat, chi, dthv, del_thv_dry, del_thv_sat, temp, qs, dq, dth) &
      !$acc& async(2)
      do k = 2, kmax
        do j = 2 , j1
          do i = 2, i1
!
!         default thv jump computed unsaturated
!
            epsilon = rd/rv
            eps_I = 1/epsilon - 1  !cstep approx 0.608

            a_dry = 1 + eps_I * qt0(i,j,k)
            b_dry = eps_I * thl0(i,j,k)

            dth = thl0(i,j,k+1)-thl0(i,j,k-1)
            dq  = qt0(i,j,k+1)-qt0(i,j,k-1)

            del_thv_dry = a_dry   * dth + b_dry * dq

            dthv = del_thv_dry

            if  (ql0(i,j,k)> 0) then  !include moist thermodynamics

               temp = thl0(i,j,k)*exnf(k)+(rlv/cp)*ql0(i,j,k)
               qs   = qt0(i,j,k) - ql0(i,j,k)

               a_moist = (1-qt0(i,j,k)+qs/epsilon*(1+rlv/(rv*temp))) &
                        /(1+rlv**2*qs/(cp*rv*temp**2))
               b_moist = a_moist*rlv/cp-temp
               c_liquid = a_dry * rlv / cp - thl0(i,j,k) / epsilon

               del_thv_sat = a_moist * dth + b_moist * dq

               chi     = 2*chi_half*(zf(k) - zf(k-1))/(dzh(k)+dzh(k+1))
               chi_sat = c_liquid * ql0(i,j,k) / (del_thv_dry - del_thv_sat)

               if (chi < chi_sat) then  !mixed parcel is saturated
                 dthv = del_thv_sat
              end if
            end if

            dthvdz(i,j,k) = dthv/(dzh(k+1)+dzh(k))
          end do
        end do
      end do

      !$acc parallel loop collapse(2) default(present) private(temp, qs, a_surf, b_surf) async(3)
      do j=2,j1
        do i=2,i1
          if(ql0(i,j,1)>0) then
            temp = thl0(i,j,1)*exnf(1)+(rlv/cp)*ql0(i,j,1)
            qs   = qt0(i,j,1) - ql0(i,j,1)
            a_surf   = (1-qt0(i,j,1)+rv/rd*qs*(1+rlv/(rv*temp))) &
                      /(1+rlv**2*qs/(cp*rv*temp**2))
            b_surf   = a_surf*rlv/(temp*cp)-1

          else
            a_surf = 1+(rv/rd-1)*qt0(i,j,1)
            b_surf = rv/rd-1

          end if
          dthvdz(i,j,1) = a_surf*dthldz(i,j) + b_surf*thl0(i,j,1)*dqtdz(i,j)
        end do
      end do

    else
      !$acc parallel loop collapse(3) default(present)
      do k = 2, k1
        do j = 2, j1
          do i = 2, i1
            thv0h(i,j,k)  = thl0h(i,j,k)
          end do
        end do
      end do

      !$acc parallel loop collapse(3) default(present)
      do k = 2, kmax
        do j = 2, j1
          do i = 2, i1
            dthvdz(i,j,k) = (thl0(i,j,k+1)-thl0(i,j,k-1))/(dzh(k+1)+dzh(k))
          end do
        end do
      end do

      !$acc parallel loop collapse(2) default(present)
      do j = 2, j1
        do i = 2, i1
          dthvdz(i,j,1) = dthldz(i,j)
        end do
      end do
    end if

    !$acc parallel loop collapse(3) default(present) async wait(2, 3)
    do k = 1, kmax
      do j = 2, j1
        do i = 2, i1
          if(abs(dthvdz(i,j,k)) < eps1) then
            dthvdz(i,j,k) = sign(eps1, dthvdz(i,j,k))
          end if
        end do
      end do
    end do

    !$acc wait

    call timer_toc('modthermodynamics/calthv')

  end subroutine calthv
!> Calculate diagnostic slab averaged fields.
!!     Calculates slab averaged fields assuming
!!     hydrostatic equilibrium for: u,v,theta_l,theta_v,
!!     qt,ql,exner,pressure and the density
!! \author      Pier Siebesma   K.N.M.I.     06/01/1995
  subroutine diagfld
  use modglobal,  only : i1,ih,j1,jh,k1,nsv,zh,zf,cu,cv,ijtot,grav,rlv,cp,rd,rv,pref0,timee,lconstexner
  use modfields,  only : u0,v0,thl0,qt0,ql0,sv0,u0av,v0av,thl0av,qt0av,ql0av,sv0av, &
                        presf,presh,exnf,exnh,rhof,thvf
  use modsurfdata,only : thls,ps
  use modmpi,     only : slabsum
  use modibm,     only : fluid_mask
  use modibmdata, only : lapply_ibm
  use modslabaverage, only : slabavg
  implicit none

  integer :: k,n

  call timer_tic('modthermodynamics/diagfld', 1)


!*********************************************************
!  1.0   calculate average profiles of u,v,thl,qt and ql *
!        assuming hydrostatic equilibrium                *
!*********************************************************

! initialise local MPI arrays

  !$acc kernels default(present)
  u0av = 0.0
  v0av = 0.0
  thl0av = 0.0
  th0av  = 0.0
  qt0av  = 0.0
  ql0av  = 0.0
  sv0av = 0.
  !$acc end kernels

  !CvH changed momentum array dimensions to same value as scalars!
  !$acc host_data use_device(u0av, u0, v0av, v0, thl0av, thl0, qt0av, qt0, &
  !$acc&                     ql0av, ql0, sv0av, sv0)
  if (.not. lapply_ibm) then
    call slabsum(u0av  ,1,k1,u0  ,2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1)
    call slabsum(v0av  ,1,k1,v0  ,2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1)
    call slabsum(thl0av,1,k1,thl0,2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1)
    call slabsum(qt0av ,1,k1,qt0 ,2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1)
    call slabsum(ql0av ,1,k1,ql0 ,2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1)
    do n=1,nsv
      call slabsum(sv0av(1:1,n),1,k1,sv0(:,:,:,n),2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1)
    end do
  else
    call slabavg(u0,fluid_mask,ih,u0av)
    call slabavg(v0,fluid_mask,ih,v0av)
    call slabavg(thl0,fluid_mask,ih,thl0av)
    call slabavg(qt0,fluid_mask,ih,qt0av)
    call slabavg(ql0,fluid_mask,ih,ql0av)
    do n=1,nsv
      call slabavg(sv0(:,:,:,n),fluid_mask,ih,sv0av(:,n))
    end do
  end if
  !$acc end host_data

  !$acc kernels default(present)
  if (.not. lapply_ibm) then
    u0av   = u0av  /ijtot + cu
    v0av   = v0av  /ijtot + cv
    thl0av = thl0av/ijtot
    qt0av  = qt0av /ijtot
    ql0av  = ql0av /ijtot
    sv0av  = sv0av /ijtot
  end if
  if (timee < 0.01 .or. .not. lconstexner) then
    exnf   = 1-grav*zf/(cp*thls)
    exnh   = 1-grav*zh/(cp*thls)
  endif
  th0av  = thl0av+ (rlv/cp)*ql0av/exnf
  !$acc end kernels

!***********************************************************
!  2.0   calculate average profile of pressure at full and *
!        half levels, assuming hydrostatic equilibrium.    *
!***********************************************************

!    2.1 Use first guess of theta, then recalculate theta

   call fromztop

   !$acc kernels default(present)
   th0av = thl0av + (rlv/cp)*ql0av/exnf
   if (timee < 0.01 .or. .not. lconstexner) then
      exnf = (presf/pref0)**(rd/cp)
   endif
   !$acc end kernels

!    2.2 Use new updated value of theta for determination of pressure

   call fromztop

!***********************************************************
!  3.0   Construct density profiles and exner function     *
!       for further use in the program                     *
!***********************************************************

!  3.1 determine exner
   if (timee < 0.01 .or. .not. lconstexner) then
     !$acc serial default(present) async(1)
     exnh(1) = (ps/pref0)**(rd/cp)
     exnf(1) = (presf(1)/pref0)**(rd/cp)
     !$acc end serial

     !$acc parallel loop default(present) async(2)
     do k=2,k1
       exnf(k) = (presf(k)/pref0)**(rd/cp)
       exnh(k) = (presh(k)/pref0)**(rd/cp)
     end do
   endif

!  3.2 determine rho
   !$acc parallel loop default(present) async wait(1, 2)
   do k=1,k1
     thvf(k) = th0av(k)*exnf(k)*(1+(rv/rd-1)*qt0av(k)-rv/rd*ql0av(k))
     rhof(k) = presf(k)/(rd*thvf(k))
   end do
   !$acc wait

   call timer_toc('modthermodynamics/diagfld')

   return
  end subroutine diagfld

!> Calculates slab averaged pressure
!!      Input :  zf,zh,theta and qt profile
!!      Output:  pressure profile at full and
!!               half levels
!!
!!      Method: Using hydrostatic equilibrium
!!
!!                              -g*pref0**(rd/cp)
!! =====>       dp**(rd/cp)/dz = --------------
!!                                 cp*thetav
!! \author Pier Siebesma   K.N.M.I.     06/01/1995
  subroutine fromztop

  use modglobal, only : k1,dzf,dzh,rv,rd,cp,zf,grav,pref0
  use modfields, only : qt0av,ql0av,presf,presh,thvh,thvf
  use modsurfdata,only : ps
  implicit none

  integer   k
  real(field_r)  rdocp

  call timer_tic('modthermodynamics/fromztop', 1)

  rdocp = rd/cp

!**************************************************
!    1.0 Determine theta and qt at half levels    *
!**************************************************

  !$acc parallel loop default(present)
  do k=2,k1
    thetah(k) = (th0av(k)*dzf(k-1) + th0av(k-1)*dzf(k))/(2*dzh(k))
    qth   (k) = (qt0av(k)*dzf(k-1) + qt0av(k-1)*dzf(k))/(2*dzh(k))
    qlh   (k) = (ql0av(k)*dzf(k-1) + ql0av(k-1)*dzf(k))/(2*dzh(k))
  end do

!**************************************************
!     2.1  calculate pressures at full levels     *
!          assuming hydrostatic equilibrium       *
!**************************************************

!     1: lowest level: use first level value for safety!

  !$acc update self(thetah, qth, qlh, th0av, qt0av, ql0av)

  thvh(1) = th0av(1)*(1+(rv/rd-1)*qt0av(1)-rv/rd*ql0av(1))
  presf(1) = ps**rdocp - grav*(pref0**rdocp)*zf(1) /(cp*thvh(1))
  presf(1) = presf(1)**(1/rdocp)

!     2: higher levels

  do k=2,k1
    thvh(k)  = thetah(k)*(1+(rv/rd-1)*qth(k)-rv/rd*qlh(k))
    presf(k) = presf(k-1)**rdocp - &
                   grav*(pref0**rdocp)*dzh(k) /(cp*thvh(k))
    presf(k) = presf(k)**(1/rdocp)
  end do

!**************************************************
!     2.2   calculate pressures at half levels    *
!           assuming hydrostatic equilibrium      *
!**************************************************

  presh(1) = ps
  thvf(1) = th0av(1)*(1+(rv/rd-1)*qt0av(1)-rv/rd*ql0av(1))

  do k=2,k1
    thvf(k)  = th0av(k)*(1+(rv/rd-1)*qt0av(k)-rv/rd*ql0av(k))
    presh(k) = presh(k-1)**rdocp - &
                   grav*(pref0**rdocp)*dzf(k-1) / (cp*thvf(k-1))
    presh(k) = presh(k)**(1/rdocp)
  end do

  !$acc update device(thvh, presf, thvf, presh)
  call timer_toc('modthermodynamics/fromztop')

  return
  end subroutine fromztop

!> Calculates liquid water content.
!!     Given theta_l and q_tot the liquid water content
!!     is calculated, making an "all-or-nothing" assumption.
!!     if lfull=true   ==> ql at full levels on output
!!     if lfull=false  ==> ql at half levels on output
!!
!! \author Hans Cuijpers   I.M.A.U.
!! \author Pier Siebesma   K.N.M.I.     06/01/1995
  subroutine thermo (thl,qt,ql,pressure,exner)



  use modglobal, only : ih,jh,i1,j1,k1,es0,at,bt,rd,rv,rlv,cp,tmelt
  implicit none

  integer i, j, k
  real tl, es, qs, qsl, b1
  real, intent(in)  :: qt(2-ih:i1+ih,2-jh:j1+jh,k1),thl(2-ih:i1+ih,2-jh:j1+jh,k1),exner(k1),pressure(k1)
  real, intent(out) :: ql(2-ih:i1+ih,2-jh:j1+jh,k1)
  real :: Tnr,qsatur,Tnr_old
  integer :: niter,nitert
    if (lqlnr) then

!mc      calculation of T with Newton-Raphson method
!mc      first guess is Tnr=tl
!mc
      nitert = 0
      do k=1,k1
        do j=2,j1
          do i=2,i1

            tl  = thl(i,j,k)*exner(k)
            Tnr=tl
            Tnr_old=0
            do while (abs(Tnr-Tnr_old)/Tnr>1e-5)
              niter = niter+1
              Tnr_old = Tnr
              es    = es0*exp(at*(Tnr-tmelt)/(Tnr-bt))
              qsatur= rd/rv*es/(pressure(k)-(1-rd/rv)*es)
              Tnr = Tnr - (Tnr+(rlv/cp)*qsatur-tl- &
                      (rlv/cp)*qt(i,j,k))/(1+(rlv**2*qsatur)/ &
                      (rv*cp*Tnr**2))
            end do
            nitert =max(nitert,niter)
            niter = 0

            ql(i,j,k) = dim(qt(i,j,k)-qsatur,0.)

          end do
        end do
      end do
    else


      do k=1,k1
        do j=2,j1
          do i=2,i1
            tl  = thl(i,j,k)*exner(k)
            es  = es0*exp(at*(tl-tmelt)/(tl-bt))
            qsl = rd/rv*es/(pressure(k)-(1-rd/rv)*es)
            b1  = rlv**2/(tl**2*cp*rv)
            qs  = qsl*(1+b1*qt(i,j,k))/(1.+b1*qsl)
            ql(i,j,k) = dim(qt(i,j,k)-qs,0.)
          end do
        end do
      end do
    end if

  return
  end subroutine thermo

!!!!!!!!! new thermo
!> Magnus formulas for q_sat over liquid and ice
!> from Huang 2018 https://doi.org/10.1175/JAMC-D-17-0334.
!> Warning: for performance, check that rd/rv etc are pre-computed
  pure function qsat_magnus(T, p) result(qsat)
    use modglobal, only : rd,rv,tup,tdn
    implicit none
    real(field_r), intent(in) :: T, p
    real :: qsat
    real ilratio, TC, esl, esi, es
    ilratio = max(0._field_r,min(1._field_r,(T-tdn)/(tup-tdn)))

    TC = T - 273.15 ! in Celcius
    esl = 610.94_field_r * exp( (17.625_field_r*TC) / (TC+243.04_field_r) ) ! Magnus
    esi = 611.21_field_r * exp( (22.587_field_r*TC) / (TC+273.86_field_r) ) ! Magnus

    ! interpolated saturation vapor pressure
    es = ilratio*esl + (1-ilratio)*esi

    ! convert saturation vapor pressure to saturation humidity
    qsat = (rd/rv) * es / (p - (1-rd/rv)*es)
  end function qsat_magnus

!> Huang's formulas for q_sat over liquid and ice
!> from Huang 2018 https://doi.org/10.1175/JAMC-D-17-0334.
!> should be more accurate than Magnus, at the cost of more divisions
!> Warning: for performance, check that rd/rv etc are pre-computed
  pure function qsat_huang(T, p) result(qsat)
    use modglobal, only : rd,rv,tup,tdn
    implicit none
    real(field_r), intent(in) :: T, p
    real :: qsat
    real ilratio, TC, esl, esi, es
    ilratio = max(0._field_r,min(1._field_r,(T-tdn)/(tup-tdn)))

    TC = T - 273.15_field_r ! in Celcius
    esl = exp(34.494_field_r - 4924.99_field_r / (TC  + 237.1_field_r)) /  (TC+105)**1.57_field_r  ! Huang
    esi = exp(43.494_field_r - 6545.8_field_r/(TC+278)) / (TC+868)**2              ! Huang

    ! interpolated saturation vapor pressure
    es = ilratio*esl + (1-ilratio)*esi

    ! convert saturation vapor pressure to saturation humidity
    qsat = (rd/rv) * es / (p - (1-rd/rv)*es)
  end function qsat_huang

  ! return esat for ice-liquid mix using table
  pure function esat_tab(T) result(es)
    use modglobal, only : esatmtab

    implicit none
    !$acc routine seq
    real(field_r), intent(in) :: T
    integer :: tlonr
    real(field_r) :: tlo, thi, es

    ! interpolated ice-liquid saturation vapor pressure from table
    ! note if imicto==imicro_bulk3, the table is for liquid only
    tlonr=int((T-150)*5)
    tlo = 150 + 0.2_field_r*tlonr
    thi = tlo + 0.2_field_r
    es = (thi-T)*5*esatmtab(tlonr)+(T-tlo)*5*esatmtab(tlonr+1)
  end function esat_tab

!> q_sat over liquid and ice, using interpolation in a table created in modglobal.
!> seems to be faster than the Magnus formula (on CPU)
  pure function qsat_tab(T, p) result(qsat)
    use modglobal, only : rd,rv
    use modglobal, only : esatmtab

    implicit none
    !$acc routine seq
    real(field_r), intent(in) :: T, p
    real(field_r) :: qsat
    integer :: tlonr
    real(field_r) :: tlo, thi, es

    ! interpolated ice-liquid saturation vapor pressure from table
    tlonr=int((T-150)*5)
    tlo = 150 + 0.2_field_r*tlonr
    thi = tlo + 0.2_field_r
    es = (thi-T)*5*esatmtab(tlonr)+(T-tlo)*5*esatmtab(tlonr+1)

    ! convert saturation vapor pressure to saturation humidity
    qsat = (rd/rv) * es / (p - (1-rd/rv)*es)
  end function qsat_tab

  subroutine icethermo0_fast
    !> Calculates liquid water content ql from thl0 and qt0.
    !> Using 2 iterations of Eq. (59) in the Heus 2010 article
    !> and e_sat interpolated between liquid and ice expressions.
    !>
    !> Given thl0 and qt0, we search for T such that
    !> (1) ql = qt - qsat(T)     (definition of ql, instant condensation if above saturation)
    !> (2) ql = cp/L * (T - Tl)  (definition of Tl)
    !> hold simultaneously, and solve for qsat(T).
    !> Tl is thl0/exnf(k) .
    !>
    !> Steps of the derivation:
    !> - 1st order Taylor expansion of qsat(T) around T = Tl
    !> - insert (1) and (2)
    !> - solve for qsat(T)
    !> - use the Clausius-Clapeyron relation for the T-derivative of qsat
    !>
    !> 2 iterations gives a good accuracy, 1 iteration is not
    !> sufficient.  Fixing the number of iterations makes the code
    !> vectorize, a variable number of iterations prevents
    !> vectorization.
    !>
    !> The procedure works also when qsat(T) is a linear interoplation
    !> between qsat_liquid and qsat_ice, with slightly reduced
    !> accuracy in the interpolation region.
    !>
    !> qsat (T) can be calculated in different ways, with different
    !> accuracy vs computing cost.  The fastest so far is to use a
    !> lookup table for esat(T), and interpolate linearly in it.
    !>
    !> \author Fredrik Jansson, Jisk Attema, Pier Siebesma

    use modglobal, only : i1,j1,k1,rv,rlv,cp,rd
    use modfields, only : qt0,thl0,exnf,presf,ql0
    use modfields, only : tmp0, qsat, esl, qvsl, qvsi          ! consider not storing these
    use modglobal, only : esatltab, esatitab

    implicit none
    integer :: i, j, k
    real(field_r) :: Tl, qsat_, qt, ql, b, T
    real(field_r) :: Tl_min, Tl_max, qt_max
    real(field_r) :: esi1, tlo, thi
    integer       :: tlonr
    do k=1,k1
       ! Optimization: if the whole horizontal slab at k is unsaturated,
       ! the calculation of this slab can be skipped.
       ! Find highest qt and lowest thl in the slab.
       ! If they in combination are not saturated, the whole slab is below saturation.
       ! Also do range checks of Tl here. Tl must be within the range of the table,
       ! and below the boiling point of water at this level.
       ! Setting the limit at 5K below the boiling point here. Crossing the boiling point
       ! is detected by esat > presf(k)
       Tl_min = minval(thl0(2:i1,2:j1,k)) * exnf(k)
       Tl_max = maxval(thl0(2:i1,2:j1,k)) * exnf(k)
       qt_max = maxval(qt0(2:i1,2:j1,k))
       if (Tl_min < 150) STOP 'icethermo0_fast: Tl_min below limit 150K'
       if (esat_tab(Tl_max + 5) > presf(k)) STOP 'icethermo0_fast: Tl_max too close to boiling point'

       qsat_ = qsat_tab(Tl_min, presf(k)) ! lowest possible qsat in this slab
       if (qt_max > qsat_) then
          do j=2,j1
             do i=2,i1
                Tl = exnf(k)*thl0(i,j,k)
                qt = qt0(i,j,k)

                ! first step
                qsat_ = qsat_tab(Tl, presf(k))
                b = rlv**2 / (rv * cp * Tl**2)
                qsat_ = qsat_ * (1 + b * qt) / (1 + b * qsat_)

                ql = max(qt0(i,j,k) - qsat_, 0._field_r)

                ! update the starting point
                Tl = Tl + (rlv/cp) * ql
                qt = qt - ql

                ! second step
                qsat_ = qsat_tab(Tl, presf(k))
                b = rlv**2 / (rv * cp * Tl**2)
                qsat_ = qsat_ * (1 + b * qt) / (1 + b * qsat_)

                ! save results
                ql = max(qt0(i,j,k) - qsat_, 0._field_r)
                ql0(i,j,k) = ql

                !!!!!!!!!!!!!!!!!
                ! The following could
                ! be done on the fly to save
                ! precious memory

                !qsat(i,j,k) = qsat_ ! qsat_ is not a good approximation when not saturated
                                     ! but ql is still good in that case.
                T = exnf(k)*thl0(i,j,k) + (rlv/cp) * ql
                tmp0(i,j,k) = T

                ! use the separate e_sat tables for liquid and ice to calculate and store esl, qvsl, qvsi
                tlonr=int((T-150)*5)
                tlo = 150 + 0.2_field_r*tlonr
                thi = tlo + 0.2_field_r
                esl(i,j,k) = (thi-T)*5*esatltab(tlonr)+(T-tlo)*5*esatltab(tlonr+1) ! saturation vapor pressure liquid
                esi1       = (thi-T)*5*esatitab(tlonr)+(T-tlo)*5*esatitab(tlonr+1) ! saturation vapor pressure ice
                qvsl(i,j,k)=rd/rv*esl(i,j,k)/(presf(k)-(1-rd/rv)*esl(i,j,k))        ! saturation humidity liquid
                qvsi(i,j,k)=rd/rv*esi1      /(presf(k)-(1-rd/rv)*esi1)              ! saturation humidity ice
                !!!!!!!!!!!
             end do
          end do
       else
          ! possibly faster option when the whole layer is below saturation
          ! If many of the els, qsvl, qsvi are stored in arrays, they still need to be saved here
          ql0(2:i1,2:j1,k) = 0
          do j=2,j1
             do i=2,i1
                T = exnf(k)*thl0(i,j,k) ! + (rlv/cp) * ql omitted because ql is 0
                tmp0(i,j,k) = T
                qsat(i,j,k) = qsat_tab(T, presf(k))

                ! use the separate e_sat tables for liquid and ice to calculate and store esl, qvsl, qvsi
                tlonr=int((T-150)*5)
                tlo = 150 + 0.2_field_r*tlonr
                thi = tlo + 0.2_field_r
                esl(i,j,k) = (thi-T)*5*esatltab(tlonr)+(T-tlo)*5*esatltab(tlonr+1) ! saturation vapor pressure liquid
                esi1       = (thi-T)*5*esatitab(tlonr)+(T-tlo)*5*esatitab(tlonr+1) ! saturation vapor pressure ice
                qvsl(i,j,k)=rd/rv*esl(i,j,k)/(presf(k)-(1-rd/rv)*esl(i,j,k))        ! saturation humidity liquid
                qvsi(i,j,k)=rd/rv*esi1      /(presf(k)-(1-rd/rv)*esi1)              ! saturation humidity ice
             end do
          end do
       end if
    end do
  end subroutine icethermo0_fast

  subroutine icethermo0_fast_gpu
    !> Calculates liquid water content ql from thl0 and qt0.
    !> Using 2 iterations of Eq. (59) in the Heus 2010 article
    !> and e_sat interpolated between liquid and ice expressions.
    !>
    !> Given thl0 and qt0, we search for T such that
    !> (1) ql = qt - qsat(T)     (definition of ql, instant condensation if above saturation)
    !> (2) ql = cp/L * (T - Tl)  (definition of Tl)
    !> hold simultaneously, and solve for qsat(T).
    !> Tl is thl0/exnf(k) .
    !>
    !> Steps of the derivation:
    !> - 1st order Taylor expansion of qsat(T) around T = Tl
    !> - insert (1) and (2)
    !> - solve for qsat(T)
    !> - use the Clausius-Clapeyron relation for the T-derivative of qsat
    !>
    !> 2 iterations gives a good accuracy, 1 iteration is not
    !> sufficient.  Fixing the number of iterations makes the code
    !> vectorize, a variable number of iterations prevents
    !> vectorization.
    !>
    !> The procedure works also when qsat(T) is a linear interoplation
    !> between qsat_liquid and qsat_ice, with slightly reduced
    !> accuracy in the interpolation region.
    !>
    !> qsat (T) can be calculated in different ways, with different
    !> accuracy vs computing cost.  The fastest so far is to use a
    !> lookup table for esat(T), and interpolate linearly in it.
    !>
    !> C. Jungbacker: This version is faster on GPU than the original
    !> because it does not check for each slab.
    !> \author Fredrik Jansson, Jisk Attema, Pier Siebesma

    use modglobal, only : i1,j1,k1,rv,rlv,cp,rd
    use modfields, only : qt0,thl0,exnf,presf,ql0
    use modfields, only : tmp0, qsat, esl, qvsl, qvsi          ! consider not storing these
    use modglobal, only : esatltab, esatitab

    implicit none
    integer :: i, j, k
    real(field_r) :: Tl, qsat_, qt, ql, b, T
    real(field_r) :: Tl_min, PrDiff_min
    real(field_r) :: esi1, tlo, thi
    integer       :: tlonr

    call timer_tic('modthermodynamics/icethermo0_fast', 1)

    ! Sanity checks
    Tl_min = 400
    PrDiff_min = 100
    !$acc parallel loop collapse(3) reduction(min:Tl_min, PrDiff_min)
    do k = 1, k1
      do j = 2, j1
        do i = 2, i1
          Tl_min = min(Tl_min,thl0(i,j,k)*exnf(k))
          PrDiff_min = min(PrDiff_min, presf(k) - esat_tab(thl0(i,j,k)*exnf(k) + 5.0_field_r))
        end do
      end do
    end do
    if (Tl_min < 150) stop 'icethermo0_fast: Tl_min below limit 150K'
    if (PrDiff_min < 0) stop 'icethermo0_fast: Tl_max too close to boiling point'

    !$acc parallel loop collapse(3) private(Tl, qsat_, qt, ql, b, T, esi1, tlo, thi, tlonr) default(present)
    do k = 1, k1
      do j = 2, j1
        do i = 2, i1
          qsat_ = qsat_tab(Tl_min, presf(k))
          Tl = exnf(k)*thl0(i,j,k)
          qt = qt0(i,j,k)

          ! first step
          qsat_ = qsat_tab(Tl, presf(k))
          b = rlv**2 / (rv * cp * Tl**2)
          qsat_ = qsat_ * (1 + b * qt) / (1 + b * qsat_)

          ql = max(qt0(i,j,k) - qsat_, 0._field_r)

          ! update the starting point
          Tl = Tl + (rlv/cp) * ql
          qt = qt - ql

          ! second step
          qsat_ = qsat_tab(Tl, presf(k))
          b = rlv**2 / (rv * cp * Tl**2)
          qsat_ = qsat_ * (1 + b * qt) / (1 + b * qsat_)

          ! save results
          ql = max(qt0(i,j,k) - qsat_, 0._field_r)
          ql0(i,j,k) = ql

          !!!!!!!!!!!!!!!!!
          ! The following could
          ! be done on the fly to save
          ! precious memory
          qsat(i,j,k) = qsat_
          T = exnf(k)*thl0(i,j,k) + (rlv/cp) * ql
          tmp0(i,j,k) = T

          ! use the separate e_sat tables for liquid and ice to calculate and store esl, qvsl, qvsi
          tlonr=int((T-150)*5)
          tlo = 150 + 0.2_field_r*tlonr
          thi = tlo + 0.2_field_r
          esl(i,j,k) = (thi-T)*5*esatltab(tlonr)+(T-tlo)*5*esatltab(tlonr+1) ! saturation vapor pressure liquid
          esi1       = (thi-T)*5*esatitab(tlonr)+(T-tlo)*5*esatitab(tlonr+1) ! saturation vapor pressure ice
          qvsl(i,j,k)=rd/rv*esl(i,j,k)/(presf(k)-(1-rd/rv)*esl(i,j,k))        ! saturation humidity liquid
          qvsi(i,j,k)=rd/rv*esi1      /(presf(k)-(1-rd/rv)*esi1)              ! saturation humidity ice
          !!!!!!!!!!!
        end do
      end do
    end do
    call timer_toc('modthermodynamics/icethermo0_fast')
  end subroutine icethermo0_fast_gpu

  subroutine icethermoh_fast
    !> Calculates liquid water content ql for halflevels
    !> Using 2 iterations of Eq. (59) in the Heus 2010 article
    !> and e_sat interpolated between liquid and ice expressions.
    !> See comments in icethermo0_fast above for more details.
    !>
    !> \author Fredrik Jansson, Jisk Attema, Pier Siebesma
    !> this could be merged with icethermo0
    !> and input and output fields are given as parameters.
    !> in: thl, qt, exner,
    !> out: ql
    !>
    !> alternatively merge with calc_halflev and calthv
    !> to eliminate qt0h, thl0h, ql0h fields

    use modglobal, only : i1,j1,k1,rv,rlv,cp
    use modfields, only : qt0h,thl0h,exnh,presh,ql0h

    implicit none
    integer :: i, j, k
    real(field_r) :: Tl, qsat, qt, ql, b
    real(field_r) :: Tl_min, Tl_max, qt_max

    do k=1,k1
       ! find highest qt and lowest thl in the slab.
       ! if they in combination are not saturated, the whole slab is below saturation
       Tl_min = minval(thl0h(2:i1,2:j1,k)) * exnh(k)
       Tl_max = maxval(thl0h(2:i1,2:j1,k)) * exnh(k)
       if (Tl_min < 150) STOP 'icethermoh_fast: Tl_min below limit 150K'
       if (esat_tab(Tl_max + 5) > presh(k)) STOP 'icethermoh_fast: Tl_max too close to boiling point'
       qt_max = maxval(qt0h(2:i1,2:j1,k))
       qsat = qsat_tab(Tl_min, presh(k))
       if (qt_max > qsat) then
          do j=2,j1
             do i=2,i1
                Tl = exnh(k)*thl0h(i,j,k)
                qt = qt0h(i,j,k)

                ! first step
                qsat = qsat_tab(Tl, presh(k))
                b = rlv**2 / (rv * cp * Tl**2)
                qsat = qsat * (1 + b * qt) / (1 + b * qsat)

                ql = max(qt0h(i,j,k) - qsat, 0._field_r)

                ! update the starting point
                Tl = Tl + (rlv/cp) * ql
                qt = qt - ql

                ! second step
                qsat = qsat_tab(Tl, presh(k))
                b = rlv**2 / (rv * cp * Tl**2)
                qsat = qsat * (1 + b * qt) / (1 + b * qsat)

                ! save results
                ql = max(qt0h(i,j,k) - qsat, 0._field_r)
                ql0h(i,j,k) = ql
             end do
          end do
       else
          ql0h(2:i1,2:j1,k) = 0
       end if
    end do
  end subroutine icethermoh_fast

  subroutine icethermoh_fast_gpu
    !> Calculates liquid water content ql for halflevels
    !> Using 2 iterations of Eq. (59) in the Heus 2010 article
    !> and e_sat interpolated between liquid and ice expressions.
    !> See comments in icethermo0_fast above for more details.
    !>
    !> C. Jungbacker: This version is faster on GPU than the original
    !> because it does not check for each slab.
    !> \author Fredrik Jansson, Jisk Attema, Pier Siebesma
    !>
    !> this could be merged with icethermo0
    !> and input and output fields are given as parameters.
    !> in: thl, qt, exner,
    !> out: ql
    !>
    !> alternatively merge with calc_halflev and calthv
    !> to eliminate qt0h, thl0h, ql0h fields

    use modglobal, only : i1, j1, k1, rv, rlv, cp
    use modfields, only : qt0h, thl0h, exnh, presh, ql0h

    implicit none
    integer :: i, j, k
    real(field_r) :: Tl, qsat, qt, ql, b
    real(field_r) :: Tl_min, PrDiff_min

    call timer_tic('modthermodynamics/icethermoh_fast', 1)

    Tl_min = 400
    PrDiff_min = 100
    !$acc parallel loop collapse(3) reduction(min:Tl_min, PrDiff_min)
    do k = 1, k1
      do j = 2, j1
        do i = 2, i1
          Tl_min = min(Tl_min,thl0h(i,j,k)*exnh(k))
          PrDiff_min = min(PrDiff_min, presh(k) - esat_tab(thl0h(i,j,k)*exnh(k) + 5.0_field_r))
        end do
      end do
    end do
    if (Tl_min < 150) stop 'icethermoh_fast: Tl_min below limit 150K'
    if (PrDiff_min < 0.0) stop 'icethermoh_fast: Tl_max too close to boiling point'

    !$acc parallel loop collapse(3) default(present) private(Tl, qt, qsat, b, ql)
    do k = 1, k1
      do j = 2, j1
        do i = 2, i1
          Tl = exnh(k)*thl0h(i,j,k)
          qt = qt0h(i,j,k)

          ! first step
          qsat = qsat_tab(Tl, presh(k))
          b = rlv**2 / (rv * cp * Tl**2)
          qsat = qsat * (1 + b * qt) / (1 + b * qsat)

          ql = max(qt0h(i,j,k) - qsat, 0._field_r)

          ! update the starting point
          Tl = Tl + (rlv/cp) * ql
          qt = qt - ql

          ! second step
          qsat = qsat_tab(Tl, presh(k))
          b = rlv**2 / (rv * cp * Tl**2)
          qsat = qsat * (1 + b * qt) / (1 + b * qsat)

          ! save results
          ql = max(qt0h(i,j,k) - qsat, 0._field_r)
          ql0h(i,j,k) = ql
        end do
      end do
    end do
    call timer_toc('modthermodynamics/icethermoh_fast')
  end subroutine icethermoh_fast_gpu
!!!!!!!!! new thermo

  subroutine icethermo0
!> Calculates liquid water content.and temperature
!! \author Steef B\"oing

  use modglobal, only : i1,j1,k1,rd,rv,rlv,tup,tdn,cp,ttab,esatltab,esatitab
  use modfields, only : qvsl,qvsi,qt0,thl0,exnf,presf,tmp0,ql0,esl,qsat
  implicit none

  integer i, j, k
  real :: ilratio, esl1,esi1,qvsl1,qvsi1,qsatur, thlguess, thlguessmin,tlo,thi,ttry
  real :: Tnr,Tnr_old
  integer :: niter,nitert,tlonr,thinr

!     calculation of T with Newton-Raphson method
!     first guess is Tnr=tl
      nitert = 0
      niter = 0
      do k=1,k1
      do j=2,j1
      do i=2,i1
            ! first guess for temperature
            Tnr=exnf(k)*thl0(i,j,k)
            ilratio = max(0.,min(1.,(Tnr-tdn)/(tup-tdn)))
            tlonr=int((Tnr-150.)*5.)
            thinr=tlonr+1
            tlo=ttab(tlonr)
            thi=ttab(thinr)
            esl1=(thi-Tnr)*5.*esatltab(tlonr)+(Tnr-tlo)*5.*esatltab(thinr)
            esi1=(thi-Tnr)*5.*esatitab(tlonr)+(Tnr-tlo)*5.*esatitab(thinr)
            qvsl1=(rd/rv)*esl1/(presf(k)-(1.-rd/rv)*esl1)
            qvsi1=(rd/rv)*esi1/(presf(k)-(1.-rd/rv)*esi1)
            qsatur = ilratio*qvsl1+(1.-ilratio)*qvsi1
            if(qt0(i,j,k)>qsatur) then
              Tnr_old=0.
              niter = 0
              thlguess = Tnr/exnf(k)-(rlv/(cp*exnf(k)))*max(qt0(i,j,k)-qsatur,0.)
              ttry=Tnr-0.002
              ilratio = max(0.,min(1.,(ttry-tdn)/(tup-tdn)))
              tlonr=int((ttry-150.)*5.)
              thinr=tlonr+1
              tlo=ttab(tlonr)
              thi=ttab(thinr)
              esl1=(thi-ttry)*5.*esatltab(tlonr)+(ttry-tlo)*5.*esatltab(thinr)
              esi1=(thi-ttry)*5.*esatitab(tlonr)+(ttry-tlo)*5.*esatitab(thinr)
              qsatur = ilratio*(rd/rv)*esl1/(presf(k)-(1.-rd/rv)*esl1)+(1.-ilratio)*(rd/rv)*esi1/(presf(k)-(1.-rd/rv)*esi1)
              thlguessmin = ttry/exnf(k)-(rlv/(cp*exnf(k)))*max(qt0(i,j,k)-qsatur,0.)

              Tnr = Tnr - (thlguess-thl0(i,j,k))/((thlguess-thlguessmin)*500.)
              do while ((abs(Tnr-Tnr_old) > 0.002).and.(niter<100))
                niter = niter+1
                Tnr_old=Tnr
                ilratio = max(0.,min(1.,(Tnr-tdn)/(tup-tdn)))
                tlonr=int((Tnr-150.)*5.)
                if(tlonr<1 .or.tlonr>1999) then
                  write(*,*) 'thermo crash: i,j,k,niter,thl0(i,j,k),qt0(i,j,k)'
                  write(*,*) i,j,k,niter,thl0(i,j,k),qt0(i,j,k)
                endif
                thinr=tlonr+1
                tlo=ttab(tlonr)
                thi=ttab(thinr)
                esl1=(thi-Tnr)*5.*esatltab(tlonr)+(Tnr-tlo)*5.*esatltab(thinr)
                esi1=(thi-Tnr)*5.*esatitab(tlonr)+(Tnr-tlo)*5.*esatitab(thinr)
                qsatur = ilratio*(rd/rv)*esl1/(presf(k)-(1.-rd/rv)*esl1)+(1.-ilratio)*(rd/rv)*esi1/(presf(k)-(1.-rd/rv)*esi1)
                thlguess = Tnr/exnf(k)-(rlv/(cp*exnf(k)))*max(qt0(i,j,k)-qsatur,0.)

                ttry=Tnr-0.002
                ilratio = max(0.,min(1.,(ttry-tdn)/(tup-tdn)))
                tlonr=int((ttry-150.)*5.)
                thinr=tlonr+1
                tlo=ttab(tlonr)
                thi=ttab(thinr)
                esl1=(thi-ttry)*5.*esatltab(tlonr)+(ttry-tlo)*5.*esatltab(thinr)
                esi1=(thi-ttry)*5.*esatitab(tlonr)+(ttry-tlo)*5.*esatitab(thinr)
                qsatur = ilratio*(rd/rv)*esl1/(presf(k)-(1.-rd/rv)*esl1)+(1.-ilratio)*(rd/rv)*esi1/(presf(k)-(1.-rd/rv)*esi1)
                thlguessmin = ttry/exnf(k)-(rlv/(cp*exnf(k)))*max(qt0(i,j,k)-qsatur,0.)

                Tnr = Tnr - (thlguess-thl0(i,j,k))/((thlguess-thlguessmin)*500.)
              enddo
              nitert =max(nitert,niter)
              tmp0(i,j,k)= Tnr
              ilratio = max(0.,min(1.,(Tnr-tdn)/(tup-tdn)))
              tlonr=int((Tnr-150.)*5.)
              thinr=tlonr+1
              tlo=ttab(tlonr)
              thi=ttab(thinr)
              esl(i,j,k)=(thi-Tnr)*5.*esatltab(tlonr)+(Tnr-tlo)*5.*esatltab(thinr)
              esi1=(thi-Tnr)*5.*esatitab(tlonr)+(Tnr-tlo)*5.*esatitab(thinr)
              qvsl(i,j,k)=rd/rv*esl(i,j,k)/(presf(k)-(1.-rd/rv)*esl(i,j,k))
              qvsi(i,j,k)=rd/rv*esi1/(presf(k)-(1.-rd/rv)*esi1)
              qsatur = ilratio*qvsl(i,j,k)+(1.-ilratio)*qvsi(i,j,k)
            else
              tmp0(i,j,k)= Tnr
              esl(i,j,k)=esl1
              esi1=esi1
              qvsl(i,j,k)=qvsl1
              qvsi(i,j,k)=qvsi1
            endif
            ql0(i,j,k) = max(qt0(i,j,k)-qsatur,0.)
            qsat(i,j,k) = qsatur
      end do
      end do
      end do
      if(nitert>99) then
        write(*,*) 'thermowarning'
      endif

  end subroutine icethermo0

  subroutine icethermoh
!> Calculates liquid water content.and temperature
!! \author Steef B\"oing

  use modglobal, only : i1,j1,k1,rd,rv,rlv,tup,tdn,cp,ttab,esatltab,esatitab
  use modfields, only : qt0h,thl0h,exnh,presh,ql0h
  implicit none

  integer i, j, k
  real :: ilratio, esl1, esi1, qvsl1,qvsi1, qsatur, thlguess, thlguessmin,tlo,thi,ttry
  real :: Tnr,Tnr_old
  integer :: niter,nitert,tlonr,thinr

!     calculation of T with Newton-Raphson method
!     first guess is Tnr=tl
      nitert = 0
      niter = 0
      do k=1,k1
      do j=2,j1
      do i=2,i1
            ! first guess for temperature
            Tnr=exnh(k)*thl0h(i,j,k)
            ilratio = max(0.,min(1.,(Tnr-tdn)/(tup-tdn)))
            tlonr=int((Tnr-150.)*5.)
            thinr=tlonr+1
            tlo=ttab(tlonr)
            thi=ttab(thinr)
            esl1=(thi-Tnr)*5.*esatltab(tlonr)+(Tnr-tlo)*5.*esatltab(thinr)
            esi1=(thi-Tnr)*5.*esatitab(tlonr)+(Tnr-tlo)*5.*esatitab(thinr)
            qvsl1=(rd/rv)*esl1/(presh(k)-(1.-rd/rv)*esl1)
            qvsi1=(rd/rv)*esi1/(presh(k)-(1.-rd/rv)*esi1)
            qsatur = ilratio*qvsl1+(1.-ilratio)*qvsi1
            if(qt0h(i,j,k)>qsatur) then
              Tnr_old=0.
              niter = 0
              thlguess = Tnr/exnh(k)-(rlv/(cp*exnh(k)))*max(qt0h(i,j,k)-qsatur,0.)
              ttry=Tnr-0.002
              ilratio = max(0.,min(1.,(ttry-tdn)/(tup-tdn)))
              tlonr=int((ttry-150.)*5.)
              thinr=tlonr+1
              tlo=ttab(tlonr)
              thi=ttab(thinr)
              esl1=(thi-ttry)*5.*esatltab(tlonr)+(ttry-tlo)*5.*esatltab(thinr)
              esi1=(thi-ttry)*5.*esatitab(tlonr)+(ttry-tlo)*5.*esatitab(thinr)
              qsatur = ilratio*(rd/rv)*esl1/(presh(k)-(1.-rd/rv)*esl1)+(1.-ilratio)*(rd/rv)*esi1/(presh(k)-(1.-rd/rv)*esi1)
              thlguessmin = ttry/exnh(k)-(rlv/(cp*exnh(k)))*max(qt0h(i,j,k)-qsatur,0.)

              Tnr = Tnr - (thlguess-thl0h(i,j,k))/((thlguess-thlguessmin)*500.)
              do while ((abs(Tnr-Tnr_old) > 0.002).and.(niter<100))
                niter = niter+1
                Tnr_old=Tnr
                ilratio = max(0.,min(1.,(Tnr-tdn)/(tup-tdn)))
                tlonr=int((Tnr-150.)*5.)
                if(tlonr<1 .or.tlonr>1999) then
                  write(*,*) 'thermo crash: i,j,k,niter,thl0h(i,j,k),qt0h(i,j,k)'
                  write(*,*) i,j,k,niter,thl0h(i,j,k),qt0h(i,j,k)
                endif
                thinr=tlonr+1
                tlo=ttab(tlonr)
                thi=ttab(thinr)
                esl1=(thi-Tnr)*5.*esatltab(tlonr)+(Tnr-tlo)*5.*esatltab(thinr)
                esi1=(thi-Tnr)*5.*esatitab(tlonr)+(Tnr-tlo)*5.*esatitab(thinr)
                qsatur = ilratio*(rd/rv)*esl1/(presh(k)-(1.-rd/rv)*esl1)+(1.-ilratio)*(rd/rv)*esi1/(presh(k)-(1.-rd/rv)*esi1)
                thlguess = Tnr/exnh(k)-(rlv/(cp*exnh(k)))*max(qt0h(i,j,k)-qsatur,0.)

                ttry=Tnr-0.002
                ilratio = max(0.,min(1.,(ttry-tdn)/(tup-tdn)))
                tlonr=int((ttry-150.)*5.)
                thinr=tlonr+1
                tlo=ttab(tlonr)
                thi=ttab(thinr)
                esl1=(thi-ttry)*5.*esatltab(tlonr)+(ttry-tlo)*5.*esatltab(thinr)
                esi1=(thi-ttry)*5.*esatitab(tlonr)+(ttry-tlo)*5.*esatitab(thinr)
                qsatur = ilratio*(rd/rv)*esl1/(presh(k)-(1.-rd/rv)*esl1)+(1.-ilratio)*(rd/rv)*esi1/(presh(k)-(1.-rd/rv)*esi1)
                thlguessmin = ttry/exnh(k)-(rlv/(cp*exnh(k)))*max(qt0h(i,j,k)-qsatur,0.)

                Tnr = Tnr - (thlguess-thl0h(i,j,k))/((thlguess-thlguessmin)*500.)
              enddo
              nitert =max(nitert,niter)
              ilratio = max(0.,min(1.,(Tnr-tdn)/(tup-tdn)))
              tlonr=int((Tnr-150.)*5.)
              thinr=tlonr+1
              tlo=ttab(tlonr)
              thi=ttab(thinr)
              esl1=(thi-Tnr)*5.*esatltab(tlonr)+(Tnr-tlo)*5.*esatltab(thinr)
              esi1=(thi-Tnr)*5.*esatitab(tlonr)+(Tnr-tlo)*5.*esatitab(thinr)
              qvsl1=rd/rv*esl1/(presh(k)-(1.-rd/rv)*esl1)
              qvsi1=rd/rv*esi1/(presh(k)-(1.-rd/rv)*esi1)
              qsatur = ilratio*qvsl1+(1.-ilratio)*qvsi1
            endif
            ql0h(i,j,k) = max(qt0h(i,j,k)-qsatur,0.)
      end do
      end do
      end do
      if(nitert>99) then
        write(*,*) 'thermowarning'
      endif

  end subroutine icethermoh

!> Calculates the scalars at half levels.
!! If the kappa advection scheme is active, interpolation needs to be done consistently.
  subroutine calc_halflev
    use modglobal, only : i1, j1, k1, dzf, dzhi, iadv_thl, iadv_qt, iadv_kappa
    use modfields, only : thl0, thl0h, qt0, qt0h
    use modsurfdata,only: qts, thls
    use advec_kappa,only: halflev_kappa
    implicit none

    integer :: i, j, k

    call timer_tic('modthermodynamics/calc_halflev', 1)

    if (iadv_thl==iadv_kappa) then
      call halflev_kappa(thl0,thl0h)
    else
      !$acc parallel loop collapse(3) default(present) async(1)
      do k = 2, k1
        do j = 2 ,j1
          do i = 2 ,i1
            thl0h(i,j,k) = (thl0(i,j,k)*dzf(k-1)+thl0(i,j,k-1)*dzf(k)) * (0.5_field_r * dzhi(k))
          end do
        end do
      end do
    end if

    !$acc parallel loop collapse(2) default(present) async(1)
    do j = 2, j1
      do i = 2, i1
        thl0h(i,j,1) = thls
      end do
    end do

    if (iadv_qt==iadv_kappa) then
      call halflev_kappa(qt0,qt0h)
    else
      !$acc parallel loop collapse(3) default(present) async(1)
      do k = 2, k1
        do j = 2, j1
          do i = 2, i1
            qt0h(i,j,k)  = (qt0(i,j,k)*dzf(k-1)+qt0(i,j,k-1)*dzf(k)) * (0.5_field_r * dzhi(k))
          end do
        end do
      end do

      !$acc parallel loop collapse(2) default(present) async(1)
      do j = 2, j1
        do i = 2, i1
          qt0h(i,j,1) = qts
        end do
      end do
    end if

    !$acc wait(1)
    call timer_toc('modthermodynamics/calc_halflev')
  end subroutine calc_halflev

end module modthermodynamics
