!> \file modspraying.f99
!! Stephan de Roode and Annelot Broerze



module modspraying
   use modsprayingdata, only: i_glob_spray,j_glob_spray,k_glob_spray,&
                              i_loc_spray,j_loc_spray,k_loc_spray,&
                              water_spray_rate,salt_spray_rate,&
                              lwater_spraying,lsalt_spraying,salinity

contains
  subroutine initspraying
  use modglobal,    only : i1,j1,imax,jmax,kmax,ifnamopt,fname_options,nsv
  use modmpi,       only : myid,myidx,myidy,comm3d, mpierr, d_mpi_bcast
  use modmicrodata, only: imicro

  integer i,j,k,iglob,jglob,ierr

  namelist/NAMSPRAYING/ lwater_spraying,lsalt_spraying,&
                        i_glob_spray,j_glob_spray,k_glob_spray,&
                        water_spray_rate,salt_spray_rate,salinity

  if(myid==0) then    !first myid
    open(ifnamopt,file=fname_options,status='old',iostat=ierr)
    read (ifnamopt,NAMSPRAYING,iostat=ierr)
    if (ierr > 0) then
      print *, 'Problem in namoptions NAMSPRAYING'
      print *, 'iostat error: ', ierr
      stop 'ERROR: Problem in namoptions NAMSPRAYING'
    endif
    write(6 ,NAMSPRAYING)
    close(ifnamopt)
  endif

  call D_MPI_BCAST(lwater_spraying ,    1,  0, comm3d, mpierr)
  call D_MPI_BCAST(lsalt_spraying  ,    1,  0, comm3d, mpierr)
  call D_MPI_BCAST(i_glob_spray    ,    1,  0, comm3d, mpierr)
  call D_MPI_BCAST(j_glob_spray    ,    1,  0, comm3d, mpierr)
  call D_MPI_BCAST(k_glob_spray    ,    1,  0, comm3d, mpierr)
  call D_MPI_BCAST(water_spray_rate,    1,  0, comm3d, mpierr)
  call D_MPI_BCAST(salt_spray_rate,    1,  0, comm3d, mpierr) 

  if (lwater_spraying) then
     lsalt_spraying  = .true. 
     salt_spray_rate = water_spray_rate * salinity   !directy coupled to water spray rate
  else
     water_spray_rate = 0.
  endif

  if (lsalt_spraying) then
     isv_salt = 1
     if (imicro.ne.0) then 
        isv_salt = 3
     endif
     if (isv_salt.gt.nsv) then
         print *, 'Problem with scalar field for salt'
         print *, 'increase number of passive scalar fields'
         stop 'ERROR: not enough passive scalar fields'
     endif
  endif

  !determine local position of spraying from global position

  i_loc_spray = -999
  j_loc_spray = -999
  k_loc_spray = -999
  do i=2,i1
    do j=2,j1
       do k=1,kmax
          iglob = i+myidx*imax
          jglob = j+myidy*jmax
          if (iglob.eq.i_glob_spray.and.jglob.eq.j_glob_spray.and.k.eq.k_glob_spray) then
             i_loc_spray = i
             j_loc_spray = j
             k_loc_spray = k

             write(6,*) 'spraying point at myid = ',myid
             write(6,*) 'global locations ',i_glob_spray,j_glob_spray,k_glob_spray
             write(6,*) 'local locations ',i_loc_spray,j_loc_spray,k_loc_spray
           
          endif 
       enddo
     enddo
  enddo
 
  if (myid==0) then 
     write(6,*) 'Spraying data used: '
     write(6,*) 'lwater_spraying  ',lwater_spraying
     write(6,*) 'lsalt_spraying   ',lsalt_spraying
     write(6,*) 'i_glob_spray     ',i_glob_spray
     write(6,*) 'j_glob_spray     ',j_glob_spray
     write(6,*) 'k_glob_spray     ',k_glob_spray
     write(6,*) 'water_spray_rate ',water_spray_rate
     write(6,*) 'salt_spray_rate  ',salt_spray_rate
     write(6,*)
  endif

  
  end subroutine initspraying

end module modspraying
