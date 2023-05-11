module module_ramses

  use module_constants, only : kB, mp, XH, planck, clight, cmtoA
  use module_domain
  use coolrates_module
  
  implicit none

  private 

  ! Stop pretending this would work in 2D 
  integer(kind=4),parameter :: ndim = 3
  integer(kind=4),parameter :: twondim = 6
  integer(kind=4),parameter :: twotondim= 8 

  ! QuadHilbert related precision
  integer,parameter::qdp=kind(1.0_16) ! real*16


  integer(kind=4)                  :: ncpu
  integer(kind=4)                  :: U_precision=8 ! hydro-precision in RAMSES output
  integer(kind=4)                  :: RT_precision=8 ! RT-precision in RAMSES output

  ! conversion factors (units)
  logical                        :: conversion_scales_are_known = .False. 
  real(kind=8)                   :: dp_scale_l,dp_scale_d,dp_scale_t
  real(kind=8)                   :: dp_scale_T2,dp_scale_zsun,dp_scale_nH
  real(kind=8)                   :: dp_scale_nHe,dp_scale_v,dp_scale_m

  ! cooling-related stuff -------------------------------------------------------------
  type cooling_table
     integer(kind=4)          :: n11
     integer(kind=4)          :: n22
     real(kind=8),allocatable :: nH(:)    
     real(kind=8),allocatable :: T2(:)    
     real(kind=8),allocatable :: metal(:,:)  
     real(kind=8),allocatable :: cool(:,:)  
     real(kind=8),allocatable :: heat(:,:)  
     real(kind=8),allocatable :: cool_com(:,:)  
     real(kind=8),allocatable :: heat_com(:,:)  
     real(kind=8),allocatable :: cool_com_prime(:,:)  
     real(kind=8),allocatable :: heat_com_prime(:,:)  
     real(kind=8),allocatable :: metal_prime(:,:)  
     real(kind=8),allocatable :: cool_prime(:,:)  
     real(kind=8),allocatable :: heat_prime(:,:)  
     real(kind=8),allocatable :: mu(:,:)  
     real(kind=8),allocatable :: spec(:,:,:)  ! last dimension (6) is n_e, n_HI, n_HII, n_HeI, n_HeII, n_HeIII
  end type cooling_table
  type(cooling_table) :: cooling
  type cool_interp
     integer(kind=4)  :: n_nH
     real(kind=8)     :: nH_start,nH_step
     integer(kind=4)  :: n_T2
     real(kind=8)     :: T2_start,T2_step
  end type cool_interp
  type(cool_interp)  :: cool_int
  logical            :: cooling_is_read = .False. 
  ! ----------------------------------------------------------------------------------

  ! particle-related stuff -----------------------------------------------------------
  character(30) :: ParticleFields(20)  ! array of particle fields (e.g. (/'pos','vel','mass','iord','level'/) for a DM-only run)
  ! conformal time things
  integer(kind=4),parameter             :: n_frw = 1000
  real(KIND=8),dimension(:),allocatable :: aexp_frw,hexp_frw,tau_frw,t_frw
  ! ----------------------------------------------------------------------------------


  ! --------------------------------------------------------------------------
  ! user-defined parameters - read from section [ramses] of the parameter file
  ! --------------------------------------------------------------------------
  ! ramses options (not guessable from outputs)
  logical                  :: self_shielding    = .true.   ! if true, reproduce self-shielding approx made in ramses to compute nHI. 
  logical                  :: ramses_rt         = .false.  ! if true, read ramses-RT output and compute nHI and T accordingly.
  logical                  :: read_rt_variables = .false.  ! if true, read RT variables (e.g. to compute heating terms)
  logical                  :: use_initial_mass  = .false.  ! if true, use initial masses of star particles instead of mass at output time
  logical                  :: cosmo             = .true.   ! if false, assume idealised simulation
  logical                  :: use_proper_time   = .false.  ! if true, use proper time instead of conformal time for cosmo runs. 
  ! miscelaneous
  logical                  :: verbose        = .false. ! display some run-time info on this module
  ! RT variable indices
  integer(kind=4) :: itemp  = 5 ! index of thermal pressure
  integer(kind=4) :: imetal = 6 ! index of metallicity 
  integer(kind=4) :: ihii   = 7 ! index of HII fraction 
  integer(kind=4) :: iheii  = 8 ! index of HeII fraction 
  integer(kind=4) :: iheiii = 9 ! index of HeIII fraction 
  ! Solar abundance ratio
  ! Si
  ! abundance_Si_mass == abundance_Si_number * 28.085
  real(kind=8),parameter    :: abundance_Si_number = 3.24d-5 ! From Scarlata (private comm.)
  ! Mg
  ! abundance_Mg_mass == abundance_Mg_number * 24.305
  real(kind=8),parameter    :: abundance_Mg_number = 3.39d-5 ! From Scarlata (private comm.)
  ! Fe
  ! abundance_Fe_mass == abundance_Fe_number * 55.845
  real(kind=8),parameter    :: abundance_Fe_number = 2.82d-5 ! From Scarlata (private comm.)
  ! --------------------------------------------------------------------------
  
  public :: ramses_get_leaf_cells_slomp, ramses_get_leaf_cells
  public :: ramses_get_leaf_cells_in_domain_slomp, get_ngridtot_cpus
  public :: ramses_get_box_size_cm, get_cpu_list, get_cpu_list_periodic, get_ncpu
  public :: ramses_get_velocity_cgs, ramses_get_T_nhi_cgs, ramses_get_metallicity,  ramses_get_nh_cgs
  public :: ramses_get_T_nSiII_cgs, ramses_get_T_nMgII_cgs, ramses_get_T_nFeII_cgs
  public :: ramses_get_nh_nhi_nhei_nehii_cgs
  public :: ramses_read_stars_in_domain
  public :: read_ramses_params, print_ramses_params, dump_ramses_info
  public :: ramses_get_LyaEmiss_HIDopwidth,ramses_get_cooling_time
  
  !==================================================================================
contains

  ! ----------------
  ! public functions 
  ! ----------------


  subroutine ramses_get_leaf_cells(repository, snapnum, ncpu_read, cpu_list, &
       & nleaftot, nvar, xleaf_all, ramses_var_all, leaf_level_all, selection_domain)
    ! non-openMP method, as in minirats...
    ! no subroutine, store cell data directly into final arrays
    
    implicit none 
    character(2000),intent(in)                :: repository
    integer(kind=4),intent(in)                :: snapnum, ncpu_read
    integer(kind=4),allocatable,intent(in)    :: cpu_list(:)    
    integer(kind=4),intent(inout)             :: nleaftot, nvar
    real(kind=8),allocatable, intent(inout)   :: ramses_var_all(:,:)
    real(kind=8),allocatable,intent(inout)    :: xleaf_all(:,:)
    integer(kind=4),allocatable,intent(inout) :: leaf_level_all(:)
    type(domain),intent(in),optional          :: selection_domain
    integer(kind=4)                           :: ileaf,nleaf,k,icpu,ivar,iloop
    real(kind=8)                              :: time1,time2,time3,rate
    integer(kind=8)                           :: c1,c2,c3,cr
    character(1000)                           :: filename 
    logical                                   :: ok_cell
    integer(kind=4)                           :: i,j,ilevel,nx,ny,nz,nlevelmax,nboundary
    integer(kind=4)                           :: idim,ind,iu1,iu2,iu3,rank
    ! stuff read from AMR files
    integer(kind=4)                           :: ngridmax,ngrid_current
    real(kind=8),allocatable                  :: xg(:,:)        ! grids position
    integer(kind=4),allocatable               :: son(:,:)       ! sons grids
    real(KIND=8),dimension(1:3)               :: xbound=(/0d0,0d0,0d0/)  
    integer(kind=4),allocatable               :: ngridfile(:,:),ngridlevel(:,:),ngridbound(:,:)
    integer(kind=4)                           :: ngrida,ncpused
    logical,allocatable                       :: ref(:,:)
    real(kind=8)                              :: dx,boxlen
    integer(kind=4)                           :: ix,iy,iz,nvarH,nvarRT
    real(kind=8),allocatable                  :: xc(:,:),xp(:,:,:)
    ! stuff read from the HYDRO files
    real(kind=8),allocatable                  :: var(:,:,:)
    real(kind=4),allocatable                  :: var_sp(:)
    logical                                   :: cellInDomain
    real(kind=8),dimension(3)                 :: xx
    logical,allocatable                       :: cpu_is_useful(:)
    
    if(verbose) print *,'Reading RAMSES cells...'

    call cpu_time(time1)
    call system_clock(count_rate=cr)
    rate = float(cr)
    call system_clock(c1)

    allocate(cpu_is_useful(ncpu_read))
    cpu_is_useful = .false.

    if(present(selection_domain))then
       call ramses_count_leaf_cells_in_domain(repository, snapnum, ncpu_read, cpu_list, &
            & selection_domain, nleaftot, cpu_is_useful)
       !print*,'nleaftot (new) in selection_domain =',nleaftot
       ncpused=0
       do k = 1,ncpu_read
          if (cpu_is_useful(k)) ncpused=ncpused+1
       end do
       !print*,'--> ncpu to really read : ',ncpused
    else
       nleaftot = get_nleaf_new(repository,snapnum,ncpu_read,cpu_list)
       !print*,'nleaftot (new) =',nleaftot
       cpu_is_useful = .true.
    endif
    
    call cpu_time(time2)
    call system_clock(c2)
   !  print '(" --> Time to get nleaf new = ",f12.3," seconds.")',time2-time1
   !  print '("         system_clock time = ",f12.3," seconds.")',(c2-c1)/rate
    
    nvar     = get_nvar(repository,snapnum)
    allocate(ramses_var_all(nvar,nleaftot), xleaf_all(nleaftot,3), leaf_level_all(nleaftot))

   ! Check whether the ramses output is in single or double precision
    U_precision = nint(get_param_real(repository,snapnum,'U_precision',default_value=8d0))
    if(read_rt_variables) then
       RT_precision = nint(get_param_real(repository,snapnum,'rtprecision' &
            ,default_value=8d0,rt_info=.true.))
       print*,'The RT precision is ',RT_precision  !JOKI
    endif

    if(verbose) print *,'-- ramses_get_leaf_cells : nleaftot(_read), nvar, ncpu(_read) =',nleaftot,nvar,ncpu_read

    
    rank = 1
    iu1 = 10+rank*3
    iu2 = 10+rank*3+1
    iu3 = 10+rank*3+2

    nleaf=0
    ileaf=1
    iloop=0
    ! loop over cpu
    do k=1,ncpu_read
       icpu=cpu_list(k)
#ifdef DISPLAY_PROGRESS_PERCENT
      !  write (*, "(A, f5.2, A, A, $)") &           ! Progress bar that works with ifort
      !       ' Reading leaves ',dble(iloop) / ncpu_read * 100,' % ',char(13)
       iloop=iloop+1
#endif
       if (.not. cpu_is_useful(k)) cycle
       ! verify AMR input file -> already done above in get_nleaf_new
       write(filename,'(a,a,i5.5,a,i5.5,a,i5.5)') trim(repository),'/output_',snapnum,'/amr_',snapnum,'.out',icpu
       ! Open AMR file and skip header
       open(unit=iu1,file=filename,form='unformatted',status='old',action='read')
       read(iu1)ncpu
       read(iu1)      !ndim
       read(iu1)nx,ny,nz
       read(iu1)nlevelmax
       read(iu1)ngridmax
       read(iu1)nboundary
       read(iu1)ngrid_current
       read(iu1)boxlen
       do i=1,13
          read(iu1)
       end do
       !twotondim=2**ndim
       xbound=(/dble(nx/2),dble(ny/2),dble(nz/2)/)
       if(allocated(ngridfile)) deallocate(ngridfile,ngridlevel)
       allocate(ngridfile(1:ncpu+nboundary,1:nlevelmax))
       allocate(ngridlevel(1:ncpu,1:nlevelmax))
       if(nboundary>0)then
          if(allocated(ngridbound)) deallocate(ngridbound)
          allocate(ngridbound(1:nboundary,1:nlevelmax))
       endif
       ! Read grid numbers
       read(iu1)ngridlevel
       ngridfile(1:ncpu,1:nlevelmax)=ngridlevel
       read(iu1)
       if(nboundary>0)then
          do i=1,2
             read(iu1)
          end do
          read(iu1)ngridbound
          ngridfile(ncpu+1:ncpu+nboundary,1:nlevelmax)=ngridbound
       endif
       read(iu1)
       ! ROM: comment the single follwing line for old stuff
       read(iu1)
       read(iu1)
       read(iu1)
       read(iu1)
       read(iu1)
       
       if(allocated(xc)) deallocate(xc)
       allocate(xc(1:twotondim,1:ndim))
       
       
       ! open hydro file and get nvarH
       write(filename,'(a,a,i5.5,a,i5.5,a,i5.5)') trim(repository),'/output_',snapnum,'/hydro_',snapnum,'.out',icpu
       open(unit=iu2,file=filename,form='unformatted',status='old',action='read')
       read(iu2)
       read(iu2)nvarH
       read(iu2)
       read(iu2)
       read(iu2)
       read(iu2)
       
       if (read_rt_variables) then
          ! Open RT file and get nvarRT
          write(filename,'(a,a,i5.5,a,i5.5,a,i5.5)') trim(repository),'/output_',snapnum,'/rt_',snapnum,'.out',icpu
          open(unit=iu3,file=filename,status='old',form='unformatted')
          read(iu3)
          read(iu3)nvarRT
          read(iu3)
          read(iu3)
          read(iu3)
          read(iu3)
       else
          nvarRT = 0
       end if
       !ncoarse = nx*ny*nz
       !ncell   = ncoarse+twotondim*ngridmax
       
       ! Loop over levels
       do ilevel=1,nlevelmax
          
          ! Geometry
          dx=0.5**ilevel
          do ind=1,twotondim
             iz=(ind-1)/4
             iy=(ind-1-4*iz)/2
             ix=(ind-1-2*iy-4*iz)
             xc(ind,1)=(dble(ix)-0.5D0)*dx
             xc(ind,2)=(dble(iy)-0.5D0)*dx
             xc(ind,3)=(dble(iz)-0.5D0)*dx
          end do
          
          ! Allocate work arrays
          if(allocated(xg)) then 
             deallocate(xg,son,var,xp,ref)
          endif
          if(allocated(var_sp)) deallocate(var_sp)
          ngrida=ngridfile(icpu,ilevel)
          if(ngrida>0)then
             allocate(xg(1:ngrida,1:ndim))
             allocate(son(1:ngrida,1:twotondim))
             allocate(var(1:ngrida,1:twotondim,1:nvarh+nvarRT))
             if((read_rt_variables .and. rt_Precision.eq.4) .or. U_precision.eq.4) allocate(var_sp(1:ngrida))
             allocate(xp(1:ngrida,1:twotondim,1:ndim))
             allocate(ref(1:ngrida,1:twotondim))
             ref=.false.
          endif
          
          
          ! Loop over domains
          do j=1,nboundary+ncpu
             
             ! Read AMR data
             if(ngridfile(j,ilevel)>0)then
                read(iu1) ! Skip grid index
                read(iu1) ! Skip next index
                read(iu1) ! Skip prev index
                ! Read grid center
                do idim=1,ndim
                   if(j.eq.icpu)then
                      read(iu1)xg(:,idim)
                   else
                      read(iu1)
                   endif
                end do
                read(iu1) ! Skip father index
                do ind=1,2*ndim
                   read(iu1) ! Skip nbor index
                end do
                ! Read son index
                do ind=1,twotondim
                   if(j.eq.icpu)then
                      read(iu1)son(:,ind)
                   else
                      read(iu1)
                   end if
                end do
                ! Skip cpu map
                do ind=1,twotondim
                   read(iu1)
                end do
                ! Skip refinement map
                do ind=1,twotondim
                   read(iu1)
                end do
             endif
             
             ! Read HYDRO data
             read(iu2)
             read(iu2)
             if(read_rt_variables)read(iu3)
             if(read_rt_variables)read(iu3)
             if(ngridfile(j,ilevel)>0)then
                ! Read hydro variables
                do ind=1,twotondim
                   do ivar=1,nvarh
                      if(j.eq.icpu)then
                         if(U_precision.eq.4) then
                            read(iu2) var_sp(:)
                            var(:,ind,ivar) = var_sp(:)
                         else
                            read(iu2)var(:,ind,ivar)
                         endif
                      else
                         read(iu2)
                      end if
                   end do
                   do ivar=1,nvarRT
                      if(j.eq.icpu)then
                         if(rt_Precision.eq.4) then
                            read(iu3) var_sp(:)
                            var(:,ind,nvarh+ivar) = var_sp(:)
                         else
                            read(iu3)var(:,ind,nvarh+ivar)
                         endif
                      else
                         read(iu3)
                      end if
                   end do
                end do
             end if
             
          enddo ! end loop over domains
          
          ! Get leaf cells and store data
          if(ngrida>0)then
             ! Loop over cells
             do ind=1,twotondim
                ! Compute cell center
                do i=1,ngrida
                   xp(i,ind,1)=(xg(i,1)+xc(ind,1)-xbound(1))
                   xp(i,ind,2)=(xg(i,2)+xc(ind,2)-xbound(2))
                   xp(i,ind,3)=(xg(i,3)+xc(ind,3)-xbound(3))
                end do
                ! Check if cell is refined
                do i=1,ngrida
                   ref(i,ind)=son(i,ind)>0.and.ilevel<nlevelmax
                end do
                ! Store leaf cells
                do i=1,ngrida
                   ok_cell= .not.ref(i,ind)
                   if(ok_cell)then
                   !if(.not.ref(i,ind))then
                      cellInDomain=.true.
                      if(present(selection_domain))then
                         !
                         cellInDomain=.false.
                         xx(1:3) = xp(i,ind,1:3)
                         dx = 0.5d0**(ilevel)
                         if (domain_contains_cell(xx,dx,selection_domain)) cellInDomain=.true.
                      endif
                      if(cellInDomain)then
                         xleaf_all(ileaf,1:3) = xp(i,ind,1:3)
                         leaf_level_all(ileaf) = ilevel
                         do ivar = 1,nvar
                            ramses_var_all(ivar,ileaf) = var(i,ind,ivar)
                         end do
                         ileaf=ileaf+1
                      endif
                   endif
                enddo
             end do
          endif
          
          
       enddo ! end loop over levels
       
       
       close(iu1)
       close(iu2)
       close(iu3)
       
    enddo ! end loop over cpu

    !!print*,'Number of leaf cells =',ileaf-1,nvarH+nvarRT
    !!print*,icpu,ileaf-1
    nleaf = ileaf-1

    call cpu_time(time3)
    call system_clock(c3)
   !  print '(" --> Time to get leaf = ",f12.3," seconds.")',time3-time2
   !  print '("    system_clock time = ",f12.3," seconds.")',(c3-c2)/rate

   !  print*,'Nleaf read = ',nleaf, nleaftot
    return
  end subroutine ramses_get_leaf_cells



  subroutine ramses_count_leaf_cells_in_domain(repository, snapnum, &
       & ncpu_read, cpu_list, selection_domain, nleafInDomain, cpu_is_useful)
    ! count leaf cells in cpu_list && in selection_domain
    
    implicit none 
    character(2000),intent(in)              :: repository
    integer(kind=4),intent(in)              :: snapnum, ncpu_read
    integer(kind=4),allocatable,intent(in)  :: cpu_list(:)
    type(domain),intent(in)                 :: selection_domain
    integer(kind=4),intent(inout)           :: nleafInDomain
    logical,allocatable,intent(inout)       :: cpu_is_useful(:)
    integer(kind=4)                         :: k,icpu,iloop
    character(1000)                         :: filename 
    logical                                 :: ok_cell
    integer(kind=4)                         :: i,j,ilevel,nx,ny,nz,nlevelmax,nboundary
    integer(kind=4)                         :: idim,ind,iu1,iu2,rank
    integer(kind=4)                         :: ngridmax,ngrid_current
    real(kind=8),allocatable                :: xg(:,:)        ! grids position
    integer(kind=4),allocatable             :: son(:,:)       ! sons grids
    real(KIND=8),dimension(1:3)             :: xbound=(/0d0,0d0,0d0/),xx  
    integer(kind=4),allocatable             :: ngridfile(:,:),ngridlevel(:,:),ngridbound(:,:)
    integer(kind=4)                         :: ngrida,nleafInCpu
    logical,allocatable                     :: ref(:,:)
    real(kind=8)                            :: dx,boxlen
    integer(kind=4)                         :: ix,iy,iz,nvarH
    real(kind=8),allocatable                :: xc(:,:),xp(:,:,:)
    
    
    rank = 1
    iu1 = 10+rank*3
    iu2 = 10+rank*3+1

    nleafInDomain=0
    iloop=0
    ! loop over cpu
    do k=1,ncpu_read
       icpu=cpu_list(k)
#ifdef DISPLAY_PROGRESS_PERCENT
       write (*, "(A, f5.2, A, A, $)") &           ! Progress bar that works with ifort
            ' Counting leaves ',dble(iloop) / ncpu_read * 100,' % ',char(13)
       iloop=iloop+1
#endif
       
       ! verify AMR input file -> already done above in get_nleaf_new
       write(filename,'(a,a,i5.5,a,i5.5,a,i5.5)') trim(repository),'/output_',snapnum,'/amr_',snapnum,'.out',icpu
       ! Open AMR file and skip header
       open(unit=iu1,file=filename,form='unformatted',status='old',action='read')
       read(iu1)ncpu
       read(iu1)      !ndim
       read(iu1)nx,ny,nz
       read(iu1)nlevelmax
       read(iu1)ngridmax
       read(iu1)nboundary
       read(iu1)ngrid_current
       read(iu1)boxlen
       do i=1,13
          read(iu1)
       end do
       !twotondim=2**ndim
       xbound=(/dble(nx/2),dble(ny/2),dble(nz/2)/)
       if(allocated(ngridfile)) deallocate(ngridfile,ngridlevel)
       allocate(ngridfile(1:ncpu+nboundary,1:nlevelmax))
       allocate(ngridlevel(1:ncpu,1:nlevelmax))
       if(nboundary>0)then
          if(allocated(ngridbound)) deallocate(ngridbound)
          allocate(ngridbound(1:nboundary,1:nlevelmax))
       endif
       ! Read grid numbers
       read(iu1)ngridlevel
       ngridfile(1:ncpu,1:nlevelmax)=ngridlevel
       read(iu1)
       if(nboundary>0)then
          do i=1,2
             read(iu1)
          end do
          read(iu1)ngridbound
          ngridfile(ncpu+1:ncpu+nboundary,1:nlevelmax)=ngridbound
       endif
       read(iu1)
       ! ROM: comment the single follwing line for old stuff
       read(iu1)
       read(iu1)
       read(iu1)
       read(iu1)
       read(iu1)
       
       if(allocated(xc)) deallocate(xc)
       allocate(xc(1:twotondim,1:ndim))
       
       
       ! open hydro file and get nvarH
       write(filename,'(a,a,i5.5,a,i5.5,a,i5.5)') trim(repository),'/output_',snapnum,'/hydro_',snapnum,'.out',icpu
       open(unit=iu2,file=filename,form='unformatted',status='old',action='read')
       read(iu2)
       read(iu2)nvarH
       read(iu2)
       read(iu2)
       read(iu2)
       read(iu2)

       nleafInCpu=0

       ! Loop over levels
       do ilevel=1,nlevelmax
          
          ! Geometry
          dx=0.5**ilevel
          do ind=1,twotondim
             iz=(ind-1)/4
             iy=(ind-1-4*iz)/2
             ix=(ind-1-2*iy-4*iz)
             xc(ind,1)=(dble(ix)-0.5D0)*dx
             xc(ind,2)=(dble(iy)-0.5D0)*dx
             xc(ind,3)=(dble(iz)-0.5D0)*dx
          end do
          
          ! Allocate work arrays
          if(allocated(xg)) then 
             deallocate(xg,son,xp,ref)
          endif
          ngrida=ngridfile(icpu,ilevel)
          if(ngrida>0)then
             allocate(xg(1:ngrida,1:ndim))
             allocate(son(1:ngrida,1:twotondim))
             allocate(xp(1:ngrida,1:twotondim,1:ndim))
             allocate(ref(1:ngrida,1:twotondim))
             ref=.false.
          endif
          
          
          ! Loop over domains
          do j=1,nboundary+ncpu
             
             ! Read AMR data
             if(ngridfile(j,ilevel)>0)then
                read(iu1) ! Skip grid index
                read(iu1) ! Skip next index
                read(iu1) ! Skip prev index
                ! Read grid center
                do idim=1,ndim
                   if(j.eq.icpu)then
                      read(iu1)xg(:,idim)
                   else
                      read(iu1)
                   endif
                end do
                read(iu1) ! Skip father index
                do ind=1,2*ndim
                   read(iu1) ! Skip nbor index
                end do
                ! Read son index
                do ind=1,twotondim
                   if(j.eq.icpu)then
                      read(iu1)son(:,ind)
                   else
                      read(iu1)
                   end if
                end do
                ! Skip cpu map
                do ind=1,twotondim
                   read(iu1)
                end do
                ! Skip refinement map
                do ind=1,twotondim
                   read(iu1)
                end do
             endif
             ! no need to read HYDRO data...
          enddo ! end loop over domains
          
          ! Get leaf cells and store data
          if(ngrida>0)then
             ! Loop over cells
             do ind=1,twotondim
                ! Compute cell center
                do i=1,ngrida
                   xp(i,ind,1)=(xg(i,1)+xc(ind,1)-xbound(1))
                   xp(i,ind,2)=(xg(i,2)+xc(ind,2)-xbound(2))
                   xp(i,ind,3)=(xg(i,3)+xc(ind,3)-xbound(3))
                end do
                ! Check if cell is refined
                do i=1,ngrida
                   ref(i,ind)=son(i,ind)>0.and.ilevel<nlevelmax
                end do
                ! check if cell is in selection_domain
                do i=1,ngrida
                   ok_cell= .not.ref(i,ind)
                   if(ok_cell)then
                   !if(.not.ref(i,ind))then
                      xx(1:3) = xp(i,ind,1:3)
                      dx = 0.5d0**(ilevel)
                      if (domain_contains_cell(xx,dx,selection_domain)) then
                         nleafInDomain=nleafInDomain+1
                         nleafInCpu=nleafInCpu+1
                      endif
                   endif
                enddo
             end do
          endif
       enddo ! end loop over levels
       
       close(iu1)
       close(iu2)
       cpu_is_useful(k) = (nLeafInCpu > 0)
       
    enddo ! end loop over cpu
    
    return
  end subroutine ramses_count_leaf_cells_in_domain
  
  
  
  subroutine ramses_get_leaf_cells_slomp(repository, snapnum, ncpu_read, cpu_list, &
       & nleaftot, nvar, xleaf_all, ramses_var_all, leaf_level_all)
    
    ! read all leaf cells from a simulation snapshot belonging to given
    ! list of cpus. Return standard ramses variables through 
    ! ramses_var(nvar,nleaftot) and positions (xleaf(3,nleaftot))
    ! and levels (leaf_level).
    ! 05-2020: corrected version that doesn't use cpu_map anymore.
    !          --> not very efficient
    
    !$ use OMP_LIB
    implicit none 
    character(2000),intent(in)                :: repository
    integer(kind=4),intent(in)                :: snapnum, ncpu_read
    integer(kind=4),dimension(:),allocatable,intent(in) :: cpu_list

    integer(kind=4),intent(inout)             :: nleaftot, nvar
    real(kind=8),allocatable, intent(inout)   :: ramses_var_all(:,:)
    real(kind=8),allocatable,intent(inout)    :: xleaf_all(:,:)
    integer(kind=4),allocatable,intent(inout) :: leaf_level_all(:)

    real(kind=8),allocatable     :: cell_pos(:,:),cell_var(:,:)
    integer(kind=4),allocatable  :: cell_lev(:)
    
    integer(kind=4) :: ileaf,nleaf,k,icpu,ilast,ncell,ivar,iloop
    real(kind=8) :: time1,time2,time3,rate,ot1,ot2,ot3
    integer(kind=8) :: c1,c2,c3,cr
    
    if(verbose) print *,'Reading RAMSES cells...'

    ot1=0.
    ot2=0.
    ot3=0.
    call system_clock(count_rate=cr)
    rate = float(cr)
    call system_clock(c1)
    call cpu_time(time1)
    !$ ot1 = omp_get_wtime()
    
    nleaftot = get_nleaf_new(repository,snapnum,ncpu_read,cpu_list)
    !print*,'nleaftot (new) =',nleaftot

    call cpu_time(time2)
    call system_clock(c2)
    !$ ot2 = omp_get_wtime()
   !  print '(" --> Time to get nleaf new = ",f12.3," seconds.")',time2-time1
   !  print '("         system_clock time = ",f12.3," seconds.")',(c2-c1)/rate
   !  print '("             omp_get_wtime = ",f12.3," seconds.")',ot2-ot1
    
    nvar     = get_nvar(repository,snapnum)
    allocate(ramses_var_all(nvar,nleaftot), xleaf_all(nleaftot,3), leaf_level_all(nleaftot))

   ! Check whether the ramses output is in single or double precision
    U_precision = nint(get_param_real(repository,snapnum,'U_precision',default_value=8d0))
    if(read_rt_variables) then
       RT_precision = nint(get_param_real(repository,snapnum,'rtprecision' &
            ,default_value=8d0,rt_info=.true.))
       print*,'The RT precision is ',RT_precision  !JOKI
    endif

    if(verbose) print *,'-- ramses_get_leaf_cells : nleaftot(_read), nvar, ncpu(_read) =',nleaftot,nvar,ncpu_read

    
    nleaf=0
    ilast=1
    iloop=0
!$OMP PARALLEL &
!$OMP DEFAULT(private) &
!$OMP SHARED(nleaf, repository, snapnum, ncpu_read, cpu_list, ilast, nvar, xleaf_all, leaf_level_all, ramses_var_all, iloop)
!!!!!$OMP DO
!$OMP DO SCHEDULE(DYNAMIC, 5) 
    do k=1,ncpu_read
       icpu=cpu_list(k)

       call get_leaf_cells_per_cpu(repository,snapnum,icpu,ileaf,ncell,cell_pos,cell_var,cell_lev)

!$OMP CRITICAL
#ifdef DISPLAY_PROGRESS_PERCENT
       write (*, "(A, f5.2, A, A, $)") &           ! Progress bar that works with ifort
            ' Reading leaves ',dble(iloop) / ncpu_read * 100,' % ',char(13)
       iloop=iloop+1
#endif
       ! ileaf is now the number of leaves on local cpu
       if(ileaf .gt. 0) then
          ! save leaf cells to return arrays
          xleaf_all(ilast:ilast-1+ileaf,1:3)  = cell_pos(1:ileaf,1:3)
          leaf_level_all(ilast:ilast-1+ileaf) = cell_lev(1:ileaf)
          !do ivar = 1,nvar
          !   ramses_var_all(ivar,ilast:ilast-1+ileaf) = cell_var(1:ileaf,ivar)
          !end do
          ramses_var_all(1:nvar,ilast:ilast-1+ileaf) = cell_var(1:nvar,1:ileaf)
       endif
       ilast=ilast+ileaf
       nleaf=nleaf+ileaf
!$OMP END CRITICAL
    end do
!$OMP END DO
!$OMP END PARALLEL

    call cpu_time(time3)
    print '(" --> Time to get leaf = ",f12.3," seconds.")',time3-time2
    call system_clock(c3)
    print '("    system_clock time = ",f12.3," seconds.")',(c3-c2)/rate
    !$ ot3 = omp_get_wtime()
    print '("        omp_get_wtime = ",f12.3," seconds.")',ot3-ot2

    print*,'Nleaf read = ',nleaf
    return
  end subroutine ramses_get_leaf_cells_slomp



  subroutine ramses_get_leaf_cells_in_domain_slomp(repository, snapnum, selection_domain, &
       & ncpu_read, cpu_list, &
       & nleaftot_all, nvar, xleaf_all, ramses_var_all, leaf_level_all)

    ! Read leaf cells from a selection of cpu files in a simulation snapshot
    ! and select cells inside selection_domain.
    !
    ! Return standard ramses variables through ramses_var(nvar,nleaftot) and
    ! positions (xleaf(3,nleaftot)) and levels (leaf_level).
    !
    ! Slower method since it needs to read 2 times each file (one to count, one to extract)
    
    implicit none 

    character(2000),intent(in)                :: repository
    integer(kind=4),intent(in)                :: snapnum, ncpu_read
    type(domain),intent(in)                   :: selection_domain
    integer(kind=4),intent(inout)             :: nleaftot_all, nvar
    integer(kind=4)                           :: nleaftot
    real(kind=8),allocatable                  :: ramses_var(:,:)
    real(kind=8),allocatable                  :: xleaf(:,:)
    integer(kind=4),allocatable               :: leaf_level(:)
    real(kind=8),allocatable, intent(inout)   :: ramses_var_all(:,:)
    real(kind=8),allocatable,intent(inout)    :: xleaf_all(:,:)
    integer(kind=4),allocatable,intent(inout) :: leaf_level_all(:)
    integer(kind=4),allocatable,intent(in)    :: cpu_list(:)
    
    logical                                   :: do_allocs
    integer(kind=4)                           :: icpu, ivar, nleaf_in_domain, k, i
    integer(kind=4)                           :: ileaf, ilast, iloop=0, nleaf, ncell
    real(kind=8),dimension(3)                 :: xtemp
    real(kind=8)                              :: dx

    real(kind=8),allocatable                  :: cell_var(:,:)
    real(kind=8),allocatable                  :: cell_pos(:,:)
    integer(kind=4),allocatable               :: cell_lev(:)
    integer(kind=4)                           :: nLeafInCpu
    logical,allocatable                       :: cpu_is_useful(:)
    
    if(verbose) print *,'Reading RAMSES cells...'

    nvar = get_nvar(repository,snapnum)
    ncpu = get_ncpu(repository,snapnum)
    ! Check whether the ramses output is in single or double precision
    U_precision = nint(get_param_real(repository,snapnum,'U_precision',default_value=8d0))
    if(read_rt_variables) then
       RT_precision = nint(get_param_real(repository,snapnum,'rtprecision' &
            ,default_value=8d0,rt_info=.true.))
       print*,'The RT precision is ',RT_precision
    endif
    allocate(cpu_is_useful(ncpu_read))
    cpu_is_useful = .false.

    ! first count leaf cells in selection_domain...
    nleaftot = 0 ; nleaf_in_domain = 0 ; iloop=0
!$OMP PARALLEL &
!$OMP REDUCTION(+:nleaftot,nleaf_in_domain) &
!$OMP DEFAULT(private) &
!$OMP SHARED(iloop,repository, snapnum, ncpu_read, cpu_list, selection_domain, cpu_is_useful)
!$OMP DO
    do k=1,ncpu_read
       icpu=cpu_list(k)
       nLeafInCpu = 0

       call get_leaf_cells_per_cpu(repository,snapnum,icpu,ileaf,ncell,cell_pos,cell_var,cell_lev)
       
       ! count leaf cells in selection_domain
       do i = 1,ileaf
          xtemp(:) = cell_pos(i,:)
          dx = 0.5d0**(cell_lev(i))
          nleaftot = nleaftot+1
          if (domain_contains_cell(xtemp,dx,selection_domain)) then
             nleaf_in_domain = nleaf_in_domain + 1
             nLeafInCpu = nLeafInCpu + 1
          end if
       end do
!$OMP CRITICAL
       cpu_is_useful(k) = (nLeafInCpu > 0)
#ifdef DISPLAY_PROGRESS_PERCENT
       !write (*, "(A, f5.2, A, A)", advance='no') &           ! Progress bar
       !     ' Counting leaves ',dble(iloop) / ncpu_read * 100,' % ',char(13)
       write (*, "(A, f5.2, A, A, $)") &           ! Progress bar that works with ifort
            ' Counting leaves ',dble(iloop) / ncpu_read * 100,' % ',char(13)
#endif
       iloop=iloop+1
!$OMP END CRITICAL 
    end do
!$OMP END DO
!$OMP END PARALLEL

    ! JB--
    nLeafInCpu = 0
    do k = 1,ncpu_read
       if (cpu_is_useful(k)) nLeafInCpu = nLeafInCpu + 1
    end do
    !print*,'--> ncpu to really read : ',nLeafInCpu
    ! --JB
    
    if(verbose)print *,'-- ramses_get_leaf_cells_in_domain: nleaftot, nleaf_in_domain, nvar, ncpu =',nleaftot, nleaf_in_domain, nvar, ncpu_read
    
    allocate(ramses_var_all(nvar,nleaf_in_domain), xleaf_all(nleaf_in_domain,3), leaf_level_all(nleaf_in_domain))
    ilast = 1
    iloop=0
!$OMP PARALLEL &
!$OMP DEFAULT(private) &
!$OMP SHARED(iloop,ilast, xleaf_all, leaf_level_all, ramses_var_all, repository, snapnum, nvar, nleaftot, ncpu_read, cpu_list, selection_domain, cpu_is_useful, nLeafInCpu)
    do_allocs=.true.
!$OMP DO
    do k=1,ncpu_read
       ! JB--
       if (.not. cpu_is_useful(k)) cycle
       !--JB
       icpu=cpu_list(k)
       call get_leaf_cells_per_cpu(repository,snapnum,icpu,nleaf,ncell,cell_pos,cell_var,cell_lev)
       if (do_allocs) allocate(ramses_var(nvar,ncell), xleaf(ncell,3), leaf_level(ncell))
       do_allocs=.false.
       ! collect leaf cells
       ileaf=0
       do i = 1,nleaf
          xtemp(:) = cell_pos(i,:)
          dx = 0.5d0**(cell_lev(i))
          if (domain_contains_cell(xtemp,dx,selection_domain)) then
             ileaf = ileaf + 1
             !do ivar = 1,nvar
             !   ramses_var(ivar,ileaf) = cell_var(i,ivar)
             !end do
             ramses_var(1:nvar,ileaf) = cell_var(1:nvar,i)
             xleaf(ileaf,:)    = cell_pos(i,:)
             leaf_level(ileaf) = cell_lev(i)
          end if
       end do
!$OMP CRITICAL
#ifdef DISPLAY_PROGRESS_PERCENT
       !write (*, "(A, f5.2, A, A)", advance='no') &           ! Progress bar
       !     ' Reading leaves ',dble(iloop) / nLeafInCpu * 100,' % ',char(13)
       write (*, "(A, f5.2, A, A, $)") &           ! Progress bar that works with ifort
            ' Reading leaves ',dble(iloop) / nLeafInCpu * 100,' % ',char(13)
#endif
       ! only one CRITICAL zone
       iloop=iloop+1
       ! ileaf is now the number of leaves on local cpu
       if(ileaf .gt. 0) then
          ! save leaf cells to return arrays
          xleaf_all(ilast:ilast-1+ileaf,1:3)  = xleaf(1:ileaf,1:3)
          leaf_level_all(ilast:ilast-1+ileaf) = leaf_level(1:ileaf)
          ramses_var_all(1:nvar,ilast:ilast-1+ileaf) = ramses_var(1:nvar,1:ileaf)
       endif
       ilast=ilast+ileaf
!$OMP END CRITICAL
    end do
!$OMP END DO
    if(allocated(ramses_var)) deallocate(ramses_var,xleaf,leaf_level)
!$OMP END PARALLEL
    
    nleaftot_all = nleaf_in_domain
    
    return

  end subroutine ramses_get_leaf_cells_in_domain_slomp



  subroutine get_cpu_list(repository, snapnum, xmin,xmax,ymin,ymax,zmin,zmax, ncpu_read, cpu_list)

    implicit none

    character(2000),intent(in) :: repository
    integer(kind=4),intent(in) :: snapnum
    real(KIND=8),intent(in) :: xmin,xmax,ymin,ymax,zmin,zmax

    integer(kind=4) :: i,j,ndom
    integer(kind=4) :: ilevel, lmax, lmin
    integer(kind=4) :: ncpu_read, ncpu
    integer(kind=4) :: imin,imax,jmin,jmax,kmin,kmax
    integer(kind=4) :: impi,bit_length,maxdom
    integer(kind=4),dimension(1:8):: idom,jdom,kdom,cpu_min,cpu_max
    real(qdp),dimension(1:8):: bounding_min,bounding_max, order_min
    real(qdp)::dkey
    !!!real(kind=8),dimension(1:8):: bounding_min_dp,bounding_max_dp, order_min_dp
    real(KIND=8)::dmax,dx

    real(qdp),dimension(:),allocatable :: bound_key
    logical,dimension(:),allocatable      :: cpu_read
    integer(kind=4),dimension(:),allocatable,intent(out)      :: cpu_list

    
    lmax = nint(get_param_real(repository,snapnum,'levelmax'))
    ncpu = get_ncpu(repository,snapnum)
    
    dmax = max(xmax-xmin,ymax-ymin,zmax-zmin)

    allocate(cpu_list(1:ncpu))
    allocate(bound_key(0:ncpu))
    allocate(cpu_read(1:ncpu))
    cpu_read=.false.
    cpu_list=0
    
    call read_hilbert_keys_raw(repository,snapnum,ncpu,bound_key)
    
    do ilevel=1,lmax
       dx=0.5d0**ilevel
       if(dx.lt.dmax)exit
    end do
    lmin=ilevel
    bit_length=lmin-1
    maxdom=2**bit_length
    imin=0; imax=0; jmin=0; jmax=0; kmin=0; kmax=0
    if(bit_length>0)then
       imin=int(xmin*dble(maxdom))
       imax=imin+1
       jmin=int(ymin*dble(maxdom))
       jmax=jmin+1
       kmin=int(zmin*dble(maxdom))
       kmax=kmin+1
    endif

    !dkey=(dble(2**(lmax+1)/dble(maxdom)))**ndim
    dkey=(real(2**(lmax+1),kind=qdp)/real(maxdom,kind=qdp))**ndim
    ndom=1
    if(bit_length>0)ndom=8
    idom(1)=imin; idom(2)=imax
    idom(3)=imin; idom(4)=imax
    idom(5)=imin; idom(6)=imax
    idom(7)=imin; idom(8)=imax
    jdom(1)=jmin; jdom(2)=jmin
    jdom(3)=jmax; jdom(4)=jmax
    jdom(5)=jmin; jdom(6)=jmin
    jdom(7)=jmax; jdom(8)=jmax
    kdom(1)=kmin; kdom(2)=kmin
    kdom(3)=kmin; kdom(4)=kmin
    kdom(5)=kmax; kdom(6)=kmax
    kdom(7)=kmax; kdom(8)=kmax

    do i=1,ndom
       if(bit_length>0)then
          call hilbert3d(idom(i),jdom(i),kdom(i),order_min(i),bit_length,1)
       else
          order_min(i)=0.0d0
       endif
       bounding_min(i)=(order_min(i))*dkey
       bounding_max(i)=(order_min(i)+1.0D0)*dkey
    end do

    cpu_min=0; cpu_max=0
    do impi=1,ncpu
       do i=1,ndom
          if (   bound_key(impi-1).le.bounding_min(i).and.&
               & bound_key(impi  ).gt.bounding_min(i))then
             cpu_min(i)=impi
          endif
          if (   bound_key(impi-1).lt.bounding_max(i).and.&
               & bound_key(impi  ).ge.bounding_max(i))then
             cpu_max(i)=impi
          endif
       end do
    end do

    ncpu_read=0
    do i=1,ndom
       do j=cpu_min(i),cpu_max(i)
          if(.not. cpu_read(j))then
             ncpu_read=ncpu_read+1
             cpu_list(ncpu_read)=j
             cpu_read(j)=.true.
          endif
       enddo
    enddo

    deallocate(bound_key,cpu_read)

    !print*,'--> nCPU to read = ',ncpu_read, cpu_list(1:ncpu_reaD)
    
    return
    
  end subroutine get_cpu_list
  
  
  subroutine get_cpu_list_periodic(repository, snapnum, xmin, xmax, ymin, &
       & ymax, zmin, zmax, ncpu_read, cpu_list)
    
    implicit none

    character(2000),intent(in) :: repository
    integer(kind=4),intent(in) :: snapnum
    real(KIND=8),intent(in) :: xmin,xmax,ymin,ymax,zmin,zmax

    ! Can have dom_min<0 and dom_max>1
    real(KIND=8),dimension(3)::     dom_min,dom_max
    ! Can NOT have dom_pmin<0 and dom_pmax>1
    ! (therefore extra domains across boundaires)
    real(KIND=8),dimension(3,1:2):: dom_pmin,dom_pmax
    integer(kind=4) :: i,j,ix,iy,iz,ndom
    integer(kind=4) :: ilevel, lmax, lmin
    integer(kind=4) :: ncpu_read, ncpu
    integer(kind=4) :: imin,imax,jmin,jmax,kmin,kmax
    integer(kind=4) :: impi,bit_length,maxdom
    integer(kind=4),dimension(1:8):: idom,jdom,kdom,cpu_min,cpu_max
    real(qdp),dimension(1:8):: bounding_min,bounding_max, order_min
    real(qdp)::dkey
    real(KIND=8)::dx,dmin,dmax

    real(qdp),dimension(:),allocatable :: bound_key
    logical,dimension(:),allocatable      :: cpu_read
    integer(kind=4),dimension(:),allocatable,intent(out)      :: cpu_list


    lmax = nint(get_param_real(repository,snapnum,'levelmax'))
    ncpu = get_ncpu(repository,snapnum)
    
    allocate(cpu_list(1:ncpu))
    allocate(bound_key(0:ncpu))
    allocate(cpu_read(1:ncpu))
    cpu_read=.false.
    cpu_list=0
    
    if(verbose) write(*,*)'Getting CPU list...'
    
    call read_hilbert_keys_raw(repository,snapnum,ncpu,bound_key)
    
    ! Set up the periodic domains
    dom_min(1) = xmin;   dom_max(1) = xmax
    dom_min(2) = ymin;   dom_max(2) = ymax
    dom_min(3) = zmin;   dom_max(3) = zmax
    dom_pmin(:,1) = dom_min ! Main domain ...
    dom_pmax(:,1) = dom_max ! ...always used
    dom_pmin(:,2) = 1.
    dom_pmax(:,2) = 0.
    do i=1,3 ! loop x,y,z
       if (dom_min(i) < 0.) then ! Lower limit left of boundary
          dom_pmin(i,1) = 0.     ! Adjust main domain
          dom_pmin(i,2) = 1. + dom_min(i) ! Extra domain
          dom_pmax(i,2) = 1.
       endif
       if (dom_max(i) > 1.) then ! Upper limit right of boundary
          dom_pmax(i,1) = 1.     ! Adjust main domain
          dom_pmax(i,2) = dom_min(i) - 1. ! Extra domain
          dom_pmin(i,2) = 0.
       endif
    end do

    ncpu_read=0
    do ix=1,2
       do iy=1,2
          do iz=1,2
             
             dmin = min(dom_pmax(1,ix)-dom_pmin(1,ix)  &
                  ,dom_pmax(2,iy)-dom_pmin(2,iy)  &
                  ,dom_pmax(3,iz)-dom_pmin(3,iz) )
             if(dmin<=0.) cycle
             !print*,'Extracting from hilbert domain '
             !print*,xmin,xmax
             !print*,ymin,ymax
             !print*,zmin,zmax
             !print*,dom_pmin(1,ix),dom_pmax(1,ix)
             !print*,dom_pmin(2,iy),dom_pmax(2,iy)
             !print*,dom_pmin(3,iz),dom_pmax(3,iz)
             !dmax = max(xmax-xmin,ymax-ymin,zmax-zmin)
             dmax = max(dom_pmax(1,ix)-dom_pmin(1,ix)  &
                  ,dom_pmax(2,iy)-dom_pmin(2,iy)  &
                  ,dom_pmax(3,iz)-dom_pmin(3,iz) )
             
             ! Set up at most 8 sub-domains in periodic box
             do ilevel=1,lmax
                dx=0.5d0**ilevel
                if(dx.lt.dmax)exit
             end do
             lmin=ilevel
             bit_length=lmin-1
             maxdom=2**bit_length
             imin=0; imax=0; jmin=0; jmax=0; kmin=0; kmax=0
             if(bit_length>0)then
                imin=int(dom_pmin(1,ix)*dble(maxdom))
                imax=imin+1
                jmin=int(dom_pmin(2,iy)*dble(maxdom))
                jmax=jmin+1
                kmin=int(dom_pmin(3,iz)*dble(maxdom))
                kmax=kmin+1
             endif

             !dkey=(dble(2**(lmax+1)/dble(maxdom)))**ndim
             dkey=(real(2**(lmax+1),kind=qdp)/real(maxdom,kind=qdp))**ndim
             ndom=1
             if(bit_length>0)ndom=8
             idom(1)=imin; idom(2)=imax
             idom(3)=imin; idom(4)=imax
             idom(5)=imin; idom(6)=imax
             idom(7)=imin; idom(8)=imax
             jdom(1)=jmin; jdom(2)=jmin
             jdom(3)=jmax; jdom(4)=jmax
             jdom(5)=jmin; jdom(6)=jmin
             jdom(7)=jmax; jdom(8)=jmax
             kdom(1)=kmin; kdom(2)=kmin
             kdom(3)=kmin; kdom(4)=kmin
             kdom(5)=kmax; kdom(6)=kmax
             kdom(7)=kmax; kdom(8)=kmax
             
             do i=1,ndom
                if(bit_length>0)then
                   call hilbert3d(idom(i),jdom(i),kdom(i),order_min(i),bit_length,1)
                else
                   order_min(i)=0.0d0
                endif
                bounding_min(i)=(order_min(i))*dkey
                bounding_max(i)=(order_min(i)+1.0D0)*dkey
             end do
             
             cpu_min=0; cpu_max=0
             do impi=1,ncpu
                do i=1,ndom
                   if (   bound_key(impi-1).le.bounding_min(i).and.&
                        & bound_key(impi  ).gt.bounding_min(i))then
                      cpu_min(i)=impi
                   endif
                   if (   bound_key(impi-1).lt.bounding_max(i).and.&
                        & bound_key(impi  ).ge.bounding_max(i))then
                      cpu_max(i)=impi
                   endif
                end do
             end do
             
             do i=1,ndom
                do j=cpu_min(i),cpu_max(i)
                   if(.not. cpu_read(j))then
                      ncpu_read=ncpu_read+1
                      cpu_list(ncpu_read)=j
                      cpu_read(j)=.true.
                   endif
                enddo
             enddo
             
          end do !ix=1,2
       end do !iy=1,2
    end do !iz=1,2

    deallocate(bound_key,cpu_read)

    !print*,'--> nCPU to read = ',ncpu_read, cpu_list(1:ncpu_reaD)
    
    return

  end subroutine get_cpu_list_periodic


  subroutine read_hilbert_keys(repository,snapnum,ncpu,bound_key)
    ! read the hilbert keys in the info file
    
    implicit none

    character(2000),intent(in)                   :: repository
    integer(kind=4),intent(in)                   :: snapnum, ncpu
    real(kind=8),dimension(0:ncpu),intent(inout) :: bound_key

    logical(kind=4)            :: not_ok
    character(512)             :: filename
    character(512)             :: line,name,value,orderingtype
    integer(kind=4)            :: i, impi
    integer(kind=4),parameter  :: param_unit = 13

    not_ok = .true.
    write(filename,'(a,a,i5.5,a,i5.5,a)') trim(repository),'/output_',snapnum,'/info_',snapnum,'.txt'
    open(unit=param_unit,file=filename,status='old',form='formatted')
    do
       read(param_unit,'(a)',end=2) line
       i = scan(line,'=')
       if (i==0 .or. line(1:1)=='#') cycle
       name=trim(adjustl(line(:i-1)))
       value=trim(adjustl(line(i+1:)))
       ! check for a comment at end of line !
       i = scan(value,'!')
       if (i /= 0) value = trim(adjustl(value(:i-1)))
       if (trim(name) .eq. 'ordering type') then
          read(value,*) orderingtype
          if (orderingtype=='hilbert') then
             ! start reading keys
             ! skip one line
             read(param_unit,*)
             do impi=1,ncpu
                read(param_unit,'(I8,1X,E23.15,1X,E23.15)')i,bound_key(impi-1),bound_key(impi)
             end do
          end if
       end if
    end do
2   close (param_unit)
    
    return
    
  end subroutine read_hilbert_keys


  subroutine read_hilbert_keys_raw(repository,snapnum,ncpu,bound_key)
    ! read hilbert keys in one amr file
    
    implicit none
    
    character(2000),intent(in)                :: repository
    integer(kind=4),intent(in)                :: snapnum, ncpu
    real(qdp),dimension(0:ncpu),intent(inout) :: bound_key
    real(kind=8),dimension(0:ncpu)            :: bound_key_dp    
    logical(kind=4)                           :: ok
    character(512)                            :: nomfich
    character(128)                            :: orderingtype
    integer(kind=4)                           :: i, ios
    integer(kind=4),parameter                 :: param_unit = 13
    
    write(nomfich,'(a,a,i5.5,a,i5.5,a)') trim(repository), '/output_', snapnum, '/amr_', snapnum, '.out00001'
    inquire(file=nomfich, exist=ok)
    if(.not. ok)then
       write(*,*)'File '//TRIM(nomfich)//' not found'    
       stop
    end if
    open(unit=param_unit,file=nomfich,form='unformatted',status='old',action='read',iostat=ios)
    do i=1,24 ! Assume that there is no "simple boundary"
       read(param_unit,iostat=ios)
       !print*,'ios =',ios
    end do
    read(param_unit,iostat=ios) orderingtype
    !print*,'ios =',ios

    if (trim(orderingtype) .ne. 'bisection') then
       read(param_unit,iostat=ios) bound_key
       !print*,'ios bk =',ios
       if(ios/=0) then ! read in dp
          print*,'Reading Hilbert keys in quad precision failed, read them in double precision...'
          close(param_unit)
          open(unit=param_unit,file=nomfich,form='unformatted',status='old',action='read',iostat=ios)
          do i=1,24 ! Assume that there is no "simple boundary"
             read(param_unit,iostat=ios)
             !print*,'ios =',ios
          end do
          read(param_unit,iostat=ios) orderingtype
          !print*,'ios =',ios
          read(param_unit,iostat=ios) bound_key_dp
          !print*,'ios bk =',ios
          if(ios/=0) then
             print*,'Reading Hilbert keys in double precision failed, read them in the info file...'
             call read_hilbert_keys(repository,snapnum,ncpu,bound_key_dp)
          end if
          bound_key = real(bound_key_dp, kind=qdp)
       end if
     end if

    close (param_unit)
    
    return
    
  end subroutine read_hilbert_keys_raw



  subroutine hilbert3d(x,y,z,order,bit_length,npoint)
    implicit none

    integer     ,INTENT(IN)                     ::bit_length,npoint
    integer     ,INTENT(IN) ,dimension(1:npoint)::x,y,z
    real(qdp),INTENT(OUT),dimension(1:npoint)::order

    logical,dimension(0:3*bit_length-1)::i_bit_mask
    logical,dimension(0:1*bit_length-1)::x_bit_mask,y_bit_mask,z_bit_mask
    integer,dimension(0:7,0:1,0:11)::state_diagram
    integer::i,ip,cstate,nstate,b0,b1,b2,sdigit,hdigit

    if(bit_length>bit_size(bit_length))then
       write(*,*)'Maximum bit length=',bit_size(bit_length)
       write(*,*)'stop in hilbert3d'
       stop
    endif

    state_diagram = RESHAPE( (/   1, 2, 3, 2, 4, 5, 3, 5,&
         &   0, 1, 3, 2, 7, 6, 4, 5,&
         &   2, 6, 0, 7, 8, 8, 0, 7,&
         &   0, 7, 1, 6, 3, 4, 2, 5,&
         &   0, 9,10, 9, 1, 1,11,11,&
         &   0, 3, 7, 4, 1, 2, 6, 5,&
         &   6, 0, 6,11, 9, 0, 9, 8,&
         &   2, 3, 1, 0, 5, 4, 6, 7,&
         &  11,11, 0, 7, 5, 9, 0, 7,&
         &   4, 3, 5, 2, 7, 0, 6, 1,&
         &   4, 4, 8, 8, 0, 6,10, 6,&
         &   6, 5, 1, 2, 7, 4, 0, 3,&
         &   5, 7, 5, 3, 1, 1,11,11,&
         &   4, 7, 3, 0, 5, 6, 2, 1,&
         &   6, 1, 6,10, 9, 4, 9,10,&
         &   6, 7, 5, 4, 1, 0, 2, 3,&
         &  10, 3, 1, 1,10, 3, 5, 9,&
         &   2, 5, 3, 4, 1, 6, 0, 7,&
         &   4, 4, 8, 8, 2, 7, 2, 3,&
         &   2, 1, 5, 6, 3, 0, 4, 7,&
         &   7, 2,11, 2, 7, 5, 8, 5,&
         &   4, 5, 7, 6, 3, 2, 0, 1,&
         &  10, 3, 2, 6,10, 3, 4, 4,&
         &   6, 1, 7, 0, 5, 2, 4, 3 /), &
         & (/8 ,2, 12 /) )

    do ip=1,npoint

       ! convert to binary
       do i=0,bit_length-1
          x_bit_mask(i)=btest(x(ip),i)
          y_bit_mask(i)=btest(y(ip),i)
          z_bit_mask(i)=btest(z(ip),i)
       enddo
       
       ! interleave bits
       do i=0,bit_length-1
          i_bit_mask(3*i+2)=x_bit_mask(i)
          i_bit_mask(3*i+1)=y_bit_mask(i)
          i_bit_mask(3*i  )=z_bit_mask(i)
       end do
       
       ! build Hilbert ordering using state diagram
       cstate=0
       do i=bit_length-1,0,-1
          b2=0 ; if(i_bit_mask(3*i+2))b2=1
          b1=0 ; if(i_bit_mask(3*i+1))b1=1
          b0=0 ; if(i_bit_mask(3*i  ))b0=1
          sdigit=b2*4+b1*2+b0
          nstate=state_diagram(sdigit,0,cstate)
          hdigit=state_diagram(sdigit,1,cstate)
          i_bit_mask(3*i+2)=btest(hdigit,2)
          i_bit_mask(3*i+1)=btest(hdigit,1)
          i_bit_mask(3*i  )=btest(hdigit,0)
          cstate=nstate
       enddo

       ! save Hilbert key as double precision real
       order(ip)=0.
       do i=0,3*bit_length-1
          b0=0 ; if(i_bit_mask(i))b0=1
          order(ip)=order(ip)+real(b0,kind=qdp)*real(2,kind=qdp)**i
       end do
       
    end do

  end subroutine hilbert3d



  function get_nGridTot_cpus(repository,snapnum,ncpu_read,cpu_list)

    ! get total number of grids in the simulation belonging to given list of cpus

    implicit none 

    integer(kind=4),intent(in)  :: snapnum, ncpu_read
    character(1000),intent(in)  :: repository
    integer(kind=4),dimension(:),allocatable,intent(in) :: cpu_list
    integer(kind=4)             :: get_nGridTot_cpus
    character(1000)             :: filename
    logical                     :: ok
    integer(kind=4)             :: k,icpu,ngrid_current

    get_nGridTot_cpus = 0
    do k = 1,ncpu_read
       icpu=cpu_list(k)
       write(filename,'(a,a,i5.5,a,i5.5,a,i5.5)') trim(repository),'/output_',snapnum,'/amr_',snapnum,'.out',icpu
       inquire(file=filename, exist=ok)
       if(.not. ok)then
          write(*,*)'File '//TRIM(filename)//' not found'    
          stop
       end if
       open(unit=10,file=filename,form='unformatted',status='old',action='read')
       read(10)
       read(10)
       read(10)
       read(10)
       read(10)
       read(10)
       read(10)ngrid_current
       close(10)
       get_nGridTot_cpus = get_nGridTot_cpus + ngrid_current
    end do

    return

  end function get_nGridTot_cpus


  subroutine ramses_get_T_nhi_cgs(repository,snapnum,nleaf,nvar,ramses_var,temp,nhi)

    implicit none 

    character(1000),intent(in)  :: repository
    integer(kind=4),intent(in)  :: snapnum
    integer(kind=4),intent(in)  :: nleaf, nvar
    real(kind=8),intent(in)     :: ramses_var(nvar,nleaf) ! one cell only
    real(kind=8),intent(inout)  :: nhi(nleaf), temp(nleaf)
    real(kind=8),allocatable    :: boost(:)
    real(kind=8)                :: xhi
    integer(kind=4) :: ihx,ihy,i
    real(kind=8)    :: xx,yy,dxx1,dxx2,dyy1,dyy2,f
    integer(kind=4) :: if1,if2,jf1,jf2

    real(kind=8),allocatable,dimension(:)    :: mu

    ! get conversion factors if necessary
    if (.not. conversion_scales_are_known) then 
       call read_conversion_scales(repository,snapnum)
       conversion_scales_are_known = .True.
    end if

    if(ramses_rt)then
       ! ramses RT
       allocate(mu(1:nleaf))
       nhi  = ramses_var(1,:) * dp_scale_nh  * (1.d0 - ramses_var(ihii,:))   ! nb of H atoms per cm^3
       mu   = XH * (1.d0+ramses_var(ihii,:)) + 0.25d0*(1.d0-XH)*(1.d0 + ramses_var(iheii,:) + 2.d0*ramses_var(iheiii,:)) ! assumes no metals
       mu   = 1.0d0 / mu   
       temp = ramses_var(itemp,:) / ramses_var(1,:) * dp_scale_T2            ! T/mu [ K ]
       temp = temp * mu                                                      ! This is now T (in K) with no bloody mu ... 
       deallocate(mu)
    else
       ! ramses standard

       if (.not. cooling_is_read) then
          call read_cooling(repository,snapnum)
          cooling_is_read = .True.
       end if

       nhi  = ramses_var(1,:) * dp_scale_nh  ! nb of H atoms per cm^3
       temp = ramses_var(itemp,:) / ramses_var(1,:) * dp_scale_T2  ! T/mu [ K ]

       allocate(boost(nleaf))
       if (self_shielding) then
          do i=1,nleaf
             boost(i)=MAX(exp(-nhi(i)/0.01),1.0D-20) ! same as hard-coded in RAMSES. 
          end do
       else
          boost = 1.0d0
       end if

    
       ! compute the ionization state and temperature using the 'cooling' tables
       do i = 1, nleaf
          xx  = min(max(log10(nhi(i)/boost(i)),cooling%nh(1)),cooling%nh(cooling%n11))
          ihx = int((xx - cool_int%nh_start)/cool_int%nh_step) + 1
          if (ihx < 1) then 
             ihx = 1 
          else if (ihx > cool_int%n_nh) then
             ihx = cool_int%n_nh
          end if
          yy  = log10(temp(i))
          ihy = int((yy - cool_int%t2_start)/cool_int%t2_step) + 1
          if (ihy < 1) then 
             ihy = 1 
          else if (ihy > cool_int%n_t2) then
             ihy = cool_int%n_t2
          end if
          ! 2D linear interpolation:
          if (ihx < cool_int%n_nh) then 
             dxx1  = max(xx - cooling%nh(ihx),0.0d0) / cool_int%nh_step 
             dxx2  = min(cooling%nh(ihx+1) - xx,cool_int%nh_step) / cool_int%nh_step
             if1  = ihx
             if2  = ihx+1
          else
             dxx1  = 0.0d0
             dxx2  = 1.0d0
             if1  = ihx
             if2  = ihx
          end if
          if (ihy < cool_int%n_t2) then 
             dyy1  = max(yy - cooling%t2(ihy),0.0d0) / cool_int%t2_step
             dyy2  = min(cooling%t2(ihy+1) - yy,cool_int%t2_step) / cool_int%t2_step
             jf1  = ihy
             jf2  = ihy + 1
          else
             dyy1  = 0.0d0
             dyy2  = 1.0d0
             jf1  = ihy
             jf2  = ihy
          end if
          if (abs(dxx1+dxx2-1.0d0) > 1.0d-6 .or. abs(dyy1+dyy2-1.0d0) > 1.0d-6) then 
             write(*,*) 'Fucked up the interpolation ... '
             print*,dxx1+dxx2,dyy1+dyy2
             stop
          end if
          ! neutral H density 
          f = dxx1 * dyy1 * cooling%spec(if2,jf2,2) + dxx2 * dyy1 * cooling%spec(if1,jf2,2) &
               & + dxx1 * dyy2 * cooling%spec(if2,jf1,2) + dxx2 * dyy2 * cooling%spec(if1,jf1,2)
          xhi = 10.0d0**(f-xx)  ! this is xHI = nHI/nH where the n's include the boost. 
          nhi(i) = xhi * nhi(i) ! nHI (cm^-3) (boost-free)
          ! GET MU to convert T/MU into T ... 
          f = dxx1 * dyy1 * cooling%mu(if2,jf2) + dxx2 * dyy1 * cooling%mu(if1,jf2) &
               & + dxx1 * dyy2 * cooling%mu(if2,jf1) + dxx2 * dyy2 * cooling%mu(if1,jf1)
          temp(i) = temp(i) * f   ! This is now T (in K) with no bloody mu ... 
       end do

       deallocate(boost)

    endif
    
    return

  end subroutine ramses_get_T_nhi_cgs

  subroutine ramses_get_velocity_cgs(repository,snapnum,nleaf,nvar,ramses_var,velocity_cgs)

    implicit none

    character(1000),intent(in)  :: repository
    integer(kind=4),intent(in)  :: snapnum
    integer(kind=4),intent(in)  :: nleaf, nvar
    real(kind=8),intent(in)     :: ramses_var(nvar,nleaf) ! one cell only
    real(kind=8),intent(inout)  :: velocity_cgs(3,nleaf)

    ! get conversion factors if necessary
    if (.not. conversion_scales_are_known) then 
       call read_conversion_scales(repository,snapnum)
       conversion_scales_are_known = .True.
    end if

    velocity_cgs = ramses_var(2:4,:) * dp_scale_v ! [ cm / s ]

    return

  end subroutine ramses_get_velocity_cgs

  subroutine ramses_get_nh_cgs(repository,snapnum,nleaf,nvar,ramses_var,nh)

    implicit none

    character(1000),intent(in)  :: repository
    integer(kind=4),intent(in)  :: snapnum
    integer(kind=4),intent(in)  :: nleaf, nvar
    real(kind=8),intent(in)     :: ramses_var(nvar,nleaf) ! one cell only
    real(kind=8),intent(inout)  :: nh(nleaf)

    ! get conversion factors if necessary
    if (.not. conversion_scales_are_known) then 
       call read_conversion_scales(repository,snapnum)
       conversion_scales_are_known = .True.
    end if

    nh = ramses_var(1,:) * dp_scale_nh ! [ H / cm^3 ]

    return

  end subroutine ramses_get_nh_cgs


  subroutine ramses_get_metallicity(nleaf,nvar,ramses_var,metallicity)

    implicit none

    integer(kind=4),intent(in)  :: nleaf, nvar
    real(kind=8),intent(in)     :: ramses_var(nvar,nleaf) ! one cell only
    real(kind=8),intent(inout)  :: metallicity(nleaf)

    if (nvar < 6) then
       print*,'No metals !!! '
       stop
    end if
    metallicity = ramses_var(imetal,:) 

    return

  end subroutine ramses_get_metallicity


  
  subroutine ramses_get_cooling_time(repository,snapnum,nleaf,nvar,ramses_var,coolingTime,sample)

    implicit none 
    character(1000),intent(in)  :: repository
    integer(kind=4),intent(in)  :: snapnum
    integer(kind=4),intent(in) :: nleaf,nvar
    real(kind=8),intent(in)    :: ramses_var(nvar,nleaf)
    real(kind=8),intent(inout) :: coolingTime(:)
    integer(kind=4),intent(in),optional :: sample(:)
    
    real(kind=8)               :: aexp,xhii,xheii,xheiii,nh,nhi,nhii,nhe,nhei,nheii,nheiii,ne,mu,T,crate,enerth,dcooldT
    integer(kind=4)            :: n,j,i
    logical                    :: subsample
    ! Heating terms (RT)
    real(kind=8),parameter::eV_to_erg=1.6022d-12  ! eV to erg conv. constant
    !!JB- should read these from files ! 
    !real(kind=8),parameter,dimension(3)::ion_egy = (/13.60d0, 24.59d0, 54.42d0/)*eV_to_erg
    real(kind=8),parameter,dimension(4)::ion_egy = (/13.60d0, 15.20d0, 24.59d0, 54.42d0/)*eV_to_erg
    !! -JB
    character(1000)            :: filename
    integer(kind=4)            :: ilun=33,nRTvar,nIons,nGroups,igroup,indexgroup,nvarH
    real(kind=8),allocatable   :: group_egy(:),group_csn(:,:),group_cse(:,:)
    real(kind=8)               :: hrate,unit_fp
    

    if (.not. ramses_rt) then
       print*,'Cooling time not implemented without RT... '
       stop
    end if
    
    ! get conversion factors if necessary
    if (.not. conversion_scales_are_known) then 
       call read_conversion_scales(repository,snapnum)
       conversion_scales_are_known = .True.
    end if

    ! init cooling tables 
    aexp    = get_param_real(repository,snapnum,'aexp')
    call init_coolrates_tables(aexp)

    ! check Photo-heating
    if (read_rt_variables) then
       write(filename,'(a,a,i5.5,a,i5.5,a)') TRIM(repository),'output_',snapnum,'/info_rt_',snapnum,'.txt'
       open(unit=ilun,file=filename,status='old',form='formatted')
       call read_int( ilun, 'nRTvar', nRTvar)
       nvarH = nvar - nRTvar
       call read_int( ilun, 'nIons', nIons)
! JB: make sure this is OK for non-Harley cases ... 
!!$       if (nIons .ne. 3) then
!!$          print*,'nIons has to be 3 with current implementation ... '
!!$          stop
!!$       end if
       call read_int( ilun, 'nGroups', nGroups)
       allocate(group_egy(nGroups),group_csn(ngroups,nions),group_cse(ngroups,nions))
       call read_real(ilun, 'unit_pf', unit_fp)
       call read_groups(ilun)
       close(ilun)
    end if

    
    ! compute cooling rate for all sample cells 
    if (present(sample)) then
       n = size(sample)
       subsample = .true.
    else
       n = nleaf
       subsample = .false. 
    end if
    do j=1,n
       if (subsample) then
          i = sample(j)
       else
          i = j
       end if
       xhii   = ramses_var(ihii,i)
       xheii  = ramses_var(iheii,i)
       xheiii = ramses_var(iheiii,i)
       nh     = ramses_var(1,i) * dp_scale_nh
       nhi    = nh * (1.0d0 - xhii)
       nhii   = nh * xhii
       nhe    = 0.25*nh*(1.0d0-XH)/XH
       nhei   = nhe * (1.0d0 - xheii - xheiii)
       nheii  = nhe * xheii
       nheiii = nhe * xheiii
       ne     = nHII + nHe * (xHeII + 2.0d0*xHeIII)
       mu     = 1./( XH*(1.+xHII) + 0.25d0*(1.0d0-XH)*(1.+xHeII+2.*xHeIII) )
       T      = ramses_var(itemp,i)/ramses_var(1,i)*mu*dp_scale_T2
       crate  = compCoolrate(T, ne, nHI, nHII, nHeI, nHeII, nHeIII, aexp, dcooldT)  ! [erg s-1 cm-3]  
       if (read_rt_variables) then
          hrate = 0.0d0
          ! JB- There has to be a way to use non-hard-coded indexes ... 
          if (nions.eq.3) then
             do igroup=1,ngroups
                indexgroup = nvarH+1+(igroup-1)*(1+ndim)
                hrate = hrate + nhi * ramses_var(indexgroup,i) * unit_fp * (group_cse(igroup,1)*group_egy(igroup) -group_csn(iGroup,1)*ion_egy(1)) &
                     & + nhei * ramses_var(indexgroup,i) * unit_fp * (group_cse(igroup,2)*group_egy(igroup) -group_csn(iGroup,2)*ion_egy(2)) &
                     & + nheii * ramses_var(indexgroup,i) * unit_fp * (group_cse(igroup,3)*group_egy(igroup) -group_csn(iGroup,3)*ion_egy(3)) 
             end do
          else !!HK addition for nIons=4 (molecular hydrogen case)
             do igroup=1,ngroups
                indexgroup = nvarH+1+(igroup-1)*(1+ndim)
                hrate = hrate + nhi * ramses_var(indexgroup,i) * unit_fp * (group_cse(igroup,2)*group_egy(igroup) -group_csn(iGroup,2)*ion_egy(2)) &
                     & + nhei * ramses_var(indexgroup,i) * unit_fp * (group_cse(igroup,3)*group_egy(igroup) -group_csn(iGroup,3)*ion_egy(3)) &
                     & + nheii * ramses_var(indexgroup,i) * unit_fp * (group_cse(igroup,4)*group_egy(igroup) -group_csn(iGroup,4)*ion_egy(4))
             end do
          endif
          ! -JB 
          crate = max(1.0d-40,crate-hrate)  ! we're only interested in relatively fast cooling rates. 
       end if
       enerth = (1.5d0 * kb / mp ) * T / mu * ramses_var(1,i) * dp_scale_d  ! [erg cm-3]
       coolingTime(j) = enerth / crate ! [s]
    end do


    if (read_rt_variables) then 
       deallocate(group_egy,group_csn,group_cse)
    end if
       
    return


  contains
    
    !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    SUBROUTINE read_int(lun, param_name, value)
      ! Try to read a parameter from lun
      !-------------------------------------------------------------------------
      integer::lun
      character(*)::param_name
      character(128)::line,tmp
      integer::value
      !-------------------------------------------------------------------------
      rewind(unit=lun)
      do
         read(lun, '(A128)', end=223) line
         if(index(line,trim(param_name)) .eq. 1) then
            read(line,'(A13,I30)') tmp, value
            return
         endif
      end do
223   return                        ! eof reached, didn't find the parameter

    END SUBROUTINE read_int
    !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    SUBROUTINE read_real_array(lun, param_name, value)

      ! Try to read a parameter array from lun
      !-------------------------------------------------------------------------
      integer::lun
      character(*)::param_name
      character(1000)::line,tmp
      real(kind=8),dimension(:)::value
      !-------------------------------------------------------------------------
      rewind(unit=lun)
      do
         read(lun, '(A1000)', end=222) line
         if(index(line,trim(param_name)) .eq. 1) then
            read(line,'(A13,100(E23.15))') tmp, value
            return
         endif
      end do
222   return                        ! eof reached, didn't find the parameter

    END SUBROUTINE read_real_array
    !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    SUBROUTINE read_real(lun, param_name, value)

      ! Try to read a parameter from lun
      !-------------------------------------------------------------------------
      integer::lun
      character(*)::param_name
      character(1000)::line,tmp
      real(kind=8)::value
      !-------------------------------------------------------------------------
      rewind(unit=lun)
      do
         read(lun, '(A128)', end=222) line
         if(index(line,trim(param_name)) .eq. 1) then
            read(line,'(A13,E23.15)') tmp, value
            return
         endif
      end do
222   return                        ! eof reached, didn't find the parameter

    END SUBROUTINE read_real
    !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    SUBROUTINE read_groups(lun)

      ! Try to read photon group properties from lun
      !-------------------------------------------------------------------------
      integer::lun,i
      character(128)::line,tmp
      real ::group_tmp(3)
      !-------------------------------------------------------------------------
      rewind(unit=lun)
      i=0 ! Read group_egy
      do
         read(lun, '(A128)', end=220) line
         if(index(line,trim('egy')) .eq. 3) then
            i=i+1
            read(line,'(A19,100F12.3)') tmp, group_egy(i)
         endif
      end do
220   continue
      group_egy=group_egy*eV_to_erg
      i=0 ! Read group_csn
      rewind(unit=lun)
      do
         read(lun, '(A128)', end=221) line
         if(index(line,trim('csn')) .eq. 3) then
            i=i+1
            read(line,'(A19,3(1pe12.3))') tmp, group_tmp
            group_csn(i,1:nIons)=group_tmp
         endif
      end do
221   continue
      i=0 ! Read group_cse
      rewind(unit=lun)
      do
         read(lun, '(A128)', end=222) line
         if(index(line,trim('cse')) .eq. 3) then
            i=i+1
            read(line,'(A19,100(1pe12.3))') tmp, group_tmp
            group_cse(i,1:nIons) = group_tmp
         endif
      end do
222   continue
      return                        ! eof reached, didn't find the parameter

    END SUBROUTINE read_groups
    !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

  end subroutine ramses_get_cooling_time


  subroutine ramses_get_LyaEmiss_HIDopwidth(repository,snapnum,nleaf,nvar,var,recomb_em,coll_em,HIDopwidth,sample)

    implicit none

    character(1000),intent(in)  :: repository
    integer(kind=4),intent(in)  :: snapnum
    integer(kind=4),intent(in) :: nleaf,nvar
    real(kind=8),intent(in)    :: var(nvar,nleaf)
    real(kind=8),intent(inout) :: recomb_em(:),coll_em(:),HIDopwidth(:)
    integer(kind=4),intent(in),optional :: sample(:)

    integer(kind=4)            :: n,j,i
    real(kind=8),parameter     :: e_lya = planck * clight / (1215.67d0/cmtoA) ! [erg] energy of a Lya photon (consistent with HI_model)
    real(kind=8)               :: xhii,xheii,xheiii,nh,nhi,nhii,n_e,mu,TK,Ta,prob_case_B,alpha_B,collExrate_HI,lambda,nhe
    logical                    :: subsample
    
    ! get conversion factors if necessary
    if (.not. conversion_scales_are_known) then 
       call read_conversion_scales(repository,snapnum)
       conversion_scales_are_known = .True.
    end if

    if (present(sample)) then
       n = size(sample)
       subsample = .true.
    else
       n = nleaf
       subsample = .false. 
    end if

    if(ramses_rt)then
       do j=1,n

          if (subsample) then
             i = sample(j)
          else
             i = j
          end if

          xhii   = var(ihii,i)
          xheii  = var(iheii,i)
          xheiii = var(iheiii,i)
          nh     = var(1,i) * dp_scale_nh
          nhe    = 0.25*nh*(1.0d0-XH)/XH
          nhii   = nh * xhii
          nhi    = nh * (1.0d0 - xhii)
          n_e    = nHII + nHe * (xHeII + 2.0d0*xHeIII)
          mu     = 1./( XH*(1.+xHII) + 0.25d0*(1.0d0-XH)*(1.+xHeII+2.*xHeIII) )
          TK     = var(itemp,i)/var(1,i)*mu*dp_scale_T2
          HIDopwidth(j) = sqrt((2.0d0*kb/mp)*TK)
          ! Cantalupo+(08)
          Ta = max(TK,100.0) ! no extrapolation..
          prob_case_B = 0.686 - 0.106*log10(Ta/1e4) - 0.009*(Ta/1e4)**(-0.44)
          ! Hui & Gnedin (1997)
          lambda = 315614.d0/TK
          alpha_B = 2.753d-14*(lambda**(1.5))/(1+(lambda/2.74)**0.407)**(2.242) ![cm3/s]
          recomb_em(j) = prob_case_B * alpha_B * n_e * nhii * e_lya ! [erg/cm3/s]
          ! collisional emission
          collExrate_HI = 2.41d-6/sqrt(TK) * (TK/1.d4)**0.22 * exp(-1.63d-11/(kB*TK))
          coll_em(j) = nHI * n_e * collExrate_HI * e_lya
       end do
    else
       print*,'Not implemented ... '
       stop
    end if
    
  end subroutine ramses_get_LyaEmiss_HIDopwidth

  
  subroutine ramses_get_T_nSiII_cgs(repository,snapnum,nleaf,nvar,ramses_var,temp,nSiII)

    implicit none 

    character(1000),intent(in)            :: repository
    integer(kind=4),intent(in)            :: snapnum
    integer(kind=4),intent(in)            :: nleaf, nvar
    real(kind=8),intent(in)               :: ramses_var(nvar,nleaf) ! one cell only
    real(kind=8),intent(inout)            :: nSiII(nleaf), temp(nleaf)
    real(kind=8),allocatable              :: boost(:)
    integer(kind=4)                       :: ihx,ihy,i
    real(kind=8)                          :: xx,yy,dxx1,dxx2,dyy1,dyy2,f
    integer(kind=4)                       :: if1,if2,jf1,jf2
    real(kind=8),allocatable,dimension(:) :: mu, nh

    ! get conversion factors if necessary
    if (.not. conversion_scales_are_known) then 
       call read_conversion_scales(repository,snapnum)
       conversion_scales_are_known = .True.
    end if


    if(ramses_rt)then
       ! ramses RT
       allocate(mu(1:nleaf))
       mu   = XH * (1.d0+ramses_var(ihii,:)) + 0.25d0*(1.d0-XH)*(1.d0 + ramses_var(iheii,:) + 2.d0*ramses_var(iheiii,:)) ! assumes no metals
       mu   = 1.0d0 / mu
       temp = ramses_var(itemp,:) / ramses_var(1,:) * dp_scale_T2            ! T/mu [ K ]
       temp = temp * mu                                                      ! This is now T (in K) with no bloody mu ... 
       deallocate(mu)
    else
       ! ramses standard
       if (.not. cooling_is_read) then
          call read_cooling(repository,snapnum)
          cooling_is_read = .True.
       end if
       allocate(nh(nleaf))
       nh   = ramses_var(1,:) * dp_scale_nh  ! nb of H atoms per cm^3
       temp = ramses_var(itemp,:) / ramses_var(1,:) * dp_scale_T2  ! T/mu [ K ]
       allocate(boost(nleaf))
       if (self_shielding) then
          do i=1,nleaf
             boost(i)=MAX(exp(-nh(i)/0.01),1.0D-20) ! same as hard-coded in RAMSES. 
          end do
       else
          boost = 1.0d0
       end if
       ! compute the ionization state and temperature using the 'cooling' tables
       do i = 1, nleaf
          xx  = min(max(log10(nh(i)/boost(i)),cooling%nh(1)),cooling%nh(cooling%n11))
          ihx = int((xx - cool_int%nh_start)/cool_int%nh_step) + 1
          if (ihx < 1) then 
             ihx = 1 
          else if (ihx > cool_int%n_nh) then
             ihx = cool_int%n_nh
          end if
          yy  = log10(temp(i))
          ihy = int((yy - cool_int%t2_start)/cool_int%t2_step) + 1
          if (ihy < 1) then 
             ihy = 1 
          else if (ihy > cool_int%n_t2) then
             ihy = cool_int%n_t2
          end if
          ! 2D linear interpolation:
          if (ihx < cool_int%n_nh) then 
             dxx1  = max(xx - cooling%nh(ihx),0.0d0) / cool_int%nh_step 
             dxx2  = min(cooling%nh(ihx+1) - xx,cool_int%nh_step) / cool_int%nh_step
             if1  = ihx
             if2  = ihx+1
          else
             dxx1  = 0.0d0
             dxx2  = 1.0d0
             if1  = ihx
             if2  = ihx
          end if
          if (ihy < cool_int%n_t2) then 
             dyy1  = max(yy - cooling%t2(ihy),0.0d0) / cool_int%t2_step
             dyy2  = min(cooling%t2(ihy+1) - yy,cool_int%t2_step) / cool_int%t2_step
             jf1  = ihy
             jf2  = ihy + 1
          else
             dyy1  = 0.0d0
             dyy2  = 1.0d0
             jf1  = ihy
             jf2  = ihy
          end if
          if (abs(dxx1+dxx2-1.0d0) > 1.0d-6 .or. abs(dyy1+dyy2-1.0d0) > 1.0d-6) then 
             write(*,*) 'Fucked up the interpolation ... '
             print*,dxx1+dxx2,dyy1+dyy2
             stop
          end if
          ! GET MU to convert T/MU into T ... 
          f = dxx1 * dyy1 * cooling%mu(if2,jf2) + dxx2 * dyy1 * cooling%mu(if1,jf2) &
               & + dxx1 * dyy2 * cooling%mu(if2,jf1) + dxx2 * dyy2 * cooling%mu(if1,jf1)
          temp(i) = temp(i) * f   ! This is now T (in K) with no bloody mu ... 
       end do
       deallocate(boost,nh)
    endif       
    
    ! from T, we can compute the SiII fraction as 100% between 6.306571e+04 K and 1.264600e+05 K, and 0% elsewhere.
    ! (These limits correspond to the following ionisation energies:
    ! - Si-Si+   : 8.15169 ev
    ! - Si+-Si++ : 16.34585 eV
    do i = 1,nleaf
       if (temp(i) >= 6.306571d4 .and. temp(i) <= 1.2646d5) then
          nSiII(i) = ramses_var(1,i) * dp_scale_d * ramses_var(imetal,i) * dp_scale_zsun * abundance_Si_number   ! [#/cm3]
       else
          nSiII(i) = 0.0d0
       end if
    end do
    
    return

  end subroutine ramses_get_T_nSiII_cgs



  subroutine ramses_get_T_nMgII_cgs(repository,snapnum,nleaf,nvar,ramses_var,temp,nMgII)

    implicit none 

    character(1000),intent(in)            :: repository
    integer(kind=4),intent(in)            :: snapnum
    integer(kind=4),intent(in)            :: nleaf, nvar
    real(kind=8),intent(in)               :: ramses_var(nvar,nleaf) ! one cell only
    real(kind=8),intent(inout)            :: nMgII(nleaf), temp(nleaf)
    real(kind=8),allocatable              :: boost(:)
    integer(kind=4)                       :: ihx,ihy,i
    real(kind=8)                          :: xx,yy,dxx1,dxx2,dyy1,dyy2,f
    integer(kind=4)                       :: if1,if2,jf1,jf2
    real(kind=8),allocatable,dimension(:) :: mu, nh

    ! get conversion factors if necessary
    if (.not. conversion_scales_are_known) then 
       call read_conversion_scales(repository,snapnum)
       conversion_scales_are_known = .True.
    end if

    if(ramses_rt)then
       ! ramses RT
       allocate(mu(1:nleaf))
       mu   = XH * (1.d0+ramses_var(ihii,:)) + 0.25d0*(1.d0-XH)*(1.d0 + ramses_var(iheii,:) + 2.d0*ramses_var(iheiii,:)) ! assumes no metals
       mu   = 1.0d0 / mu
       temp = ramses_var(itemp,:) / ramses_var(1,:) * dp_scale_T2                ! T/mu [ K ]
       temp = temp * mu                                                      ! This is now T (in K) with no bloody mu ... 
       deallocate(mu)
    else
       ! ramses standard
       if (.not. cooling_is_read) then
          call read_cooling(repository,snapnum)
          cooling_is_read = .True.
       end if
       allocate(nh(nleaf))
       nh   = ramses_var(1,:) * dp_scale_nh  ! nb of H atoms per cm^3
       temp = ramses_var(itemp,:) / ramses_var(1,:) * dp_scale_T2  ! T/mu [ K ]
       allocate(boost(nleaf))
       if (self_shielding) then
          do i=1,nleaf
             boost(i)=MAX(exp(-nh(i)/0.01),1.0D-20) ! same as hard-coded in RAMSES. 
          end do
       else
          boost = 1.0d0
       end if
       ! compute the ionization state and temperature using the 'cooling' tables
       do i = 1, nleaf
          xx  = min(max(log10(nh(i)/boost(i)),cooling%nh(1)),cooling%nh(cooling%n11))
          ihx = int((xx - cool_int%nh_start)/cool_int%nh_step) + 1
          if (ihx < 1) then 
             ihx = 1 
          else if (ihx > cool_int%n_nh) then
             ihx = cool_int%n_nh
          end if
          yy  = log10(temp(i))
          ihy = int((yy - cool_int%t2_start)/cool_int%t2_step) + 1
          if (ihy < 1) then 
             ihy = 1 
          else if (ihy > cool_int%n_t2) then
             ihy = cool_int%n_t2
          end if
          ! 2D linear interpolation:
          if (ihx < cool_int%n_nh) then 
             dxx1  = max(xx - cooling%nh(ihx),0.0d0) / cool_int%nh_step 
             dxx2  = min(cooling%nh(ihx+1) - xx,cool_int%nh_step) / cool_int%nh_step
             if1  = ihx
             if2  = ihx+1
          else
             dxx1  = 0.0d0
             dxx2  = 1.0d0
             if1  = ihx
             if2  = ihx
          end if
          if (ihy < cool_int%n_t2) then 
             dyy1  = max(yy - cooling%t2(ihy),0.0d0) / cool_int%t2_step
             dyy2  = min(cooling%t2(ihy+1) - yy,cool_int%t2_step) / cool_int%t2_step
             jf1  = ihy
             jf2  = ihy + 1
          else
             dyy1  = 0.0d0
             dyy2  = 1.0d0
             jf1  = ihy
             jf2  = ihy
          end if
          if (abs(dxx1+dxx2-1.0d0) > 1.0d-6 .or. abs(dyy1+dyy2-1.0d0) > 1.0d-6) then 
             write(*,*) 'Fucked up the interpolation ... '
             print*,dxx1+dxx2,dyy1+dyy2
             stop
          end if
          ! GET MU to convert T/MU into T ... 
          f = dxx1 * dyy1 * cooling%mu(if2,jf2) + dxx2 * dyy1 * cooling%mu(if1,jf2) &
               & + dxx1 * dyy2 * cooling%mu(if2,jf1) + dxx2 * dyy2 * cooling%mu(if1,jf1)
          temp(i) = temp(i) * f   ! This is now T (in K) with no bloody mu ... 
       end do
       deallocate(boost,nh)
    endif       
    
    ! from T, we can compute the MgII fraction as 100% between 5.915393e+04 K and 1.163181e+05 K, and 0% elsewhere.
    ! (These limits correspond to the following ionisation energies:
    ! - Mg  - Mg+  :  7.646235 ev
    ! - Mg+ - Mg++ : 15.03527 eV
    do i = 1,nleaf
       if (temp(i) >= 5.915393d4 .and. temp(i) <= 1.163181d5) then
          nMgII(i) = ramses_var(1,i) * dp_scale_d * ramses_var(imetal,i) * dp_scale_zsun * abundance_Mg_number   ! [#/cm3]
       else
          nMgII(i) = 0.0d0
       end if
    end do
    
    return

  end subroutine ramses_get_T_nMgII_cgs

  ! Return, nHI, nHeI, nHeII in cells ------------------------------------
  subroutine ramses_get_nhi_nhei_nehii_cgs(repository,snapnum,nleaf,nvar  &
                                               ,ramses_var,nhi,nhei,nheii)

    implicit none 

    character(1000),intent(in)  :: repository
    integer(kind=4),intent(in)  :: snapnum
    integer(kind=4),intent(in)  :: nleaf, nvar
    real(kind=8),intent(in)     :: ramses_var(nvar,nleaf) ! one cell only
    real(kind=8),intent(inout)  :: nhi(nleaf), nhei(nleaf), nheii(nleaf)
    real(kind=8),allocatable    :: nHe(:)

    ! get conversion factors if necessary
    if (.not. conversion_scales_are_known) then 
       call read_conversion_scales(repository,snapnum)
       conversion_scales_are_known = .True.
    end if

    if(ramses_rt)then
       ! ramses RT
       allocate(nHe(nleaf))
       nHe  = ramses_var(1,:) * dp_scale_nhe
       nhi  = ramses_var(1,:) * dp_scale_nh  &      ! nb of HI atoms per cm^3
                              * max(0.d0,(1.d0 - ramses_var(ihii,:)))
       nhei = nHe * max(0.0d0,(1.d0 - ramses_var(iheii,:)  &
                          - ramses_var(iheiii,:)))   ! nb of HeI atoms per cm^3
       nheii  = nHe * (ramses_var(iheii,:))         ! nb of HeII atoms per cm^3
       deallocate(nHe)
    else
       ! ramses standard...nothing to do for now
    endif
    
    return

  end subroutine ramses_get_nhi_nhei_nehii_cgs

  ! Return, nHI, nHeI, nHeII in cells ------------------------------------
  subroutine ramses_get_nh_nhi_nhei_nehii_cgs(repository,snapnum,nleaf,nvar  &
                                               ,ramses_var,nh,nhi,nhei,nheii)

    implicit none 

    character(1000),intent(in)  :: repository
    integer(kind=4),intent(in)  :: snapnum
    integer(kind=4),intent(in)  :: nleaf, nvar
    real(kind=8),intent(in)     :: ramses_var(nvar,nleaf) ! one cell only
    real(kind=8),intent(inout)  :: nh(nleaf), nhi(nleaf), nhei(nleaf), nheii(nleaf)
    real(kind=8),allocatable    :: nHe(:)

    ! get conversion factors if necessary
    if (.not. conversion_scales_are_known) then 
       call read_conversion_scales(repository,snapnum)
       conversion_scales_are_known = .True.
    end if

    if(ramses_rt)then
       ! ramses RT
       allocate(nHe(nleaf))
       nHe  = ramses_var(1,:) * dp_scale_nhe
       nh   = ramses_var(1,:) * dp_scale_nh         ! nb of H atoms per cm^3
       nhi  = ramses_var(1,:) * dp_scale_nh  &      ! nb of HI atoms per cm^3
                              * max(0.d0,(1.d0 - ramses_var(ihii,:)))
       nhei = nHe * max(0.0d0,(1.d0 - ramses_var(iheii,:)  &
                          - ramses_var(iheiii,:)))   ! nb of HeI atoms per cm^3
       nheii  = nHe * (ramses_var(iheii,:))         ! nb of HeII atoms per cm^3
       deallocate(nHe)
    else
       ! ramses standard...nothing to do for now
    endif
    
    return

  end subroutine ramses_get_nh_nhi_nhei_nehii_cgs


  subroutine ramses_get_T_nFeII_cgs(repository,snapnum,nleaf,nvar,ramses_var,temp,nFeII)

    implicit none 

    character(1000),intent(in)            :: repository
    integer(kind=4),intent(in)            :: snapnum
    integer(kind=4),intent(in)            :: nleaf, nvar
    real(kind=8),intent(in)               :: ramses_var(nvar,nleaf) ! one cell only
    real(kind=8),intent(inout)            :: nFeII(nleaf), temp(nleaf)
    real(kind=8),allocatable              :: boost(:)
    integer(kind=4)                       :: ihx,ihy,i
    real(kind=8)                          :: xx,yy,dxx1,dxx2,dyy1,dyy2,f
    integer(kind=4)                       :: if1,if2,jf1,jf2
    real(kind=8),allocatable,dimension(:) :: mu, nh

    ! get conversion factors if necessary
    if (.not. conversion_scales_are_known) then 
       call read_conversion_scales(repository,snapnum)
       conversion_scales_are_known = .True.
    end if

    if(ramses_rt)then
       ! ramses RT
       allocate(mu(1:nleaf))
       mu   = XH * (1.d0+ramses_var(ihii,:)) + 0.25d0*(1.d0-XH)*(1.d0 + ramses_var(iheii,:) + 2.d0*ramses_var(iheiii,:)) ! assumes no metals
       mu   = 1.0d0 / mu
       temp = ramses_var(itemp,:) / ramses_var(1,:) * dp_scale_T2                ! T/mu [ K ]
       temp = temp * mu                                                      ! This is now T (in K) with no bloody mu ... 
       deallocate(mu)
    else
       ! ramses standard
       if (.not. cooling_is_read) then
          call read_cooling(repository,snapnum)
          cooling_is_read = .True.
       end if
       allocate(nh(nleaf))
       nh   = ramses_var(1,:) * dp_scale_nh  ! nb of H atoms per cm^3
       temp = ramses_var(itemp,:) / ramses_var(1,:) * dp_scale_T2  ! T/mu [ K ]
       allocate(boost(nleaf))
       if (self_shielding) then
          do i=1,nleaf
             boost(i)=MAX(exp(-nh(i)/0.01),1.0D-20) ! same as hard-coded in RAMSES. 
          end do
       else
          boost = 1.0d0
       end if
       ! compute the ionization state and temperature using the 'cooling' tables
       do i = 1, nleaf
          xx  = min(max(log10(nh(i)/boost(i)),cooling%nh(1)),cooling%nh(cooling%n11))
          ihx = int((xx - cool_int%nh_start)/cool_int%nh_step) + 1
          if (ihx < 1) then 
             ihx = 1 
          else if (ihx > cool_int%n_nh) then
             ihx = cool_int%n_nh
          end if
          yy  = log10(temp(i))
          ihy = int((yy - cool_int%t2_start)/cool_int%t2_step) + 1
          if (ihy < 1) then 
             ihy = 1 
          else if (ihy > cool_int%n_t2) then
             ihy = cool_int%n_t2
          end if
          ! 2D linear interpolation:
          if (ihx < cool_int%n_nh) then 
             dxx1  = max(xx - cooling%nh(ihx),0.0d0) / cool_int%nh_step 
             dxx2  = min(cooling%nh(ihx+1) - xx,cool_int%nh_step) / cool_int%nh_step
             if1  = ihx
             if2  = ihx+1
          else
             dxx1  = 0.0d0
             dxx2  = 1.0d0
             if1  = ihx
             if2  = ihx
          end if
          if (ihy < cool_int%n_t2) then 
             dyy1  = max(yy - cooling%t2(ihy),0.0d0) / cool_int%t2_step
             dyy2  = min(cooling%t2(ihy+1) - yy,cool_int%t2_step) / cool_int%t2_step
             jf1  = ihy
             jf2  = ihy + 1
          else
             dyy1  = 0.0d0
             dyy2  = 1.0d0
             jf1  = ihy
             jf2  = ihy
          end if
          if (abs(dxx1+dxx2-1.0d0) > 1.0d-6 .or. abs(dyy1+dyy2-1.0d0) > 1.0d-6) then 
             write(*,*) 'Fucked up the interpolation ... '
             print*,dxx1+dxx2,dyy1+dyy2
             stop
          end if
          ! GET MU to convert T/MU into T ... 
          f = dxx1 * dyy1 * cooling%mu(if2,jf2) + dxx2 * dyy1 * cooling%mu(if1,jf2) &
               & + dxx1 * dyy2 * cooling%mu(if2,jf1) + dxx2 * dyy2 * cooling%mu(if1,jf1)
          temp(i) = temp(i) * f   ! This is now T (in K) with no bloody mu ... 
       end do
       deallocate(boost,nh)
    endif       
    
    ! from T, we can compute the FeII fraction as 100% between 9.1704362 d+4 K and 1.87983 d+5 K, and 0% elsewhere.
    ! (These limits correspond to the following ionisation energies:
    ! - Fe - Fe+   : 7.9024678 ev = 1,1215236d-18 J  (NIST data)
    ! - Fe+ - Fe++ : 16.19920 eV = 2,56348d-18 J (NIST data)
    do i = 1,nleaf
       if (temp(i) >= 9.1704362d4 .and. temp(i) <= 1.87983d5) then
          nFeII(i) = ramses_var(1,i) * dp_scale_d * ramses_var(imetal,i) * dp_scale_zsun * abundance_Fe_number   ! [#/cm3]
       else
          nFeII(i) = 0.0d0
       end if
    end do
    
    return

  end subroutine ramses_get_T_nFeII_cgs




  function ramses_get_box_size_cm(repository,snapnum)

    implicit none

    character(1000),intent(in)  :: repository
    integer(kind=4),intent(in)  :: snapnum
    real(kind=8)                :: ramses_get_box_size_cm 

    ! get conversion factors if necessary
    if (.not. conversion_scales_are_known) then 
       call read_conversion_scales(repository,snapnum)
       conversion_scales_are_known = .True.
    end if

    ramses_get_box_size_cm = get_param_real(repository,snapnum,'boxlen') * dp_scale_l  ! [ cm ] 

    return

  end function ramses_get_box_size_cm


  subroutine dump_ramses_info(repository,snapnum,ilun)

    integer(kind=4),intent(in) :: ilun
    integer(kind=4),intent(in) :: snapnum
    character(1000),intent(in) :: repository
    real(kind=8) :: unit_t, unit_l, unit_d, boxlen, time, aexp, H0, omega_m, omega_l, omega_k, omega_b
    
    ! get data
    unit_l  = get_param_real(repository,snapnum,'unit_l')
    unit_d  = get_param_real(repository,snapnum,'unit_d')
    unit_t  = get_param_real(repository,snapnum,'unit_t')
    boxlen  = get_param_real(repository,snapnum,'boxlen')
    time    = get_param_real(repository,snapnum,'time')
    aexp    = get_param_real(repository,snapnum,'aexp')
    H0      = get_param_real(repository,snapnum,'H0')
    omega_m = get_param_real(repository,snapnum,'omega_m')
    omega_l = get_param_real(repository,snapnum,'omega_l')
    omega_k = get_param_real(repository,snapnum,'omega_k')
    omega_b = get_param_real(repository,snapnum,'omega_b')

    ! dump data in open file using unit
    ! Write physical parameters as in RAMSES
    write(ilun,'(" boxlen      =",E23.15)')boxlen
    write(ilun,'(" time        =",E23.15)')time
    write(ilun,'(" aexp        =",E23.15)')aexp
    write(ilun,'(" H0          =",E23.15)')H0
    write(ilun,'(" omega_m     =",E23.15)')omega_m
    write(ilun,'(" omega_l     =",E23.15)')omega_l
    write(ilun,'(" omega_k     =",E23.15)')omega_k
    write(ilun,'(" omega_b     =",E23.15)')omega_b
    write(ilun,'(" unit_l      =",E23.15)')unit_l
    write(ilun,'(" unit_d      =",E23.15)')unit_d
    write(ilun,'(" unit_t      =",E23.15)')unit_t
    write(ilun,*)
  
  end subroutine dump_ramses_info


  !==================================================================================
  ! ----------------
  ! private functions 
  ! ----------------



  subroutine get_leaf_cells_per_cpu(repository,snapnum,icpu,&
       & nleaf,ncell,cell_pos,cell_var,cell_lev)
    ! purpose: use only local variables for OMP
    !          read leaf cell as in amr2map
    !$ use OMP_LIB
    implicit none 

    integer(kind=4),intent(in)              :: snapnum,icpu
    character(1000),intent(in)              :: repository
    integer(kind=4),intent(out)             :: ncell,nleaf
    real(kind=8),allocatable,intent(out)    :: cell_pos(:,:),cell_var(:,:)
    integer(kind=4),allocatable,intent(out) :: cell_lev(:)

    character(1000)                         :: filename 
    logical                                 :: ok,ok_cell
    integer(kind=4)                         :: i,j,ilevel,ileaf,nx,ny,nz,nlevelmax,nboundary
    integer(kind=4)                         :: idim,ind,iu1,iu2,iu3,rank
    ! stuff read from AMR files
    integer(kind=4)                         :: ncoarse,ngridmax,ngrid_current
    real(kind=8),allocatable                :: xg(:,:)        ! grids position
    integer(kind=4),allocatable             :: son(:,:)       ! sons grids
    real(KIND=8),dimension(1:3)             :: xbound=(/0d0,0d0,0d0/)  
    integer(kind=4),allocatable             :: ngridfile(:,:),ngridlevel(:,:),ngridbound(:,:)
    integer(kind=4)                         :: ngrida
    logical,allocatable                     :: ref(:,:)
    real(kind=8)                            :: dx,boxlen
    integer(kind=4)                         :: ix,iy,iz,ivar,nvarH,nvarRT
    real(kind=8),allocatable                :: xc(:,:),xp(:,:,:)
    ! stuff read from the HYDRO files
    real(kind=8),allocatable                :: var(:,:,:)
    real(kind=4),allocatable                :: var_sp(:)
    
    rank = 1
    !$ rank = OMP_GET_THREAD_NUM()
    iu1 = 10+rank*3
    iu2 = 10+rank*3+1
    iu3 = 10+rank*3+2
    
    ! verify AMR input file
    write(filename,'(a,a,i5.5,a,i5.5,a,i5.5)') trim(repository),'/output_',snapnum,'/amr_',snapnum,'.out',icpu
    inquire(file=filename, exist=ok)
    if(.not. ok)then
       write(*,*)'File '//TRIM(filename)//' not found'    
       stop
    end if

    ! Open AMR file and skip header
    open(unit=iu1,file=filename,form='unformatted',status='old',action='read')
    read(iu1)ncpu
    read(iu1)      !ndim
    read(iu1)nx,ny,nz
    read(iu1)nlevelmax
    read(iu1)ngridmax
    read(iu1)nboundary
    read(iu1)ngrid_current
    read(iu1)boxlen
    do i=1,13
       read(iu1)
    end do
    !twotondim=2**ndim
    xbound=(/dble(nx/2),dble(ny/2),dble(nz/2)/)
    allocate(ngridfile(1:ncpu+nboundary,1:nlevelmax))
    allocate(ngridlevel(1:ncpu,1:nlevelmax))
    if(nboundary>0)allocate(ngridbound(1:nboundary,1:nlevelmax))

    ! Read grid numbers
    read(iu1)ngridlevel
    ngridfile(1:ncpu,1:nlevelmax)=ngridlevel
    read(iu1)
    if(nboundary>0)then
       do i=1,2
          read(iu1)
       end do
       read(iu1)ngridbound
       ngridfile(ncpu+1:ncpu+nboundary,1:nlevelmax)=ngridbound
    endif
    read(iu1)
    ! ROM: comment the single follwing line for old stuff
    read(iu1)
    read(iu1)
    read(iu1)
    read(iu1)
    read(iu1)

    allocate(xc(1:twotondim,1:ndim))


    ! open hydro file and get nvarH
    write(filename,'(a,a,i5.5,a,i5.5,a,i5.5)') trim(repository),'/output_',snapnum,'/hydro_',snapnum,'.out',icpu
    open(unit=iu2,file=filename,form='unformatted',status='old',action='read')
    read(iu2)
    read(iu2)nvarH
    read(iu2)
    read(iu2)
    read(iu2)
    read(iu2)
    
    if (read_rt_variables) then
       ! Open RT file and get nvarRT
       write(filename,'(a,a,i5.5,a,i5.5,a,i5.5)') trim(repository),'/output_',snapnum,'/rt_',snapnum,'.out',icpu
       open(unit=iu3,file=filename,status='old',form='unformatted')
       read(iu3)
       read(iu3)nvarRT
       read(iu3)
       read(iu3)
       read(iu3)
       read(iu3)
    else
       nvarRT = 0
    end if
    
    ncoarse = nx*ny*nz
    ncell   = ncoarse+twotondim*ngridmax

    !allocate(cell_var(1:ncell,1:nvarH+nvarRT))
    allocate(cell_var(1:nvarH+nvarRT,1:ncell))
    allocate(cell_pos(1:ncell,1:3))
    allocate(cell_lev(1:ncell))

    cell_pos = 0.0d0
    cell_var = 0.0d0
    cell_lev = -1
    ileaf=1
    
    ! Loop over levels
    do ilevel=1,nlevelmax
       
       ! Geometry
       dx=0.5**ilevel
       do ind=1,twotondim
          iz=(ind-1)/4
          iy=(ind-1-4*iz)/2
          ix=(ind-1-2*iy-4*iz)
          xc(ind,1)=(dble(ix)-0.5D0)*dx
          xc(ind,2)=(dble(iy)-0.5D0)*dx
          xc(ind,3)=(dble(iz)-0.5D0)*dx
       end do
       
       ! Allocate work arrays
       if(allocated(xg)) then 
          deallocate(xg,son,var,xp,ref)
       endif
       if(allocated(var_sp)) deallocate(var_sp)
       ngrida=ngridfile(icpu,ilevel)
       if(ngrida>0)then
          allocate(xg(1:ngrida,1:ndim))
          allocate(son(1:ngrida,1:twotondim))
          allocate(var(1:ngrida,1:twotondim,1:nvarh+nvarRT))
          if((read_rt_variables .and. rt_Precision.eq.4) .or. U_precision.eq.4) allocate(var_sp(1:ngrida))
          allocate(xp(1:ngrida,1:twotondim,1:ndim))
          allocate(ref(1:ngrida,1:twotondim))
          ref=.false.
       endif
       
       
       ! Loop over domains
       do j=1,nboundary+ncpu
          
          ! Read AMR data
          if(ngridfile(j,ilevel)>0)then
             read(iu1) ! Skip grid index
             read(iu1) ! Skip next index
             read(iu1) ! Skip prev index
             ! Read grid center
             do idim=1,ndim
                if(j.eq.icpu)then
                   read(iu1)xg(:,idim)
                else
                   read(iu1)
                endif
             end do
             read(iu1) ! Skip father index
             do ind=1,2*ndim
                read(iu1) ! Skip nbor index
             end do
             ! Read son index
             do ind=1,twotondim
                if(j.eq.icpu)then
                   read(iu1)son(:,ind)
                else
                   read(iu1)
                end if
             end do
             ! Skip cpu map
             do ind=1,twotondim
                read(iu1)
             end do
             ! Skip refinement map
             do ind=1,twotondim
                read(iu1)
             end do
          endif

          ! Read HYDRO data
          read(iu2)
          read(iu2)
          if(read_rt_variables)read(iu3)
          if(read_rt_variables)read(iu3)
          if(ngridfile(j,ilevel)>0)then
             ! Read hydro variables
             do ind=1,twotondim
                do ivar=1,nvarh
                   if(j.eq.icpu)then
                      if(U_precision.eq.4) then
                         read(iu2) var_sp(:)
                         var(:,ind,ivar) = var_sp(:)
                      else
                         read(iu2)var(:,ind,ivar)
                      endif
                   else
                      read(iu2)
                   end if
                end do
                do ivar=1,nvarRT
                   if(j.eq.icpu)then
                      if(rt_Precision.eq.4) then
                         read(iu3) var_sp(:)
                         var(:,ind,nvarh+ivar) = var_sp(:)
                      else
                         read(iu3)var(:,ind,nvarh+ivar)
                      endif
                   else
                      read(iu3)
                   end if
                end do
             end do
          end if

       enddo

       ! Get leaf cells and store data
       if(ngrida>0)then
          ! Loop over cells
          do ind=1,twotondim
             ! Compute cell center
             do i=1,ngrida
                xp(i,ind,1)=(xg(i,1)+xc(ind,1)-xbound(1))
                xp(i,ind,2)=(xg(i,2)+xc(ind,2)-xbound(2))
                xp(i,ind,3)=(xg(i,3)+xc(ind,3)-xbound(3))
             enddo
             ! Check if cell is refined
             do i=1,ngrida
                ref(i,ind)=son(i,ind)>0.and.ilevel<nlevelmax
             enddo
             ! Store leaf cells
             do i=1,ngrida
                ok_cell= .not.ref(i,ind)
                if(ok_cell)then
                !if(.not.ref(i,ind))then
                   cell_pos(ileaf,1:3) = xp(i,ind,1:3)
                   cell_lev(ileaf) = ilevel
                   !cell_var(ileaf,:) = var(i,ind,:)
                   do ivar=1,nvarh+nvarRT
                      cell_var(ivar,ileaf) = var(i,ind,ivar)
                   enddo
                   ileaf=ileaf+1
                endif
             end do
          end do
       endif
       
       
    enddo

    !!print*,'Number of leaf cells =',ileaf-1,nvarH+nvarRT
    !!print*,icpu,ileaf-1
    nleaf = ileaf-1
    
    close(iu1)
    close(iu2)
    close(iu3)

    return
  end subroutine get_leaf_cells_per_cpu



  function get_nleaf_new(repository,snapnum,ncpu_read,cpu_list)
    ! get the total number of leaf cells
    implicit none
    integer(kind=4),intent(in)                          :: snapnum
    character(1000),intent(in)                          :: repository
    integer(kind=4),intent(in)                          :: ncpu_read
    integer(kind=4),dimension(:),allocatable,intent(in) :: cpu_list
    integer(kind=4)                                     :: get_nleaf_new
    integer(kind=4)                                     :: icpu,k,nleaf,nleaftot

    get_nleaf_new = 0
    nleaftot = 0
    
!$OMP PARALLEL &
!$OMP DEFAULT(private) &
!$OMP SHARED(nleaftot, ncpu_read, cpu_list, repository, snapnum)
!$OMP DO SCHEDULE(DYNAMIC, 10) 
    do k=1,ncpu_read
       icpu=cpu_list(k)
       call get_nleaf_per_cpu(repository,snapnum,icpu,nleaf)
!$OMP ATOMIC
       nleaftot = nleaftot + nleaf
    end do
!$OMP END DO
!$OMP END PARALLEL
    
    get_nleaf_new = nleaftot
    
    return
  end function get_nleaf_new



  subroutine get_nleaf_per_cpu(repository,snapnum,icpu,nleaf)
    ! purpose: get the number of leaf cells in this icpu file

    !$ use OMP_LIB
    implicit none 

    integer(kind=4),intent(in)  :: snapnum,icpu
    character(1000),intent(in)  :: repository
    integer(kind=4),intent(out) :: nleaf

    character(1000)             :: filename 
    logical                     :: ok
    integer(kind=4)             :: i,j,ilevel,nx,ny,nz,nlevelmax,nboundary
    integer(kind=4)             :: idim,ind,iu1,rank
    integer(kind=4)             :: ngridmax,ngrid_current,ngrida
    integer,allocatable         :: son(:,:)       ! sons grids
    real(KIND=8),dimension(1:3) :: xbound=(/0d0,0d0,0d0/)  
    integer, allocatable        :: ngridfile(:,:),ngridlevel(:,:),ngridbound(:,:)
    logical, allocatable        :: ref(:,:)
    real(kind=8)                :: boxlen
    
    rank = 1
    !$ rank = OMP_GET_THREAD_NUM()
    iu1 = 10+rank*3
    
    ! verify AMR input file
    write(filename,'(a,a,i5.5,a,i5.5,a,i5.5)') trim(repository),'/output_',snapnum,'/amr_',snapnum,'.out',icpu
    inquire(file=filename, exist=ok)
    if(.not. ok)then
       write(*,*)'File '//TRIM(filename)//' not found'    
       stop
    end if
    
    ! Open AMR file and skip header
    open(unit=iu1,file=filename,form='unformatted',status='old',action='read')
    read(iu1)ncpu
    read(iu1)      !ndim
    read(iu1)nx,ny,nz
    read(iu1)nlevelmax
    read(iu1)ngridmax
    read(iu1)nboundary
    read(iu1)ngrid_current
    read(iu1)boxlen
    do i=1,13
       read(iu1)
    end do
    !twotondim=2**ndim
    xbound=(/dble(nx/2),dble(ny/2),dble(nz/2)/)
    allocate(ngridfile(1:ncpu+nboundary,1:nlevelmax))
    allocate(ngridlevel(1:ncpu,1:nlevelmax))
    if(nboundary>0)allocate(ngridbound(1:nboundary,1:nlevelmax))
    
    ! Read grid numbers
    read(iu1)ngridlevel
    ngridfile(1:ncpu,1:nlevelmax)=ngridlevel
    read(iu1)
    if(nboundary>0)then
       do i=1,2
          read(iu1)
       end do
       read(iu1)ngridbound
       ngridfile(ncpu+1:ncpu+nboundary,1:nlevelmax)=ngridbound
    endif
    read(iu1)
    ! ROM: comment the single follwing line for old stuff
    read(iu1)
    read(iu1)
    read(iu1)
    read(iu1)
    read(iu1)
    
    nleaf=0
    
    ! Loop over levels
    do ilevel=1,nlevelmax
       
       ! Allocate work arrays
       if(allocated(son)) then 
          deallocate(son,ref)
       endif
       ngrida=ngridfile(icpu,ilevel)
       if(ngrida>0)then
          allocate(son(1:ngrida,1:twotondim))
          allocate(ref(1:ngrida,1:twotondim))
          ref=.false.
       endif
       
       ! Loop over domains
       do j=1,nboundary+ncpu
          ! Read AMR data
          if(ngridfile(j,ilevel)>0)then
             read(iu1) ! Skip grid index
             read(iu1) ! Skip next index
             read(iu1) ! Skip prev index
             do idim=1,ndim ! Skip grid center
                read(iu1)
             end do
             read(iu1) ! Skip father index
             do ind=1,2*ndim
                read(iu1) ! Skip nbor index
             end do
             ! Read son index
             do ind=1,twotondim
                if(j.eq.icpu)then
                   read(iu1)son(:,ind)
                else
                   read(iu1)
                end if
             end do
             ! Skip cpu map
             do ind=1,twotondim
                read(iu1)
             end do
             ! Skip refinement map
             do ind=1,twotondim
                read(iu1)
             end do
          endif
       enddo

       ! Count leaf cells
       if(ngrida>0)then
          ! Loop over cells
          do ind=1,twotondim
             ! Check if cell is refined
             do i=1,ngrida
                ref(i,ind)=son(i,ind)>0.and.ilevel<nlevelmax
                if (.not.ref(i,ind))then
                   nleaf=nleaf+1
                endif
             end do
          end do
       endif
       
    enddo

    close(iu1) 
    return
  end subroutine get_nleaf_per_cpu



  ! JB--
  function get_ncell(repository,snapnum)
    !$ use OMP_LIB
    implicit none
    integer(kind=4),intent(in)  :: snapnum
    character(1000),intent(in)  :: repository
    character(1000)             :: filename
    integer(kind=4)             :: get_ncell
    integer(kind=4)             :: icpu,iunit,rank,nx,ny,nz,ngridmax_l

    ! OMP-safe read ... 
    rank = 1
    !$ rank = OMP_GET_THREAD_NUM()
    iunit = 100 + rank
    
    icpu = 1
    write(filename,'(a,a,i5.5,a,i5.5,a,i5.5)') trim(repository),'/output_',snapnum,'/amr_',snapnum,'.out',icpu
    open(unit=iunit,file=filename,form='unformatted',status='old',action='read')
    read(iunit)
    read(iunit)
    read(iunit) nx,ny,nz
    read(iunit)
    read(iunit) ngridmax_l
    close(iunit)
    get_ncell = nx*ny*nz + twotondim*ngridmax_l
    
    return
  end function get_ncell
  ! --JB
  
  function get_nvar(repository,snapnum)
    implicit none 
    integer(kind=4),intent(in)  :: snapnum
    character(1000),intent(in)  :: repository
    character(1000)             :: filename
    integer(kind=4)             :: get_nvar,icpu,nvarRT

    icpu = 1
    write(filename,'(a,a,i5.5,a,i5.5,a,i5.5)') trim(repository),'/output_',snapnum,'/hydro_',snapnum,'.out',icpu
    open(unit=10,file=filename,form='unformatted',status='old',action='read')
    read(10)
    read(10)get_nvar
    close(10)
    if (read_rt_variables) then
       ! Open RT file and get nvarRT
       write(filename,'(a,a,i5.5,a,i5.5,a,i5.5)') trim(repository),'/output_',snapnum,'/rt_',snapnum,'.out',icpu
       open(unit=12,file=filename,status='old',form='unformatted')
       read(12)
       read(12)nvarRT
       get_nvar = get_nvar + nvarRT
       close(12)
    end if

    return
  end function get_nvar


  function get_ncpu(repository,snapnum)

    implicit none

    integer(kind=4)            :: get_ncpu
    character(512),intent(in)  :: repository
    integer(kind=4),intent(in) :: snapnum

    get_ncpu = nint(get_param_real(repository,snapnum,'ncpu'))

    return
  end function get_ncpu



  function get_param_real(repository,snapnum,param,default_value,rt_info)

    implicit none 

    real(kind=8)                    :: get_param_real
    character(512),intent(in)       :: repository
    integer(kind=4),intent(in)      :: snapnum
    character(*),intent(in)         :: param
    real(kind=8),optional,intent(in):: default_value
    logical,optional,intent(in)     :: rt_info
    logical(kind=4)                 :: not_ok
    character(512)                  :: filename, infofile
    character(512)                  :: line,name,value
    integer(kind=4)                 :: i
    integer(kind=4),parameter       :: param_unit = 13

    not_ok = .true.
    infofile='/info_'
    if(present(rt_info)) then ! read info_rt file
       if(rt_info)  infofile='/info_rt_'
    endif
    write(filename,'(a,a,i5.5,a,i5.5,a)') trim(repository),'/output_',snapnum,trim(infofile),snapnum,'.txt'
    open(unit=param_unit,file=filename,status='old',form='formatted')
    do 
       read(param_unit,'(a)',end=2) line
       i = scan(line,'=')
       if (i==0 .or. line(1:1)=='#') cycle
       name=trim(adjustl(line(:i-1)))
       value=trim(adjustl(line(i+1:)))
       ! check for a comment at end of line !
       i = scan(value,'!')
       if (i /= 0) value = trim(adjustl(value(:i-1)))

       if (trim(name) .eq. trim(param)) then 
          read(value,*) get_param_real
          not_ok = .false.
          exit
       end if

    end do
2   close (param_unit)  

    if (not_ok) then 
       !write(6,*) '--> parameter not found in infoxxx.txt : ',trim(param)
       if(present(default_value)) then
          get_param_real = default_value
       else
          stop
       endif
    end if

    return

  end function get_param_real

  
  subroutine read_conversion_scales(repository,snapnum)

    implicit none 

    character(512),intent(in)  :: repository
    integer(kind=4),intent(in) :: snapnum

    ! set global variables 
    dp_scale_l    = get_param_real(repository,snapnum,'unit_l')
    dp_scale_d    = get_param_real(repository,snapnum,'unit_d')
    dp_scale_t    = get_param_real(repository,snapnum,'unit_t')
    dp_scale_nH   = XH/mp * dp_scale_d      ! convert mass density (code units) to numerical density of H atoms [/cm3]
    dp_scale_nHe  = (1d0-XH)/4.0d0/mp * dp_scale_d ! mass dens to He/cm3
    dp_scale_v    = dp_scale_l/dp_scale_t   ! -> converts velocities into cm/s
    dp_scale_T2   = mp/kB * dp_scale_v**2   ! -> converts P/rho to T/mu, in K
    dp_scale_zsun = 1.d0/0.0127     
    dp_scale_m    = dp_scale_d * dp_scale_l**3 ! convert mass in code units to cgs. 

    return

  end subroutine read_conversion_scales

  subroutine read_cosmo_params(repository,snapnum,omega_0,lambda_0,little_h)

    implicit none 

    character(512),intent(in)  :: repository
    integer(kind=4),intent(in) :: snapnum
    real(kind=8),intent(out)   :: omega_0,lambda_0,little_h

    omega_0  = get_param_real(repository,snapnum,'omega_m')
    lambda_0 = get_param_real(repository,snapnum,'omega_l')
    little_h = get_param_real(repository,snapnum,'H0') / 100.0d0

    return

  end subroutine read_cosmo_params

  subroutine get_fields_from_header(dir,ts,nfields)

    implicit none

    character(1000),intent(in)  :: dir
    integer(kind=4),intent(in)  :: ts
    integer(kind=4),intent(out) :: nfields
    character(2000)             :: filename,line
    integer(kind=4) :: i

    write(filename,'(a,a,i5.5,a,i5.5,a)') trim(dir),'/output_',ts,'/header_',ts,'.txt'
    open(unit=50,file=filename,status='old',action='read',form='formatted')
    read(50,*) ! total nb of particles
    read(50,*)
    read(50,*) ! nb of DM particles
    read(50,*)
    read(50,*) ! nb of star particles
    read(50,*)
    read(50,*) ! nb of sinks
    read(50,*)
    read(50,*) ! nb of clouds
    read(50,*)
    read(50,*) ! Field list
    read(50,'(a)') line
    close(50)

    ! parse the Field list ...
    nfields = 0
    do
       i    = scan(line,' ') ! find a blank
       nfields = nfields + 1
       ParticleFields(nfields) = trim(adjustl(line(:i)))
       line = trim(adjustl(line(i:)))
       if (len_trim(line) == 0) exit
    end do

    return

  end subroutine get_fields_from_header

  !*****************************************************************************************************************

  subroutine read_cooling(repository,snapnum)

    implicit none

    character(1000),intent(in) :: repository
    integer(kind=4),intent(in) :: snapnum
    character(1024)            :: filename
    integer(kind=4)            :: n1,n2

    ! initialize cooling variables
    call clear_cooling

    ! read cooling variables from cooling.out file
    write(filename,'(a,a,i5.5,a,i5.5,a)') trim(repository),'/output_',snapnum,"/cooling_",snapnum,".out"
    open(unit=44,file=filename,form='unformatted')
    read(44) n1,n2
    cooling%n11 = n1
    cooling%n22 = n2
    allocate(cooling%nH(n1),cooling%T2(n2))
    allocate(cooling%cool_com(n1,n2),cooling%heat_com(n1,n2),cooling%cool_com_prime(n1,n2),cooling%heat_com_prime(n1,n2))
    allocate(cooling%cool(n1,n2),cooling%heat(n1,n2),cooling%mu(n1,n2))
    allocate(cooling%cool_prime(n1,n2),cooling%heat_prime(n1,n2),cooling%metal_prime(n1,n2))
    allocate(cooling%metal(n1,n2),cooling%spec(n1,n2,6))
    read(44)cooling%nH
    read(44)cooling%T2
    read(44)cooling%cool
    read(44)cooling%heat
    read(44)cooling%cool_com
    read(44)cooling%heat_com
    read(44)cooling%metal
    read(44)cooling%cool_prime
    read(44)cooling%heat_prime
    read(44)cooling%cool_com_prime
    read(44)cooling%heat_com_prime
    read(44)cooling%metal_prime
    read(44)cooling%mu
    read(44)cooling%spec
    close(44)

    ! define useful quantities for interpolation 
    cool_int%n_nh     = n1
    cool_int%nh_start = minval(cooling%nh)
    cool_int%nh_step  = cooling%nh(2) - cooling%nh(1)
    cool_int%n_t2     = n2
    cool_int%t2_start = minval(cooling%t2)
    cool_int%t2_step  = cooling%t2(2) - cooling%t2(1)

    return

  end subroutine read_cooling

  !*****************************************************************************************************************

  subroutine clear_cooling

    implicit none

    if (cooling%n11 > 0 .or. cooling%n22 > 0) then 
       cooling%n11 = 0
       cooling%n22 = 0
       deallocate(cooling%nH,cooling%T2,cooling%metal)
       deallocate(cooling%heat_com,cooling%cool_com,cooling%heat_com_prime,cooling%cool_com_prime)
       deallocate(cooling%cool,cooling%heat,cooling%metal_prime,cooling%cool_prime)
       deallocate(cooling%heat_prime,cooling%mu,cooling%spec)
    end if

    cool_int%n_nh = 0
    cool_int%nh_start = 0.0d0
    cool_int%nh_step  = 0.0d0
    cool_int%n_t2 = 0
    cool_int%t2_start = 0.0d0
    cool_int%t2_step  = 0.0d0

    return

  end subroutine clear_cooling



  !==================================================================================
  ! STARS utilities 


  subroutine ramses_read_stars_in_domain(repository,snapnum,comp_domp,comp_doms, & 
   star_pos_all,star_age_all,star_mass_all,star_vel_all,star_met_all,cpu_list, & 
   ncpu_read,domdir)

    !$ use OMP_LIB

    implicit none
    real(kind=8),dimension(3),intent(in) :: comp_domp   
    real(kind=8),intent(in) :: comp_doms  
    character(1000),intent(in)             :: repository
    character(2000),intent(in)             :: domdir
    integer(kind=4),intent(in)             :: snapnum,ncpu_read
    integer(kind=4),dimension(:),allocatable,intent(in) :: cpu_list
    real(kind=8)              :: star_pos(3),star_vel(3)
    real(kind=8)                           :: star_age,star_mass,star_met
    real(kind=8),allocatable,intent(inout) :: star_pos_all(:,:),star_age_all(:),star_mass_all(:),star_vel_all(:,:),star_met_all(:)
    integer(kind=8)                        :: nstars
    real(kind=8)                           :: omega_0,lambda_0,little_h,omega_k,H0
    real(kind=8)                           :: aexp,stime,time_cu,boxsize
    integer(kind=4)                        :: ncpu,ilast,icpu,i,ifield,nfields,k
    integer(kind=4)                        :: npart ! changed for HR5
    character(1000)                        :: filename
    integer(kind=8),allocatable            :: id(:)!modified for HR5
    integer,allocatable,dimension(:)       ::lev !modified for HR5
    
    real(kind=8),allocatable               :: age(:),m(:),x(:,:),v(:,:),mets(:),imass(:)
    real(kind=8)                           :: temp(3)
    integer(kind=4)                        :: rank, iunit, ilast_all
    logical::ok,ok_part
    character(2000) :: fichier2

    real(KIND=8)::xmin=0,xmax=1,ymin=0,ymax=1,zmin=0,zmax=1

    ! get cosmological parameters to convert conformal time into ages
    call read_cosmo_params(repository,snapnum,omega_0,lambda_0,little_h)
    omega_k = 0.0d0
    h0      = little_h * 100.0d0
    call ct_init_cosmo(omega_0,lambda_0,omega_k,h0)
    ! compute cosmic time of simulation output (Myr)
    aexp  = get_param_real(repository,snapnum,'aexp') ! exp. factor of output
    stime = ct_aexp2time(aexp) ! cosmic time

    ! read units
    if (.not. conversion_scales_are_known) then 
       call read_conversion_scales(repository,snapnum)
       conversion_scales_are_known = .True.
    end if

    if(.not.cosmo)then
       ! read time
       time_cu = get_param_real(repository,snapnum,'time') ! code unit
       write(*,*)'Time simu [Myr] =',time_cu, time_cu*dp_scale_t/(365.*24.*3600.*1d6)
       boxsize = get_param_real(repository,snapnum,'boxlen') !!!* dp_scale_l  ! [ cm ]
       write(*,*)'boxlen =',boxsize
    endif

    ! read stars 
    nstars = get_tot_nstars(repository,snapnum)
    if (nstars.le.0) then
       write(*,*) 'ERROR : no star particles in output '
       stop
    end if
    !print *, 'nstars:',nstars
    
   xmin = comp_domp(1) - comp_doms/2
   xmax = comp_domp(1) + comp_doms/2
   ymin = comp_domp(2) - comp_doms/2
   ymax = comp_domp(2) + comp_doms/2
   zmin = comp_domp(3) - comp_doms/2
   zmax = comp_domp(3) + comp_doms/2
  
    ! get list of particle fields in outputs 
    call get_fields_from_header(repository,snapnum,nfields)
    ncpu  = get_ncpu(repository,snapnum)
    ilast_all = 1

    fichier2 = "stellar_particles.dat"
    open(unit=199, file=trim(domdir)//trim(fichier2))
    
 
    do k=1,ncpu_read
      icpu=cpu_list(k)
      iunit = 10+k*2 
       write(filename,'(a,a,i5.5,a,i5.5,a,i5.5)') trim(repository), '/output_', snapnum, '/part_', snapnum, '.out', icpu
       open(unit=iunit,file=filename,status='old',form='unformatted')
       !write(*,*)'Processing file '//TRIM(filename),icpu
       read(iunit) !ncpu
       read(iunit) !ndim
       read(iunit)npart !npart
       read(iunit) !localseed
       read(iunit) !nstar_tot
       read(iunit) !mstar_tot
       read(iunit) !mstar_lost
       read(iunit) !nsink
       allocate(age(1:npart))
       allocate(x(1:npart,1:ndim),m(npart),imass(npart))
       allocate(id(1:npart),lev(npart))
       allocate(mets(1:npart))
       allocate(v(1:npart,1:ndim))
       do ifield = 1,nfields
          !print *, ifield, trim(ParticleFields(ifield))
          select case(trim(ParticleFields(ifield)))
          case('pos')
             do i = 1,ndim
                read(iunit) x(1:npart,i)
             end do
          case('vel')
             do i = 1,ndim 
                read(iunit) v(1:npart,i)
             end do
          case('mass')
             read(iunit) m(1:npart)
          case('iord') 
             read(iunit) id(1:npart)
          case('level')
             read(iunit) lev(1:npart)
          case('phi')
              read(iunit) ! changed for HR5 (skipping particle potential) 
          case('tform')
             read(iunit) age(1:npart)
          case('metal')
             read(iunit) mets(1:npart)
          case('propertime')
             read(iunit) ! changed for HR5 (tform_prop)
          case('initialmass')
             read(iunit) imass(1:npart)
          case default
             read(iunit)
             print*,'Error, Field unknown: ',trim(ParticleFields(ifield))
          end select
       end do
       close(iunit)

       if(.not.cosmo)then
          x=x/boxsize
       endif
       !allocate(,star_mass(npart),star_vel(3,npart),star_met(npart))
       ! save star particles within selection region
       ilast = 0
          
       do i=1,npart
               ok_part=(x(i,1)>=xmin.and.x(i,1)<=xmax.and. &
                    &   x(i,2)>=ymin.and.x(i,2)<=ymax.and. &
                    &   x(i,3)>=zmin.and.x(i,3)<=zmax)
       
               ok_part=ok_part.and.((age(i).lt.0.0d0).and.(id(i).gt.0))
               if(ok_part)then
                  !print*,"I am ok",xmin,xmax,ymin,ymax,zmin,zmax
                  ilast = ilast + 1
                if(cosmo)then
                   if (use_proper_time) then
                      star_age = (stime - ct_proptime2time(age(i),h0))*1.d-6 ! Myr
                   else
                      ! Convert from conformal time to age in Myr
                      star_age = (stime - ct_conftime2time(age(i)))*1.d-6 ! Myr
                   end if
                else
                   ! convert from tborn to age in Myr
                   star_age   = max(0.d0, (time_cu - age(i)) * dp_scale_t / (365.d0*24.d0*3600.d0*1.d6))
                endif
                if (use_initial_mass) then 
                   star_mass = imass(i) * dp_scale_m ! [g]
                else
                   star_mass = m(i)     * dp_scale_m ! [g]
                end if
                star_pos = x(i,:)              ! [code units]
                star_vel = v(i,:) * dp_scale_v ! [cm/s]
                star_met = mets(i) 

                write(199,'(I12,4X,6(1PE15.8,4X),I2,4X)')id(i),star_mass, & 
                star_pos(1),star_pos(2),star_pos(3), &
                & star_age,star_met,lev(i)

                
       
               end if
      end do
      !print*,"Ilast final",ilast,npart
      
       ilast_all = ilast_all + ilast

    deallocate(x,v,mets,imass,lev,id,age,m)

    enddo


    return
  end subroutine ramses_read_stars_in_domain


  function get_tot_nstars(dir,ts)

    implicit none 

    integer(kind=4),intent(in) :: ts
    character(1000),intent(in) :: dir
    character(2000)            :: filename
    integer(kind=8)            :: get_tot_nstars ! change for HR5

    get_tot_nstars = 0
    write(filename,'(a,a,i5.5,a,i5.5,a)') trim(dir),'/output_',ts,'/header_',ts,'.txt'
    open(unit=50,file=filename,status='old',action='read',form='formatted')
    read(50,*) ! total nb of particles
    read(50,*)
    read(50,*) ! nb of DM particles
    read(50,*)
    read(50,*) ! nb of star particles 
    read(50,*) get_tot_nstars
    close(50)

    return

  end function get_tot_nstars

  ! conformal time utils :
    function ct_conftime2time(tau)

    ! return look-back time in yr
    
    implicit none 
    
    real(kind=8),intent(in) :: tau
    real(kind=8)            :: ct_conftime2time
    integer(kind=4)         :: i
    
    
    ! locate bracketing conf. times
    i = 1
    do while(tau_frw(i) > tau .and. i < n_frw)
       i = i + 1
    end do
    ! Interploate time
    ct_conftime2time = t_frw(i) * (tau-tau_frw(i-1))/(tau_frw(i)-tau_frw(i-1))+ &
         & t_frw(i-1)        * (tau-tau_frw(i))/(tau_frw(i-1)-tau_frw(i))
    
    return

  end function ct_conftime2time
  
  function ct_proptime2time(tau,h0)
    ! return look-back time in yr
    
    implicit none 
    real(kind=8),intent(in) :: tau,h0
    real(kind=8)            :: ct_proptime2time
    
    ct_proptime2time = tau / (h0 / 3.08d19) / (365.25*24.*3600.)
    return
  end function ct_proptime2time
  
  
  function ct_aexp2time(aexp)

    ! return look-back time in yr

    implicit none

    real(kind=8),intent(in) :: aexp
    real(kind=8)            :: ct_aexp2time
    integer(kind=4)         :: i

    ! find bracketting aexp's 
     i = 1
     do while(aexp_frw(i)>aexp.and.i<n_frw)
        i = i + 1
     end do
     ! Interploate time
     ct_aexp2time = t_frw(i) * (aexp-aexp_frw(i-1))/(aexp_frw(i)-aexp_frw(i-1))+ &
          & t_frw(i-1)    * (aexp-aexp_frw(i))/(aexp_frw(i-1)-aexp_frw(i))

    return
    
  end function ct_aexp2time

  
  subroutine ct_init_cosmo(omega_m,omega_l,omega_k,h0)
    
    ! h0 is in km/s/Mpc

    implicit none 
    real(kind=8),intent(in) :: omega_m,omega_l,omega_k,h0
    real(kind=8)            :: time_tot

    allocate(aexp_frw(0:n_frw),hexp_frw(0:n_frw))
    allocate(tau_frw(0:n_frw),t_frw(0:n_frw))
    call ct_friedman(omega_m,omega_l,omega_k,1.d-6,1.d-3,aexp_frw,hexp_frw,tau_frw,t_frw,n_frw,time_tot)
    ! convert time to yr
    t_frw = t_frw / (h0 / 3.08d19) / (365.25*24.*3600.)

    return
    
  end subroutine ct_init_cosmo


  subroutine ct_clear_cosmo
    
    implicit none
    
    deallocate(aexp_frw,hexp_frw,tau_frw,t_frw)

    return

  end subroutine ct_clear_cosmo


  subroutine ct_friedman(O_mat_0,O_vac_0,O_k_0,alpha,axp_min, &
       & axp_out,hexp_out,tau_out,t_out,ntable,age_tot)

    implicit none
    integer::ntable
    real(kind=8)::O_mat_0, O_vac_0, O_k_0
    real(kind=8)::alpha,axp_min,age_tot
    real(kind=8),dimension(0:ntable)::axp_out,hexp_out,tau_out,t_out
    ! ######################################################!
    ! This subroutine assumes that axp = 1 at z = 0 (today) !
    ! and that t and tau = 0 at z = 0 (today).              !
    ! axp is the expansion factor, hexp the Hubble constant !
    ! defined as hexp=1/axp*daxp/dtau, tau the conformal    !
    ! time, and t the look-back time, both in unit of 1/H0. !
    ! alpha is the required accuracy and axp_min is the     !
    ! starting expansion factor of the look-up table.       !
    ! ntable is the required size of the look-up table.     !
    ! ######################################################!
    real(kind=8)::axp_tau, axp_t
    real(kind=8)::axp_tau_pre, axp_t_pre
    real(kind=8)::dtau,dt
    real(kind=8)::tau,t
    integer::nstep,nout,nskip

    !  if( (O_mat_0+O_vac_0+O_k_0) .ne. 1.0D0 )then
    !     write(*,*)'Error: non-physical cosmological constants'
    !     write(*,*)'O_mat_0,O_vac_0,O_k_0=',O_mat_0,O_vac_0,O_k_0
    !     write(*,*)'The sum must be equal to 1.0, but '
    !     write(*,*)'O_mat_0+O_vac_0+O_k_0=',O_mat_0+O_vac_0+O_k_0
    !     stop
    !  end if

    axp_tau = 1.0D0
    axp_t = 1.0D0
    tau = 0.0D0
    t = 0.0D0
    nstep = 0

    do while ( (axp_tau .ge. axp_min) .or. (axp_t .ge. axp_min) ) 

       nstep = nstep + 1
       dtau = alpha * axp_tau / dadtau(axp_tau,O_mat_0,O_vac_0,O_k_0)
       axp_tau_pre = axp_tau - dadtau(axp_tau,O_mat_0,O_vac_0,O_k_0)*dtau/2.d0
       axp_tau = axp_tau - dadtau(axp_tau_pre,O_mat_0,O_vac_0,O_k_0)*dtau
       tau = tau - dtau

       dt = alpha * axp_t / dadt(axp_t,O_mat_0,O_vac_0,O_k_0)
       axp_t_pre = axp_t - dadt(axp_t,O_mat_0,O_vac_0,O_k_0)*dt/2.d0
       axp_t = axp_t - dadt(axp_t_pre,O_mat_0,O_vac_0,O_k_0)*dt
       t = t - dt

    end do

    age_tot=-t
!!$    write(*,666)-t
!!$666 format(' Age of the Universe (in unit of 1/H0)=',1pe10.3)

    nskip=nstep/ntable

    axp_t = 1.d0
    t = 0.d0
    axp_tau = 1.d0
    tau = 0.d0
    nstep = 0
    nout=0
    t_out(nout)=t
    tau_out(nout)=tau
    axp_out(nout)=axp_tau
    hexp_out(nout)=dadtau(axp_tau,O_mat_0,O_vac_0,O_k_0)/axp_tau

    do while ( (axp_tau .ge. axp_min) .or. (axp_t .ge. axp_min) ) 

       nstep = nstep + 1
       dtau = alpha * axp_tau / dadtau(axp_tau,O_mat_0,O_vac_0,O_k_0)
       axp_tau_pre = axp_tau - dadtau(axp_tau,O_mat_0,O_vac_0,O_k_0)*dtau/2.d0
       axp_tau = axp_tau - dadtau(axp_tau_pre,O_mat_0,O_vac_0,O_k_0)*dtau
       tau = tau - dtau

       dt = alpha * axp_t / dadt(axp_t,O_mat_0,O_vac_0,O_k_0)
       axp_t_pre = axp_t - dadt(axp_t,O_mat_0,O_vac_0,O_k_0)*dt/2.d0
       axp_t = axp_t - dadt(axp_t_pre,O_mat_0,O_vac_0,O_k_0)*dt
       t = t - dt

       if(mod(nstep,nskip)==0)then
          nout=nout+1
          t_out(nout)=t
          tau_out(nout)=tau
          axp_out(nout)=axp_tau
          hexp_out(nout)=dadtau(axp_tau,O_mat_0,O_vac_0,O_k_0)/axp_tau
       end if

    end do
    t_out(ntable)=t
    tau_out(ntable)=tau
    axp_out(ntable)=axp_tau
    hexp_out(ntable)=dadtau(axp_tau,O_mat_0,O_vac_0,O_k_0)/axp_tau

  contains
    function dadtau(axp_tau,O_mat_0,O_vac_0,O_k_0) 
      implicit none 
      real(kind=8)::dadtau,axp_tau,O_mat_0,O_vac_0,O_k_0
      dadtau = axp_tau*axp_tau*axp_tau *  &
           &   ( O_mat_0 + &
           &     O_vac_0 * axp_tau*axp_tau*axp_tau + &
           &     O_k_0   * axp_tau )
      dadtau = sqrt(dadtau)
      return
    end function dadtau
    
    function dadt(axp_t,O_mat_0,O_vac_0,O_k_0)
      implicit none
      real(kind=8)::dadt,axp_t,O_mat_0,O_vac_0,O_k_0
      dadt   = (1.0D0/axp_t)* &
           &   ( O_mat_0 + &
           &     O_vac_0 * axp_t*axp_t*axp_t + &
           &     O_k_0   * axp_t )
      dadt = sqrt(dadt)
      return
    end function dadt
    
  end subroutine ct_friedman


  subroutine read_ramses_params(pfile)
    
    ! ---------------------------------------------------------------------------------
    ! subroutine which reads parameters of current module in the parameter file pfile
    !
    ! default parameter values are set at declaration (head of module)
    ! ---------------------------------------------------------------------------------

    character(*),intent(in) :: pfile
    character(1000) :: line,name,value
    integer(kind=4) :: err,i
    logical         :: section_present

    section_present = .false.
    open(unit=10,file=trim(pfile),status='old',form='formatted')
    ! search for section start
    do
       read (10,'(a)',iostat=err) line
       if(err/=0) exit
       if (line(1:8) == '[ramses]') then
          section_present = .true.
          exit
       end if
    end do
    ! read section if present
    if (section_present) then 
       do
          read (10,'(a)',iostat=err) line
          if(err/=0) exit
          if (line(1:1) == '[') exit ! next section starting... -> leave
          i = scan(line,'=')
          if (i==0 .or. line(1:1)=='#' .or. line(1:1)=='!') cycle  ! skip blank or commented lines
          name=trim(adjustl(line(:i-1)))
          value=trim(adjustl(line(i+1:)))
          i = scan(value,'!')
          if (i /= 0) value = trim(adjustl(value(:i-1)))
          select case (trim(name))
          case ('self_shielding')
             read(value,*) self_shielding
          case ('ramses_rt')
             read(value,*) ramses_rt
          case ('read_rt_variables')
             read(value,*) read_rt_variables
          case ('verbose')
             read(value,*) verbose
          case ('use_initial_mass')
             read(value,*) use_initial_mass
          case ('cosmo')
             read(value,*) cosmo
          case ('use_proper_time')
             read(value,*) use_proper_time
          case('itemp') ! index of thermal pressure
             read(value,*) itemp
          case('imetal')! index of metallicity  
             read(value,*) imetal
          case('ihii') ! index of HII fraction 
             read(value,*) ihii
          case ('iheii') ! index of HeII fraction 
             read(value,*) iheii
          case('iheiii') ! index of HeIII fraction 
             read(value,*) iheiii
          end select
       end do
    end if
    close(10)


    return

  end subroutine read_ramses_params


  
  subroutine print_ramses_params(unit)
    
    ! ---------------------------------------------------------------------------------
    ! write parameter values to std output or to an open file if argument unit is
    ! present.
    ! ---------------------------------------------------------------------------------

    integer(kind=4),optional,intent(in) :: unit

    if (present(unit)) then 
       write(unit,'(a,a,a)') '[ramses]'
       write(unit,'(a,L1)') '  self_shielding    = ',self_shielding
       write(unit,'(a,L1)') '  ramses_rt         = ',ramses_rt
       write(unit,'(a,L1)') '  read_rt_variables = ',read_rt_variables
       write(unit,'(a,L1)') '  use_initial_mass  = ',use_initial_mass
       write(unit,'(a,L1)') '  cosmo             = ',cosmo
       write(unit,'(a,L1)') '  use_proper_time   = ',use_proper_time
       write(unit,'(a,L1)') '  verbose           = ',verbose
       write(unit,'(a,i2)') '  itemp             = ', itemp
       write(unit,'(a,i2)') '  imetal            = ', imetal
       write(unit,'(a,i2)') '  ihii              = ', ihii
       write(unit,'(a,i2)') '  iheii             = ', iheii
       write(unit,'(a,i2)') '  iheiii            = ', iheiii
    else
       write(*,'(a,a,a)') '[ramses]'
       write(*,'(a,L1)') '  self_shielding    = ',self_shielding
       write(*,'(a,L1)') '  ramses_rt         = ',ramses_rt
       write(*,'(a,L1)') '  read_rt_variables = ',read_rt_variables
       write(*,'(a,L1)') '  use_initial_mass  = ',use_initial_mass
       write(*,'(a,L1)') '  cosmo             = ',cosmo
       write(*,'(a,L1)') '  use_proper_time   = ',use_proper_time
       write(*,'(a,L1)') '  verbose           = ',verbose
       write(*,'(a,i2)') '  itemp             = ', itemp
       write(*,'(a,i2)') '  imetal            = ', imetal
       write(*,'(a,i2)') '  ihii              = ', ihii
       write(*,'(a,i2)') '  iheii             = ', iheii
       write(*,'(a,i2)') '  iheiii            = ', iheiii
    end if

    return
  end subroutine print_ramses_params


end module module_ramses
!==================================================================================
!==================================================================================
