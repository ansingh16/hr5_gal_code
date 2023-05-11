module module_mesh

  use module_domain
  use module_gas_composition

  implicit none

  type, public :: mesh
     type(domain)                               :: domain
     integer(kind=4)                            :: nCoarse,nOct,nLeaf,nCell
     integer(kind=4),allocatable,dimension(:)   :: octlevel,son,father
     integer(kind=4),allocatable,dimension(:,:) :: nbor
     real(kind=8),dimension(:,:),allocatable    :: xoct
     type(gas),allocatable,dimension(:)         :: gas
  end type mesh

  private

  integer(kind=4), dimension(1:3,1:8),parameter,private :: nbortest = &
       reshape( (/1,3,5, 2,3,5, 1,4,5, 2,4,5, 1,3,6, 2,3,6, 1,4,6, 2,4,6/), (/3,8/) )  

  real(kind=8), dimension(1:3,1:8), parameter :: offset = &
       reshape((/-.5d0,-.5d0,-.5d0, +.5d0,-.5d0,-.5d0, -.5d0,+.5d0,-.5d0, +.5d0,+.5d0,-.5d0, &
                 -.5d0,-.5d0,+.5d0, +.5d0,-.5d0,+.5d0, -.5d0,+.5d0,+.5d0, +.5d0,+.5d0,+.5d0/), (/3,8/))

  ! TODO: change shape of xoct,xleaf,nbor to optimize access xoct(1:noct,1:3) => xoct(1:3,1:noct)

  ! WORKING VARIABLES (internal to this module, could be renamed)
  
  integer(kind=8)                            :: nOct,nCell
  integer(kind=4)                            :: nCoarse,nleaf
  integer(kind=4)                            :: countleaf,countempty

  ! oct-tree data
  integer(kind=4),allocatable,dimension(:)   :: octlevel,son,father
  integer(kind=4),allocatable,dimension(:,:) :: nbor
  real(kind=8),dimension(:,:),allocatable    :: xoct

  real(kind=8),dimension(:,:),allocatable    :: xleaf
  integer(kind=4),dimension(:),allocatable   :: leaflevel

  ! --------------------------------------------------------------------------
  ! user-defined parameters - read from section [mesh] in the parameter file 
  ! --------------------------------------------------------------------------
  logical :: verbose = .true.      ! set to true to display floods of messages ... 
  ! --------------------------------------------------------------------------
  
  public :: mesh_from_leaves, mesh_from_file, mesh_destructor, dump_mesh, whereisphotongoing, dig_in_cell, &
       in_cell_finder, get_cell_corner, read_mesh_params, print_mesh_params
  private :: add_oct_domain, make_nbor_array, xcson, get_nleaflocal, get_ileaflocal, icell2icell

  contains

    !<><><><><><><><><><><><><><><><><><><><><>
    !
    ! mesh constructors
    !------------------
    ! subroutine mesh_from_leaves(nOctMax,dom,nleaves,leaves,xleaves,leaveslevel,mesh)
    ! subroutine mesh_from_file(file,mesh)
    !
    ! public use
    !------------
    ! subroutine mesh_destructor(mesh)
    ! subroutine dump_mesh(mesh,file)
    ! subroutine WhereIsPhotonGoing(mesh,icellold,xpnew,icellnew,flagoutvol)
    ! function dig_in_cell(mesh,xp,icell)
    ! function in_cell_finder(mesh,xp)
    ! function get_cell_corner(xoct, ind, cell_level) 
    !
    ! private use
    !-------------
    ! recursive subroutine add_oct_domain
    ! subroutine make_nbor_array
    ! function xcson
    ! function get_nleaflocal
    ! function get_ileaflocal
    ! function icell2icell
    ! some routines for checking & verbosity
    !
    !<><><><><><><><><><><><><><><><><><><><><>


    !===============================================================================================
    ! mesh constructors
    !===============================================================================================
    
    subroutine mesh_from_leaves(nOctMax,dom,nleaves,leaves,xleaves,leaveslevel,m)

      type(domain),intent(in)                          :: dom
      type(gas),dimension(:),intent(in)                :: leaves
      integer(kind=4),intent(in)                       :: noctmax,nleaves
      real(kind=8),dimension(1:nleaves,1:3),intent(in) :: xleaves
      integer(kind=4),dimension(1:nleaves),intent(in)  :: leaveslevel
      type(mesh),intent(out)                           :: m

      integer(kind=8)                                  :: ilastoct,ifathercell
      integer(kind=4)                                  :: lfather,i
      real(kind=8),dimension(3)                        :: xnew
      integer(kind=4), dimension(:),allocatable        :: ileaf
      integer(kind=4)                                  :: nocttrue

      ! INIT WORKING VARIABLES
      ncoarse = 1 ! by construction
      noct = nOctMax
      ncell = ncoarse + 8*noct
      nleaf = nleaves

      ! allocate & initialize working arrays
      allocate(xleaf(1:nleaf,1:3),leaflevel(1:nleaf))
      xleaf = xleaves
      leaflevel = leaveslevel
      allocate(ileaf(1:nLeaf))
      allocate(son(1:nCell))
      allocate(father(1:nOct))
      allocate(octlevel(1:nOct))
      allocate(xoct(1:nOct,1:3))
      ! allocate(nbor(1:nOct,1:6)) ! will be done later
      son=0
      father=-1
      octlevel=-1
      xoct=-99.
      !nbor=0
      if(size(son).le.0)then !sanity check
         print*,'problem in module_mesh, mesh_from_leaves ', size(son), ncell, nOctMax
         print*,'cannot allocate arrays of this size, decrease nOctMax...'
         stop
      endif
      
      countleaf=0
      countempty=0

      ! first pass in the coarse cell
      ilastoct=0
      xnew=(/0.5d0,0.5d0,0.5d0/)
      lfather=0
      ifathercell=1
      do i=1,nLeaf
         ileaf(i)=i
      enddo
      son(1)=1

      if (verbose) write(*,*)'Recursive loop, work in progress...'

      ! build the mesh
      call add_oct_domain(ilastoct,nLeaf,ileaf,xnew,lfather,ifathercell)

      ! resize arrays
      nOctTrue = int(ilastoct,4)
      call resize_octtree(nOct,nOctTrue)
      nOct  = nOctTrue
      nCell = nCoarse + 8*nOct

      if (verbose) then 
         write(*,*) ' nCoarse = ',nCoarse
         write(*,*) ' nOct    = ',nOct
         write(*,*) ' nLeaf   = ',nLeaf
         write(*,*) ' nCell   = ',nCell
      end if
      
      ! make nbor array
      if (verbose) then 
         !write(*,*)
         write(*,*)'Building nbor array...'
      end if
      allocate(nbor(nOct,6))
      nbor = 0
      call make_nbor_array

      ! work should be done
      ! == son, father, xoct, octlevel should be filled
      if (countleaf /= nleaves) then
         write(*,*)'mesh_from_leaves: number of leaves counted vs. init =',countleaf,nleaves
         stop
      end if
      ! do some checks if you want
      call check_octtree


      ! fill the mesh type structure
      allocate(m%son(ncell),m%father(noct),m%octlevel(noct),m%xoct(noct,3),m%nbor(noct,6))

      ! allocate & fill m%domain%stuff

      m%domain = dom

      m%ncoarse = ncoarse
      m%noct    = int(noct,4)
      m%ncell   = int(ncell,4)
      m%nleaf   = nleaf

      m%son = son
      m%father = father
      m%octlevel = octlevel
      m%xoct = xoct
      m%nbor = nbor

      ! allocate & fill m%gas
      !JB- the allocation seems to be necessary (fixes a crash).
      allocate(m%gas(nleaf))
      m%gas = leaves
      
      
      ! check code
      call check_struct(nleaf,ileaf,m)
      ! clean-up
      ! deallocate arrays
      deallocate(father,son,nbor,octlevel,xoct,ileaf,xleaf,leaflevel)

    end subroutine mesh_from_leaves



    subroutine mesh_from_file(file,m)

      character(2000),intent(in)                 :: file
      type(mesh),intent(out)                     :: m
      integer(kind=4),allocatable,dimension(:)   :: iarr1d
      integer(kind=4),allocatable,dimension(:,:) :: iarr2d
      real(kind=8),dimension(:,:),allocatable    :: farr2d
      type(domain)                               :: dom
      type(gas),dimension(:),allocatable         :: g

      open(unit=13, file=file, status='old', form='unformatted', action='read')

      ! + read type domain
      call read_domain_bin(unit=13,dom=dom)
      m%domain = dom

      read(13) m%ncoarse, m%noct, m%ncell, m%nleaf

      ! allocate m
      allocate(m%son(m%ncell),m%father(m%noct),m%octlevel(m%noct),m%xoct(m%noct,3),m%nbor(m%noct,6))

      ! father
      allocate(iarr1d(m%noct))
      read(13) iarr1d
      m%father = iarr1d
      deallocate(iarr1d)

      ! son
      allocate(iarr1d(m%ncell))
      read(13) iarr1d
      m%son = iarr1d
      deallocate(iarr1d)

      ! nbor
      allocate(iarr2d(m%noct,6))
      read(13) iarr2d
      m%nbor = iarr2d
      deallocate(iarr2d)

      ! octlevel
      allocate(iarr1d(m%noct))
      read(13) iarr1d
      m%octlevel = iarr1d
      deallocate(iarr1d)

      ! xoct
      allocate(farr2d(m%noct,3))
      read(13) farr2d
      m%xoct = farr2d
      deallocate(farr2d)

      ! + read type gas
      call read_gas(13,m%nleaf,g)
      allocate(m%gas(1:m%nleaf))
      m%gas = g

      close(13)
      
    end subroutine mesh_from_file


    !===============================================================================================
    ! public use
    !===============================================================================================

    subroutine whereIsPhotonGoing(m,icellold,xpnew,icellnew,flagoutvol)

      type(mesh),intent(in)                :: m
      integer(kind=4), intent(in)          :: icellold
      real(kind=8),dimension(3),intent(in) :: xpnew
      integer(kind=4), intent(out)         :: icellnew
      logical, intent(out)                 :: flagoutvol
      logical                              :: inside
      integer(kind=4)                      :: ind,ioct,ifather,level,ison,indleaf,ioctleaf,icell,i
      real(kind=8),dimension(3)            :: xc
      real(kind=8)                         :: dx

      ! 1/ test if photon in father(myoct)
      !     if yes -> dig
      !
      ! 2/ if no, test 3 neighbours according to ind in oct
      !
      ! 3/ if no, (photon probably escaped by edges or corners) then find where is photon by digging from the top

      flagoutvol=.false.

      ! find ind and ioct
      ind  = (icellold-m%ncoarse-1)/m%noct + 1
      ioct = icellold - m%ncoarse - (ind-1)*m%noct

      ! test if photon is in father(ioct)
      ifather = m%father(ioct)    ! is a cell
      xc      = m%xoct(ioct,1:3)
      level   = m%octlevel(ioct)  ! level l
      dx      = 0.5d0**(level)    ! cell size in oct, so half size of father cell at level l-1 
      inside  = ((abs(xc(1)-xpnew(1))<dx).and.(abs(xc(2)-xpnew(2))<dx).and.(abs(xc(3)-xpnew(3))<dx))
      
      if(inside)then
         ! if yes -> dig
         icellnew = dig_in_cell(m,xpnew,ifather)
         ! check if icellnew is inside the domain
         if(m%son(icellnew)==0) flagoutvol=.true.
         return

      else ! if no -> test 3 neighbors
         do i=1,3
            icell = m%nbor(ioct,nbortest(i,ind))
            ! The nbor of the oct are the neighboring cells of the father of the oct. 
            ! 3 possible cases: 
            ! 2a/ if son(icell)>0 -> this is an oct 
            !                     -> we have xoct and we can test if photon is inside, if yes, dig 
            ! 2b/ if son(icell)<0 -> this is a leaf cell 
            !                     -> we have ind and ioct, with xoct and offset array we get xleaf and we can test if photon is inside this leaf
            ! 2c/ if son(icell)=0 -> this is a leaf cell 
            !                     -> we have ind and ioct, with xoct and offset array we get xleaf and we can test if photon is inside. 
            !                        If inside, return flag saying that photon is moving to another domain. 
            ! Note for cases b & c: level is l+1 compared to level of case a, then we take 2*dx 

            ison=m%son(icell)
            if(ison>0)then
               xc     = m%xoct(ison,1:3)
               !level  = m%octlevel(ison)    ! = level l as defined above
               !dx     = 0.5d0**(level)
               inside = ((abs(xc(1)-xpnew(1))<dx).and.(abs(xc(2)-xpnew(2))<dx).and.(abs(xc(3)-xpnew(3))<dx))      
               if(inside)then
                  ! if yes -> dig
                  icellnew = dig_in_cell(m,xpnew,icell)
                  ! check if icellnew is inside the domain
                  if(m%son(icellnew)==0) flagoutvol=.true.
                  return
               endif
            else
               indleaf  = (icell-m%ncoarse-1)/m%noct + 1
               ioctleaf = icell - m%ncoarse - (indleaf-1)*m%noct
               xc       = m%xoct(ioctleaf,1:3)
               !level    = m%octlevel(ioctleaf)           ! = level l+1 compared to level defined above
               !dx       = 0.5d0**(level)                 ! we keep dx from level l 
               xc       = xc + offset(1:3,indleaf)*dx*2   ! but then 2*dx for level l+1 
               inside   = ((abs(xc(1)-xpnew(1))<dx).and.(abs(xc(2)-xpnew(2))<dx).and.(abs(xc(3)-xpnew(3))<dx))      
               if(inside)then
                  icellnew=icell
                  if(ison==0)flagoutvol=.true.
                  return
               endif
            endif
         enddo
         ! if not in the 3 neighbors, then dig from the top (escaped by an edge or a corner...)
         icellnew = in_cell_finder(m,xpnew)
         ! test if son(icellnew)==0 -> outside of domain
         if(m%son(icellnew)==0) flagoutvol=.true.
         return
      endif

      print *,'ERROR: where is the photon...'
      stop

    end subroutine whereIsPhotonGoing



    function dig_in_cell(m,xp,icell)

      type(mesh), intent(in)               :: m
      real(kind=8),dimension(3),intent(in) :: xp
      integer(kind=4),intent(in)           :: icell
      integer(kind=4)                      :: dig_in_cell
      integer(kind=4)                      :: ison,ix,iy,iz,ind,icellint
      real(kind=8),dimension(3)            :: x

      ! find the leaf cell containing xp, starting from icell
      ! this assumes that xp is contained inside icell... 

      ison = m%son(icell)               ! its (oct) son
      do while (ison > 0)               ! there is an oct in current cell 
         x  = m%xoct(ison,1:3)          ! oct position
         ix = merge(0,1,xp(1) < x(1)) 
         iy = merge(0,1,xp(2) < x(2))
         iz = merge(0,1,xp(3) < x(3))
         ind = ix + 2*iy + 4*iz 
         icellint = m%ncoarse + ind * m%nOct + ison
         ison  = m%son(icellint)
      end do
      dig_in_cell=icellint

      return
    end function dig_in_cell
    


    function in_cell_finder(m,xp)

      type(mesh),intent(in)                :: m
      real(kind=8),dimension(3),intent(in) :: xp
      integer(kind=4)                      :: in_cell_finder
      integer(kind=4)                      :: ison,ix,iy,iz,ind,icell
      real(kind=8),dimension(3)            :: x

      ! find leaf cell containing point xp
      icell = 1 ! the coarse cell containing the whole box is always nb 1
      ison  = m%son(icell) ! its (oct) son
      do while (ison > 0)  ! there is an oct in current cell 
         x  = m%xoct(ison,1:3)  ! oct position
         ix = merge(0,1,xp(1) < x(1)) 
         iy = merge(0,1,xp(2) < x(2))
         iz = merge(0,1,xp(3) < x(3))
         ind = ix + 2*iy + 4*iz
         icell = m%ncoarse + ind * m%nOct + ison
         ison  = m%son(icell)
      end do
      !in_cell_finder=-ison
      in_cell_finder=icell

      return
    end function in_cell_finder



    function get_cell_corner(xoct, ind, cell_level) 
      ! return position of cell corner, in box units.
      
      real(kind=8),dimension(3),intent(in) :: xoct
      integer(kind=4),intent(in)           :: ind, cell_level
      real(kind=8),dimension(3)            :: get_cell_corner, xcell
      real(kind=8)                         :: dx

      dx              = 0.5d0**cell_level
      xcell           = xoct + offset(1:3,ind)*dx
      get_cell_corner = xcell - dx*0.5d0

      return
    end function get_cell_corner



    subroutine dump_mesh(m,file)

      type(mesh),intent(in)    :: m
      character(2000),intent(in) :: file

      ! dump data in binary format
      if (verbose) write(*,*) 'Writing mesh in file: ',trim(file)

      open(unit=13, file=file, status='unknown', form='unformatted', action='write')

      ! dump domain
      call dump_domain_bin(unit=13,dom=m%domain)

      write(13) m%nCoarse, m%nOct, m%nCell, m%nLeaf
      write(13) m%father(1:nOct)
      write(13) m%son(1:nCell)
      write(13) m%nbor(1:nOct,1:6)
      write(13) m%octlevel(1:nOct)
      write(13) m%xoct(1:nOct,1:3)

      ! dump gas
      call dump_gas(13,m%gas)

      close(13)

    end subroutine dump_mesh



    subroutine mesh_destructor(m)
      type(mesh),intent(inout) :: m
      call gas_destructor(m%gas)
      ! maybe not needed...
      !call domain_destructor(m%domain)
      deallocate(m%father,m%son,m%nbor,m%xoct,m%octlevel)
    end subroutine mesh_destructor


    !===============================================================================================
    ! PRIVATE routines
    !===============================================================================================


    recursive subroutine add_oct_domain(ilastoct,n,id,fcx,fcl,ifc)
    
      integer(kind=8),intent(inout)            :: ilastoct
      integer(kind=4),intent(in)               :: n
      integer(kind=4),dimension(n),intent(in)  :: id
      integer(kind=4),intent(in)               :: fcl
      integer(kind=8),intent(in)               :: ifc
      real(kind=8),dimension(3),intent(in)     :: fcx
      integer(kind=8)                          :: icell,ilastoctLevel
      integer(kind=4)                          :: ind,level,nleaflocal
      integer(kind=4),dimension(:),allocatable :: ileaflocal
      real(kind=8),dimension(3)                :: xcentre
      logical                                  :: same_level

      ilastoct = ilastoct + 1 
      octlevel(ilastoct)=fcl+1
      xoct(ilastoct,1:3) = fcx(1:3)
      father(ilastoct) = int(ifc,4)
      level=octlevel(ilastoct)
      ilastoctLevel = ilastoct

      if(ilastoct .gt. size(father))then
         print*,'problem in module_mesh, add_oct_domain ',ilastoct,size(father)
         print*,'nOctMax used to allocate arrays is probably too small, try to increase it...'
         stop
      endif
      
      do ind=1,8

         if(level==1 .and. verbose)then
            write(*,'(a,i1,a)')'        \_level=1 cell # ',ind,'/8 done'
         endif

         ! define cell centre and dx 
         xcentre = xcson(ind,fcx,level)

         nleaflocal = get_nleaflocal(n,id,xcentre,level)

         icell = nCoarse+ilastoctLevel+(ind-1)*nOct
         if(icell .gt. size(son)) then !JOKI DEBUG!!
            print*,'problem in module_mesh, add_oct_domain ',icell,nCoarse,ilastoctLevel,nOct,ind,nleaflocal,level, size(son)
            stop
         endif

         ! 3 possible cases according to the value of nleaflocal:
         ! 1/ nleaflocal > 1 -> add a new oct
         ! 2/ nleaflocal = 1 -> this is a true leaf OR we have only one of the leaves in the domain...
         ! 3/ nleaflocal = 0 -> no leaf

         if(nleaflocal>1)then

            ! add a new oct
            allocate(ileaflocal(nleaflocal))
            ileaflocal = get_ileaflocal(n,id,xcentre,level,nleaflocal)                   ! indices dans xl(:,1:3)
            son(icell) = int(ilastoct+1,4)
            call add_oct_domain(ilastoct,nleaflocal,ileaflocal,xcentre,level,icell)
            deallocate(ileaflocal)

         else

            if(nleaflocal==1)then

               ! this is the true leaf OR one of the leaves (the other are not present in the selection)
               allocate(ileaflocal(nleaflocal))
               ileaflocal = get_ileaflocal(n,id,xcentre,level,nleaflocal)
               same_level=all(level==leaflevel(ileaflocal))
               if(same_level)then
                  ! true leaf
                  son(icell)=-ileaflocal(1)
                  countleaf = countleaf+1
                  deallocate(ileaflocal)
               else
                  ! add a new oct
                  son(icell) = int(ilastoct+1,4)
                  call add_oct_domain(ilastoct,nleaflocal,ileaflocal,xcentre,level,icell)
                  deallocate(ileaflocal)
               endif

            else
               ! nleaflocal=0
               ! there is no leaf, flag son to 0

               if(nleaflocal/=0)then
                  print*,'ERROR: Oh no, what happens',nleaflocal
                  stop
               endif

               son(icell)=0
               countempty=countempty+1

            endif

         endif

      enddo

    end subroutine add_oct_domain

    
    
    function xcson(ind,fx,level)

      integer(kind=4),intent(in)            :: ind,level
      real(kind=8),dimension(3),intent(in)  :: fx
      real(kind=8)                          :: dx
      real(kind=8),dimension(3)             :: xcson

      !dx = 0.5d0**level / 2.
      dx = 0.5d0**(level+1)
      select case(ind) 
      case(1) ! 1st cell in oct
         xcson(1) = fx(1) - dx
         xcson(2) = fx(2) - dx
         xcson(3) = fx(3) - dx
      case(2) ! 2nd cell in oct
         xcson(1) = fx(1) + dx
         xcson(2) = fx(2) - dx
         xcson(3) = fx(3) - dx
      case(3) ! 3rd cell in oct
         xcson(1) = fx(1) - dx
         xcson(2) = fx(2) + dx
         xcson(3) = fx(3) - dx
      case(4) ! 4th cell in oct
         xcson(1) = fx(1) + dx
         xcson(2) = fx(2) + dx
         xcson(3) = fx(3) - dx
      case(5) ! 5th cell in oct
         xcson(1) = fx(1) - dx
         xcson(2) = fx(2) - dx
         xcson(3) = fx(3) + dx
      case(6) ! 6th cell in oct
         xcson(1) = fx(1) + dx
         xcson(2) = fx(2) - dx
         xcson(3) = fx(3) + dx
      case(7) ! 7th cell in oct
         xcson(1) = fx(1) - dx
         xcson(2) = fx(2) + dx
         xcson(3) = fx(3) + dx
      case(8) ! 8th cell in oct
         xcson(1) = fx(1) + dx
         xcson(2) = fx(2) + dx
         xcson(3) = fx(3) + dx
      end select

      return

    end function xcson


    
    function get_nleaflocal(n,ileaf,xc,level)

      integer(kind=4),intent(in)              :: n,level
      integer(kind=4),dimension(n),intent(in) :: ileaf
      real(kind=8),dimension(3),intent(in)    :: xc
      integer(kind=4)                         :: icell,jcell,kcell,ii,jj,kk,i,count,get_nleaflocal

      count=0
      icell=int(xc(1)*2**level)+1 
      jcell=int(xc(2)*2**level)+1 
      kcell=int(xc(3)*2**level)+1 

      do i=1,n
         ii=int(xleaf(ileaf(i),1)*2**level)+1
         jj=int(xleaf(ileaf(i),2)*2**level)+1
         kk=int(xleaf(ileaf(i),3)*2**level)+1
         if(ii==icell.and.jj==jcell.and.kk==kcell)then
            ! ileaf(i) in cell
            count=count+1
         endif
      enddo
      get_nleaflocal=count

      return

    end function get_nleaflocal



    function get_ileaflocal(n,ileaf,xc,level,nleaflocal)

      integer(kind=4),intent(in)              :: n,level,nleaflocal
      integer(kind=4),dimension(n),intent(in) :: ileaf
      real(kind=8),dimension(3),intent(in)    :: xc
      integer(kind=4)                         :: icell,jcell,kcell,ii,jj,kk,i,count
      integer(kind=4),dimension(nleaflocal)   :: get_ileaflocal

      count=0
      icell=int(xc(1)*2**level)+1 
      jcell=int(xc(2)*2**level)+1 
      kcell=int(xc(3)*2**level)+1 

      do i=1,n
         ii=int(xleaf(ileaf(i),1)*2**level)+1
         jj=int(xleaf(ileaf(i),2)*2**level)+1
         kk=int(xleaf(ileaf(i),3)*2**level)+1
         if(ii==icell.and.jj==jcell.and.kk==kcell)then
            ! ileaf(i) in cell
            count=count+1
            get_ileaflocal(count)=ileaf(i)
         endif
      enddo

      return
    end function get_ileaflocal


    
    function icell2icell(icell,ngridmax,ngridmin,ncoarse)
      ! copyright JB
      ! returns icell indexed with ngridmin instead of ngridmax

      integer(kind=4),intent(in) :: icell,ngridmin,ncoarse
      integer(kind=8),intent(in) :: ngridmax
      integer(kind=4)            :: icell2icell
      integer(kind=4)            :: ind,jcell

      ind   = (icell - ncoarse - 1) / int(ngridmax,4) + 1
      jcell = icell - ncoarse - (ind - 1) * int(ngridmax,4)

      icell2icell = ncoarse + jcell + (ind-1) * ngridmin

      return

    end function icell2icell



    subroutine make_nbor_array
      ! copyright JB
      
      ! need octlevel, xoct, ncoarse, noct, son 
      integer(kind=4) :: inbor,level,ix,iy,iz,icell,l,ison
      integer(kind=4) :: i,ind,ioct
      real(kind=8)    :: xc(3),xnbor(3),x(3),dx
      integer(kind=4) :: noops1

      if (ncoarse /= 1) then 
         Print*,'ERROR: Oh no, ncoarse /= 1... '
         stop
      end if
      noops1=0
      do ioct=1,int(nOct,4)  ! loop over all octs. 
         level = octlevel(ioct)   ! level of oct -> 6 neighbors are cells at level-1
         dx    = 0.5d0**(level-1) ! size of neighbor cells (or cell containing the oct)
         xc    = xoct(ioct,1:3)   ! position of current oct.
         if (level == 1) then     ! this is the oct within the unique coarse cell (at l=0)/
            nbor(ioct,1:6) = 1    ! -> define its neighbors as the coarse cell itself.
            cycle
         end if
         if (level <= 0) then     ! this should not happen, but it does : what are these octs ? 
            if (level == 0) noops1 = noops1 + 1
            if (level < -1) print*,'oops'
            cycle
         end if
         do inbor=1,6   ! loop over 6 neighbor cells and compute the coords. of their center
            select case(inbor) 
            case (1) ! 1st neighbour (-dx)
               xnbor(1) = xc(1) - dx
               xnbor(2) = xc(2)
               xnbor(3) = xc(3)
            case (2) ! 2nd neighbour (+dx)
               xnbor(1) = xc(1) + dx
               xnbor(2) = xc(2)
               xnbor(3) = xc(3)
            case (3) ! 3rd neighbour (-dy)
               xnbor(1) = xc(1)
               xnbor(2) = xc(2) - dx
               xnbor(3) = xc(3)
            case (4) ! 4th neighbour (+dy)
               xnbor(1) = xc(1)
               xnbor(2) = xc(2) + dx
               xnbor(3) = xc(3)
            case (5) ! 5th neighbour (-dz)
               xnbor(1) = xc(1)
               xnbor(2) = xc(2)
               xnbor(3) = xc(3) - dx
            case (6) ! 6th neighbour (+dz)
               xnbor(1) = xc(1)
               xnbor(2) = xc(2)
               xnbor(3) = xc(3) + dx
            end select
            ! correct for periodic boundary conditions :
            do i = 1,3
               if (xnbor(i) > 1.0d0) xnbor(i) = xnbor(i) - 1.0d0
               if (xnbor(i) < 0.0d0) xnbor(i) = xnbor(i) + 1.0d0
            end do
            ! find cell of level "level-1" containing xnbor
            icell = 1 ! the coarse cell containing the whole box is always nb 1
            l     = 0 ! its level 
            ison  = son(icell) ! its (oct) son
            do while (ison > 0)  ! there is an oct in current cell 
               x  = xoct(ison,1:3)  ! oct position
               ix = merge(0,1,xnbor(1) < x(1)) 
               iy = merge(0,1,xnbor(2) < x(2))
               iz = merge(0,1,xnbor(3) < x(3))
               ind = ix + 2*iy + 4*iz 
               icell = ncoarse + ind * int(nOct,4) + ison
               ison  = son(icell)
               l   = l + 1
               if (l == level-1) exit
            end do
            nbor(ioct,inbor) = icell 
            if((ison/=0).and.(l/=level-1))then
               print*,'ERROR: Ouh la la',ioct,ison,l,level,inbor
               stop
            end if
         end do
      end do
      if (noops1 /= 0) print*,'-- make_nbor_array: error flag = ',noops1
      !-JB

    end subroutine make_nbor_array



    subroutine resize_octtree(nold,nnew)

      integer(kind=8), intent(in)              :: nold
      integer(kind=4), intent(in)              :: nnew
      real(kind=8),dimension(:,:),allocatable  :: tmpr
      integer(kind=4),dimension(:),allocatable :: tmpi
      integer(kind=8)                          :: ncellnew,i

      ncellnew = nCoarse+8*nnew

      if (verbose) write(*,*) 'Resizing arrays...'
      
      ! transform & resize father
      allocate(tmpi(nnew))
      do i=1,nnew
         tmpi(i) = icell2icell(father(i),nold,nnew,nCoarse) ! change father cell num from indexation with nOctTot octs to indexation with ilastoct octs.
      end do
      deallocate(father)
      allocate(father(nnew))
      father = tmpi
      deallocate(tmpi)

      ! transform & resize son
      allocate(tmpi(ncellnew))
      tmpi(1)=son(1)
      do i=1,8  ! gather all cells in begining of array (don't change values which are oct indexes).
         tmpi(nCoarse+1+(i-1)*nnew:nCoarse+i*nnew) = son(nCoarse+1+(i-1)*nold:nCoarse+(i-1)*nold+nnew)
      end do
      deallocate(son)
      allocate(son(ncellnew))
      son = tmpi
      deallocate(tmpi)

      ! resize octlevel
      allocate(tmpi(1:nnew))
      tmpi = octlevel(1:nnew)
      deallocate(octlevel)
      allocate(octlevel(nnew))
      octlevel = tmpi
      deallocate(tmpi)
      
      ! resize nbor
      !deallocate(nbor)
      !allocate(nbor(nOctTot,6))
      !nbor = 0

      ! resize xoct
      allocate(tmpr(nnew,3))
      tmpr = xoct(1:nnew,:)
      deallocate(xoct)
      allocate(xoct(nnew,3))
      xoct = tmpr
      deallocate(tmpr)
      
    end subroutine resize_octtree
    


    subroutine check_octtree

      integer(kind=4)            :: oct_status,i,ind,ioct,icount,icount2

      ! check results
      if (verbose) then 
         !write(*,*)
         !write(*,*)'...checking arrays...'
         write(*,*)'min max son          ',minval(son(:)),maxval(son(:))
         write(*,*)'min max father       ',minval(father(:)),maxval(father(:))
         write(*,*)'min max octlevel     ',minval(octlevel(:)),maxval(octlevel(:))
         write(*,'(a,f12.8,f12.8)')' min max xoct         ',minval(xoct(:,1)),maxval(xoct(:,1))
         write(*,'(a,f12.8,f12.8)')' min max yoct         ',minval(xoct(:,2)),maxval(xoct(:,2))
         write(*,'(a,f12.8,f12.8)')' min max zoct         ',minval(xoct(:,3)),maxval(xoct(:,3))
         !write(*,*)'min max octlevel (>0)',minval(octlevel, mask=(octlevel >= 0)), maxval(octlevel, mask=(octlevel >= 0))
         !write(*,*)'min max xoct (>0)    ',minval(xoct(:,1), mask=(xoct(:,1)>=0)),maxval(xoct(:,1), mask=(xoct(:,1)>=0))
         !write(*,*)'min max yoct (>0)    ',minval(xoct(:,2), mask=(xoct(:,2)>=0)),maxval(xoct(:,2), mask=(xoct(:,2)>=0))
         !write(*,*)'min max zoct (>0)    ',minval(xoct(:,3), mask=(xoct(:,3)>=0)),maxval(xoct(:,3), mask=(xoct(:,3)>=0))
         !write(*,*)
      end if

      countleaf=0
      icount=0
      do i=1,nCoarse
         icount=icount+1
         ioct=i
         if(son(ioct)<0)then
            countleaf=countleaf+1
         endif
      enddo

      if (ncoarse /= 1) then
         write(*,*)'check_octtree: problem with ncoarse',ncoarse
         stop
      else
         if (count(mask=(son==1))/=1) then
            write(*,*)'check_octtree: problem with son(ncoarse)',son(1:ncoarse)
            stop
         endif
      endif
      
      do i=1,int(nOct,4)
         do ind=1,8
            icount=icount+1
            ioct = nCoarse + i + (ind-1)*int(nOct,4)
            if(son(ioct)<0)then
               countleaf=countleaf+1
            endif
         enddo
      enddo
      
      if (icount /= ncell) then
         write(*,*)'check_octtree: problem with son array -> # of elements = ',icount,ncell
         stop
      endif
      if (countleaf /= nleaf) then
         write(*,*)'check_octtree: problem with son array -> # of leaves = ',countleaf,nleaf
         stop
      endif
      
      
      ! count incomplete octs...
      icount=0
      icount2=0
      do i=1,int(nOct,4)
         oct_status=0
         do ind=1,8
            ioct = nCoarse + i + (ind-1)*int(nOct,4)
            if(son(ioct)/=0)then
               oct_status=oct_status+1
            endif
         enddo
         if (oct_status/=8)then
            icount=icount+1
         endif
         if (oct_status==0)then
            icount2=icount2+1
         endif
      enddo
      if (icount2 /= 0) then
         write(*,*)'check_octtree: problem with son array -> # of empty oct = ',icount2
         stop
      endif

      ! father array
      if (count(mask=(father==-1)) /= 0) then
         write(*,*)'check_octtree: problem with father array -> # of -1 = ',count(mask=(father==-1))
         stop
      endif
      if (count(mask=(father==0)) /= 0) then
         write(*,*)'check_octtree: problem with father array -> # of  0 = ',count(mask=(father==0))
         stop
      endif
      if (count(mask=(father>0)) /= noct) then
         write(*,*)'check_octtree: problem with father array -> # of >0 = ',count(mask=(father>0)),noct
         stop
      endif
      
      ! if no stop, check ok
      
    end subroutine check_octtree



    subroutine check_struct(n,ileaf,m)
      ! simple test to check the structure
      ! for each leaf cell, take its central position, find in which cell it belongs, and check that son(cell_found)=-indice_cell

      integer(kind=4), intent(in)               :: n
      integer(kind=4),dimension(1:n),intent(in) :: ileaf
      type(mesh),intent(in)                     :: m
      integer(kind=4)                           :: i,ic
      real(kind=8),dimension(3)                 :: xc

      do i=1,n
         xc(1:3) = xleaf(i,1:3)
         ic = in_cell_finder(m,xc)
         ! some check
         if(m%son(ic)/=-ileaf(i))then
            print*,'check_struct: problem with in_cell_finder...',i,ic,ileaf(i),m%son(ic)
            stop
         endif
      end do

    end subroutine check_struct

    

  subroutine read_mesh_params(pfile)
    
    ! ---------------------------------------------------------------------------------
    ! subroutine which reads parameters of current module in the parameter file pfile
    ! default parameter values are set at declaration (head of module)
    !
    ! ALSO read parameter form used modules (gas_composition)
    ! ---------------------------------------------------------------------------------

    character(*),intent(in) :: pfile
    character(1000)         :: line,name,value
    integer(kind=4)         :: err,i
    logical                 :: section_present

    section_present = .false.
    open(unit=10,file=trim(pfile),status='old',form='formatted')
    ! search for section start
    do
       read (10,'(a)',iostat=err) line
       if(err/=0) exit
       if (line(1:6) == '[mesh]') then
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
          case ('verbose')
             read(value,*) verbose
          end select
       end do
    end if
    close(10)

    call read_gas_composition_params(pfile)

    return
    
  end subroutine read_mesh_params


  
  subroutine print_mesh_params(unit)
    
    ! ---------------------------------------------------------------------------------
    ! write parameter values to std output or to an open file if argument unit is
    ! present.
    ! ---------------------------------------------------------------------------------

    integer(kind=4),optional,intent(in) :: unit

    if (present(unit)) then 
       write(unit,'(a,a,a)') '[mesh]'
       write(unit,'(a,L1)')  '  verbose    = ',verbose
       write(unit,'(a)')             ' '
       call print_gas_composition_params(unit)
    else
       write(*,'(a,a,a)') '[mesh]'
       write(*,'(a,L1)')  '  verbose    = ',verbose
       write(*,'(a)')             ' '
       call print_gas_composition_params
    end if
       
    return
    
  end subroutine print_mesh_params
  


end module module_mesh
