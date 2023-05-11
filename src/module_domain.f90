module module_domain

  use module_utils, only: path
  
  implicit none

  public
  
  type shell
     real(kind=8),dimension(3) :: center
     real(kind=8)              :: r_inbound,r_outbound
  end type shell

  type cube
     real(kind=8),dimension(3) :: center
     real(kind=8)              :: size            ! convention: size is the full size of the cube, corners are x+-size/2
  end type cube

  type sphere
     real(kind=8),dimension(3) :: center
     real(kind=8)              :: radius
  end type sphere

  type slab                                   ! infinite slab in xy direction (make use of periodic boundaries)
     real(kind=8)              :: zc          ! position of the slab in the z direction
     real(kind=8)              :: thickness   ! thickness of the slab in the z direction 
  end type slab

  type domain
     character(10) :: type  ! one of the possible shapes ('shell', 'cube', 'sphere', 'slab')
     type(shell)   :: sh
     type(cube)    :: cu
     type(sphere)  :: sp
     type(slab)    :: sl
  end type domain
  

contains
  
  !--------------------------------------------------------------------------------------------------
  ! domain constructors
  !--------------------------------------------------------------------------------------------------

  subroutine domain_constructor_from_scratch(dom,type,xc,yc,zc,r,r_inbound,r_outbound,size,thickness)

    implicit none
    character(10),intent(in)         :: type
    real(kind=8),intent(in),optional :: xc,yc,zc
    real(kind=8),intent(in),optional :: r                     ! parameters for sphere
    real(kind=8),intent(in),optional :: r_inbound,r_outbound  ! parameters for shell
    real(kind=8),intent(in),optional :: size                  ! parameters for cube
    real(kind=8),intent(in),optional :: thickness             ! parameters for slab
    type(domain),intent(out)         :: dom
    logical                          :: ok

    select case(trim(type))

    case('sphere')
       ! check if optional argument required for sphere are present
       ok = present(xc).and.present(yc).and.present(zc).and.present(r)
       if (.not.ok) then
          print *,'ERROR: arguments to construct a sphere domain are missing...'
          stop
       endif
       dom%type=type
       dom%sp%center(1)=xc
       dom%sp%center(2)=yc
       dom%sp%center(3)=zc
       dom%sp%radius=r

    case('shell')
       ! check if optional argument required for sphere are present
       ok = present(xc).and.present(yc).and.present(zc).and.present(r_inbound).and.present(r_outbound)
       if (.not.ok) then
          print *,'ERROR: arguments to construct a shell domain are missing...'
          stop
       endif
       dom%type=type
       dom%sh%center(1)=xc
       dom%sh%center(2)=yc
       dom%sh%center(3)=zc
       dom%sh%r_inbound=r_inbound
       dom%sh%r_outbound=r_outbound

    case('cube')
       ! check if optional argument required for sphere are present
       ok = present(xc).and.present(yc).and.present(zc).and.present(size)
       if (.not.ok) then
          print *,'ERROR: arguments to construct a cube domain are missing...'
          stop
       endif
       dom%type=type
       dom%cu%center(1)=xc
       dom%cu%center(2)=yc
       dom%cu%center(3)=zc
       dom%cu%size=size


    case('slab')
       ! check if optional argument required for sphere are present
       ok = present(zc).and.present(thickness)
       if (.not.ok) then
          print *,'ERROR: arguments to construct a slab domain are missing...'
          stop
       endif
       dom%type=type
       dom%sl%zc=zc
       dom%sl%thickness=thickness

    case default
       print *,'ERROR: type not defined ',trim(type)
       stop
    end select

    return

  end subroutine domain_constructor_from_scratch


  
  subroutine domain_constructor_from_file(filename,dom)
    ! read a domain file (filename) and initialise domain dom from it

    implicit none
    character(2000),intent(in) :: filename
    type(domain),intent(inout) :: dom

    open(unit=13, file=trim(filename), status='old', form='formatted', action='read')

    call read_domain(unit=13,dom=dom)

    close(13)

  end subroutine domain_constructor_from_file


  !--------------------------------------------------------------------------------------------------
  ! public use
  !--------------------------------------------------------------------------------------------------

  subroutine select_cells_in_domain(dom,n,xp,level,indsel)
    
    implicit none
    integer(kind=4),intent(in)                           :: n
    type(domain),intent(in)                              :: dom 
    real(kind=8),dimension(1:n,1:3),intent(in)           :: xp
    integer(kind=4),dimension(1:n),intent(in)            :: level
    integer(kind=4),dimension(:),allocatable,intent(out) :: indsel
    integer(kind=4)                                      :: i,ii,nsel
    integer(kind=4),dimension(:),allocatable             :: tmpi 
    real(kind=8),dimension(1:3)                          :: xc
    real(kind=8)                                         :: dx
    
    allocate(indsel(1:n))
    indsel=0
    ii=0
    do i=1,n
       xc(1:3) = xp(i,1:3)
       dx = 0.5d0**(level(i))
       if(level(i)>0)then
          if(domain_contains_cell(xc,dx,dom))then
             ii=ii+1
             indsel(ii)=i
          endif
       endif
    enddo
    nsel=ii
    allocate(tmpi(1:nsel))
    tmpi(1:nsel) = indsel(1:nsel)
    deallocate(indsel)
    allocate(indsel(1:nsel))
    indsel=tmpi
    deallocate(tmpi)
    
  end subroutine select_cells_in_domain


  subroutine read_domain(unit,dom)

    integer(kind=4),intent(in)  :: unit
    type(domain),intent(out)    :: dom 

    read(unit,*) dom%type
    select case(trim(dom%type))
    case('sphere')
       read(unit,*) dom%sp%center(:)
       read(unit,*) dom%sp%radius
    case('shell')
       read(unit,*) dom%sh%center(:)
       read(unit,*) dom%sh%r_inbound
       read(unit,*) dom%sh%r_outbound
    case('cube')
       read(unit,*) dom%cu%center(:)
       read(unit,*) dom%cu%size
    case('slab')
       read(unit,*) dom%sl%zc
       read(unit,*) dom%sl%thickness
    case default
       print *,'ERROR: type not defined ',trim(dom%type)
       stop
    end select

  end subroutine read_domain



  subroutine read_domain_bin(unit,dom)

    integer(kind=4),intent(in)  :: unit
    type(domain),intent(out)    :: dom 

    read(unit) dom%type
    select case(trim(dom%type))
    case('sphere')
       read(unit) dom%sp%center(:)
       read(unit) dom%sp%radius
    case('shell')
       read(unit) dom%sh%center(:)
       read(unit) dom%sh%r_inbound
       read(unit) dom%sh%r_outbound
    case('cube')
       read(unit) dom%cu%center(:)
       read(unit) dom%cu%size
    case('slab')
       read(unit) dom%sl%zc
       read(unit) dom%sl%thickness
    case default
       print *,'ERROR: type not defined ',trim(dom%type)
       stop
    end select

  end subroutine read_domain_bin



  subroutine dump_domain_bin(unit,dom)

    integer(kind=4),intent(in) :: unit
    type(domain),intent(in)    :: dom 

    write(unit) dom%type
    select case(trim(dom%type))
    case('sphere')
       write(unit) dom%sp%center(:)
       write(unit) dom%sp%radius
    case('shell')
       write(unit) dom%sh%center(:)
       write(unit) dom%sh%r_inbound
       write(unit) dom%sh%r_outbound
    case('cube')
       write(unit) dom%cu%center(:)
       write(unit) dom%cu%size
    case('slab')
       write(unit) dom%sl%zc
       write(unit) dom%sl%thickness
    case default 
       print *,'ERROR: type not defined ',trim(dom%type)
       stop
    end select

  end subroutine dump_domain_bin

  
  
  subroutine dump_domain_form(unit,dom)

    integer(kind=4),intent(in) :: unit
    type(domain),intent(in)    :: dom 

    write(unit,'(a)') trim(dom%type)
    select case(trim(dom%type))
    case('sphere')
       write(unit,'(3(f14.10))') dom%sp%center(:)
       write(unit,'(f14.10)') dom%sp%radius
    case('shell')
       write(unit,'(3(f14.10))') dom%sh%center(:)
       write(unit,'(f14.10)') dom%sh%r_inbound
       write(unit,'(f14.10)') dom%sh%r_outbound
    case('cube')
       write(unit,'(3(f14.10))') dom%cu%center(:)
       write(unit,'(f14.10)') dom%cu%size
    case('slab')
       write(unit,'(f14.10)') dom%sl%zc
       write(unit,'(f14.10)') dom%sl%thickness
    case default 
       print *,'ERROR: type not defined ',trim(dom%type)
       stop
    end select

  end subroutine dump_domain_form



  subroutine domain_write_file(file,dom)

    character(*),intent(in) :: file
    type(domain),intent(in) :: dom 

    open(unit=14, file=trim(file), status='unknown', form='formatted', action='write')
    call dump_domain_form(unit=14,dom=dom)
    close(14)

  end subroutine domain_write_file



  function domain_contains_point(x,dom)
    ! -> returns T/F if point xyz is in domain dom.
    type(domain),intent(in)              :: dom
    real(kind=8),dimension(3),intent(in) :: x
    logical                              :: domain_contains_point
    real(kind=8)                         :: rr,dx,dy,dz
    domain_contains_point=.false.
    select case(trim(dom%type))

    case('sphere')
       ! correct cell's position for periodic boundaries 
       dx = x(1)-dom%sp%center(1)
       if (dx > 0.5d0) then 
          dx = dx -1.0d0 
       else if (dx < -0.5d0) then 
          dx = dx + 1.0d0
       end if
       dy = x(2)-dom%sp%center(2)
       if (dy > 0.5d0) then 
          dy = dy -1.0d0 
       else if (dy < -0.5d0) then 
          dy = dy + 1.0d0
       end if
       dz = x(3)-dom%sp%center(3)
       if (dz > 0.5d0) then 
          dz = dz -1.0d0 
       else if (dz < -0.5d0) then 
          dz = dz + 1.0d0
       end if
       !rr = (x(1)-dom%sp%center(1))**2 + (x(2)-dom%sp%center(2))**2 + (x(3)-dom%sp%center(3))**2
       rr = dx**2 + dy**2 + dz**2
       if(rr<dom%sp%radius*dom%sp%radius)domain_contains_point=.true.

    case('shell')
       ! correct positions for periodic boundaries 
       dx = x(1)-dom%sh%center(1)
       if (dx > 0.5d0) then 
          dx = dx -1.0d0 
       else if (dx < -0.5d0) then 
          dx = dx + 1.0d0
       end if
       dy = x(2)-dom%sh%center(2)
       if (dy > 0.5d0) then 
          dy = dy -1.0d0 
       else if (dy < -0.5d0) then 
          dy = dy + 1.0d0
       end if
       dz = x(3)-dom%sh%center(3)
       if (dz > 0.5d0) then 
          dz = dz -1.0d0 
       else if (dz < -0.5d0) then 
          dz = dz + 1.0d0
       end if
       rr = dx*dx + dy*dy + dz*dz
       !!rr = (x(1)-dom%sh%center(1))**2 + (x(2)-dom%sh%center(2))**2 + (x(3)-dom%sh%center(3))**2
       if((rr>dom%sh%r_inbound*dom%sh%r_inbound).and.(rr<dom%sh%r_outbound*dom%sh%r_outbound))domain_contains_point=.true.

    case('cube')
       ! correct positions for periodic boundaries 
       dx = x(1)-dom%cu%center(1)
       if (dx > 0.5d0) then 
          dx = dx -1.0d0 
       else if (dx < -0.5d0) then 
          dx = dx + 1.0d0
       end if
       if (abs(dx) < dom%cu%size*0.5d0) then 
          dx = x(2)-dom%cu%center(2)
          if (dx > 0.5d0) then 
             dx = dx -1.0d0 
          else if (dx < -0.5d0) then 
             dx = dx + 1.0d0
          end if
          if (abs(dx) < dom%cu%size*0.5d0) then 
             dx = x(3)-dom%cu%center(3)
             if (dx > 0.5d0) then 
                dx = dx -1.0d0 
             else if (dx < -0.5d0) then 
                dx = dx + 1.0d0
             end if
             if (abs(dx) < dom%cu%size*0.5d0) then
                domain_contains_point=.true.
             end if
          end if
       end if
       
    case('slab')
       dz = x(3) - dom%sl%zc
       if (dz > 0.5d0) then 
          dz = dz - 1.0d0
       else if (dz < -0.5d0) then 
          dz = dz + 1.0d0
       end if
       if(abs(dz) < dom%sl%thickness*0.5d0)domain_contains_point=.true.
    end select
    return
  end function domain_contains_point

  

  function domain_fully_contains_cell(x,dx,dom)
    ! -> returns T/F if the full cell at x, of size dx, is in domain dom.
    type(domain),intent(in)              :: dom
    real(kind=8),dimension(3),intent(in) :: x
    real(kind=8),intent(in)              :: dx
    logical                              :: domain_fully_contains_cell
    real(kind=8)                         :: rr,xc,dd,ddx,ddy,ddz
    real(kind=8),parameter :: sqrt3over2 = sqrt(3.0d0)*0.5d0
    
    domain_fully_contains_cell=.false.

    select case(trim(dom%type))

    case('sphere')
       ! correct cell's position for periodic boundaries
       ddx = x(1)-dom%sp%center(1)
       if (ddx > 0.5d0) then 
          ddx = ddx -1.0d0 
       else if (ddx < -0.5d0) then 
          ddx = ddx + 1.0d0
       end if
       ddy = x(2)-dom%sp%center(2)
       if (ddy > 0.5d0) then 
          ddy = ddy -1.0d0 
       else if (ddy < -0.5d0) then 
          ddy = ddy + 1.0d0
       end if
       ddz = x(3)-dom%sp%center(3)
       if (ddz > 0.5d0) then 
          ddz = ddz -1.0d0 
       else if (ddz < -0.5d0) then 
          ddz = ddz + 1.0d0
       end if
       rr = sqrt(ddx**2 + ddy**2 + ddz**2)
       rr = rr + dx*sqrt3over2
       if (rr <= dom%sp%radius) domain_fully_contains_cell=.true.

    case('shell')
       ! correct cell's position for periodic boundaries
       ddx = x(1)-dom%sh%center(1)
       if (ddx > 0.5d0) then 
          ddx = ddx -1.0d0 
       else if (ddx < -0.5d0) then 
          ddx = ddx + 1.0d0
       end if
       ddy = x(2)-dom%sh%center(2)
       if (ddy > 0.5d0) then 
          ddy = ddy -1.0d0 
       else if (ddy < -0.5d0) then 
          ddy = ddy + 1.0d0
       end if
       ddz = x(3)-dom%sh%center(3)
       if (ddz > 0.5d0) then 
          ddz = ddz -1.0d0 
       else if (ddz < -0.5d0) then 
          ddz = ddz + 1.0d0
       end if
       rr = sqrt(ddx*ddx + ddy*ddy + ddz*ddz)
       !!!rr = sqrt((x(1)-dom%sh%center(1))**2 + (x(2)-dom%sh%center(2))**2 + (x(3)-dom%sh%center(3))**2)
       if(( (rr-dx*sqrt3over2)>=dom%sh%r_inbound) .and. ((rr+dx*sqrt3over2)<=dom%sh%r_outbound) ) then
          domain_fully_contains_cell=.true.
       end if

    case('cube')
       ! correct cell's position for periodic boundaries 
       dd = x(1) - dom%cu%center(1)
       xc = x(1)
       if (dd > 0.5d0) then 
          xc = xc -1.0d0 
       else if (dd < -0.5d0) then 
          xc = xc + 1.0d0
       end if
       if ((xc+dx*0.5d0 <= dom%cu%center(1)+dom%cu%size*0.5d0).and. &
            (xc-dx*0.5d0 >= dom%cu%center(1)-dom%cu%size*0.5d0)) then
          dd = x(2) - dom%cu%center(2)
          xc = x(2)
          if (dd > 0.5d0) then 
             xc = xc -1.0d0 
          else if (dd < -0.5d0) then 
             xc = xc + 1.0d0
          end if
          if ((xc+dx*0.5d0 <= dom%cu%center(2)+dom%cu%size*0.5d0).and. &
               (xc-dx*0.5d0 >= dom%cu%center(2)-dom%cu%size*0.5d0)) then
             dd = x(3) - dom%cu%center(3)
             xc = x(3)
             if (dd > 0.5d0) then 
                xc = xc -1.0d0 
             else if (dd < -0.5d0) then 
                xc = xc + 1.0d0
             end if
             if ((xc+dx*0.5d0 <= dom%cu%center(3)+dom%cu%size*0.5d0).and. &
                  (xc-dx*0.5d0 >= dom%cu%center(3)-dom%cu%size*0.5d0)) then
                domain_fully_contains_cell=.true.
             end if
          end if
       end if
       
    case('slab')
       dd = x(3) - dom%sl%zc
       xc = x(3)
       if (dd > 0.5d0) then 
          xc = xc -1.0d0 
       else if (dd < -0.5d0) then 
          xc = xc + 1.0d0
       end if
       if((xc+dx*0.5d0 <= dom%sl%zc+dom%sl%thickness*0.5d0).and. &
            (xc-dx*0.5d0 >= dom%sl%zc-dom%sl%thickness*0.5d0)) domain_fully_contains_cell=.true.
       
    end select

    return
  end function domain_fully_contains_cell



  function domain_contains_cell(x,dx,dom)
    ! -> returns T/F if cell at x, of size dx, is (partially or fully) in domain dom.
    ! warning: it is not a strict condition, some cell may be completely outside the domain. 
    ! purpose of this function: selecting all cells belonging to a domain, if some cells are
    !    counted in and should not, we don't care.
    type(domain),intent(in)              :: dom
    real(kind=8),dimension(3),intent(in) :: x
    real(kind=8),intent(in)              :: dx
    logical                              :: domain_contains_cell
    real(kind=8)                         :: rr,dd,ddx,ddy,ddz,rcell
    real(kind=8),parameter :: sqrt3over2 = sqrt(3.0d0)*0.5d0
    
    domain_contains_cell=.false.
    
    select case(trim(dom%type))
       
    case('sphere')
       ! correct cell's position for periodic boundaries
       ddx = x(1)-dom%sp%center(1)
       if (ddx > 0.5d0) then 
          ddx = ddx -1.0d0 
       else if (ddx < -0.5d0) then 
          ddx = ddx + 1.0d0
       end if
       ddy = x(2)-dom%sp%center(2)
       if (ddy > 0.5d0) then 
          ddy = ddy -1.0d0 
       else if (ddy < -0.5d0) then 
          ddy = ddy + 1.0d0
       end if
       ddz = x(3)-dom%sp%center(3)
       if (ddz > 0.5d0) then 
          ddz = ddz -1.0d0 
       else if (ddz < -0.5d0) then 
          ddz = ddz + 1.0d0
       end if
       rr = sqrt(ddx**2 + ddy**2 + ddz**2)
       rcell = dx*sqrt3over2
       if (rr <= (dom%sp%radius + rcell)) domain_contains_cell=.true.
       
    case('shell')
       ! correct cell's position for periodic boundaries
       ddx = x(1)-dom%sh%center(1)
       if (ddx > 0.5d0) then 
          ddx = ddx -1.0d0 
       else if (ddx < -0.5d0) then 
          ddx = ddx + 1.0d0
       end if
       ddy = x(2)-dom%sh%center(2)
       if (ddy > 0.5d0) then 
          ddy = ddy -1.0d0 
       else if (ddy < -0.5d0) then 
          ddy = ddy + 1.0d0
       end if
       ddz = x(3)-dom%sh%center(3)
       if (ddz > 0.5d0) then 
          ddz = ddz -1.0d0 
       else if (ddz < -0.5d0) then 
          ddz = ddz + 1.0d0
       end if
       rr = sqrt(ddx*ddx + ddy*ddy + ddz*ddz)
       rcell = dx*sqrt3over2
       if((rr>=(dom%sh%r_inbound-rcell)) .and. (rr<=(dom%sh%r_outbound+rcell))) then
          domain_contains_cell=.true.
       end if
       
    case('cube')
       ! correct cell's position for periodic boundaries 
       dd = x(1) - dom%cu%center(1)
       if (dd > 0.5d0) then 
          dd = dd - 1.0d0 
       else if (dd < -0.5d0) then 
          dd = dd + 1.0d0
       end if
       if (abs(dd)<=(dx*0.5d0 + dom%cu%size*0.5d0)) then
          dd = x(2) - dom%cu%center(2)
          if (dd > 0.5d0) then 
             dd = dd - 1.0d0 
          else if (dd < -0.5d0) then 
             dd = dd + 1.0d0
          end if
          if (abs(dd)<=(dx*0.5d0 + dom%cu%size*0.5d0)) then
             dd = x(3) - dom%cu%center(3)
             if (dd > 0.5d0) then 
                dd = dd - 1.0d0 
             else if (dd < -0.5d0) then 
                dd = dd + 1.0d0
             end if
             if (abs(dd)<=(dx*0.5d0 + dom%cu%size*0.5d0)) then
                domain_contains_cell=.true.
             end if
          end if
       end if
       
    case('slab')
       dd = x(3) - dom%sl%zc
       if (dd > 0.5d0) then 
          dd = dd - 1.0d0 
       else if (dd < -0.5d0) then 
          dd = dd + 1.0d0
       end if
       if(abs(dd) <= (dx*0.5d0 + dom%sl%thickness*0.5d0)) then
          domain_contains_cell=.true.
       endif
       
    end select
    return
  end function domain_contains_cell



  function get_my_new_domain(x,liste_domaines)
    ! Return the domain from list of domains “liste_domaines” in which point x is most deeply embedded. 
    ! (i.e. such that the distance of x to domain border is the largest). 
    type(domain),intent(in),dimension(:) :: liste_domaines
    real(kind=8),dimension(3),intent(in) :: x
    integer(kind=4)                      :: get_my_new_domain, ndom, i, imax
    real(kind=8)                         :: dmax,db

    ndom = size(liste_domaines)
    dmax = 0.0d0
    imax = 0
    do i=1,ndom
       db = domain_distance_to_border(x,liste_domaines(i))
       ! if db < 0 point outside domain
       if (db>=dmax)then
          imax = i
          dmax = db
       endif
    enddo
    if(dmax<0.0d0 .or. imax==0)then
       print *,'ERROR: problem with get_my_new_domain'
       stop
    endif
    get_my_new_domain = imax

    return
  end function get_my_new_domain


  
  function domain_distance_to_border(x,dom)
    ! return distance of point x to the closest border of domain dom
    ! convention: negative distance means outside domain
    real(kind=8),dimension(3),intent(in) :: x
    type(domain),intent(in)              :: dom
    real(kind=8)                         :: domain_distance_to_border, rr, ddx, ddy, ddz,xc,yc,zc

    select case(trim(dom%type))

    case('sphere')
       ! correct position for periodic boundaries
       ddx = x(1)-dom%sp%center(1)
       if (ddx > 0.5d0) then 
          ddx = ddx -1.0d0 
       else if (ddx < -0.5d0) then 
          ddx = ddx + 1.0d0
       end if
       ddy = x(2)-dom%sp%center(2)
       if (ddy > 0.5d0) then 
          ddy = ddy -1.0d0 
       else if (ddy < -0.5d0) then 
          ddy = ddy + 1.0d0
       end if
       ddz = x(3)-dom%sp%center(3)
       if (ddz > 0.5d0) then 
          ddz = ddz -1.0d0 
       else if (ddz < -0.5d0) then 
          ddz = ddz + 1.0d0
       end if
       !rr = sqrt((x(1)-dom%sp%center(1))**2 + (x(2)-dom%sp%center(2))**2 + (x(3)-dom%sp%center(3))**2)
       rr = sqrt(ddx**2 + ddy**2 + ddz**2)
       domain_distance_to_border = dom%sp%radius - rr
       
    case('shell')
       ! correct position for periodic boundaries
       ddx = x(1)-dom%sh%center(1)
       if (ddx > 0.5d0) then 
          ddx = ddx -1.0d0 
       else if (ddx < -0.5d0) then 
          ddx = ddx + 1.0d0
       end if
       ddy = x(2)-dom%sh%center(2)
       if (ddy > 0.5d0) then 
          ddy = ddy -1.0d0 
       else if (ddy < -0.5d0) then 
          ddy = ddy + 1.0d0
       end if
       ddz = x(3)-dom%sh%center(3)
       if (ddz > 0.5d0) then 
          ddz = ddz -1.0d0 
       else if (ddz < -0.5d0) then 
          ddz = ddz + 1.0d0
       end if
       rr = sqrt(ddx*ddx + ddy*ddy + ddz*ddz)
       !!!rr = sqrt((x(1)-dom%sh%center(1))**2 + (x(2)-dom%sh%center(2))**2 + (x(3)-dom%sh%center(3))**2)
       domain_distance_to_border = min((rr-dom%sh%r_inbound),(dom%sh%r_outbound-rr))
       
    case('cube')
       ! correct cell's position for periodic boundaries 
       xc = x(1)
       ddx = xc - dom%cu%center(1)
       if (ddx > 0.5d0) then 
          xc = xc -1.0d0 
       else if (ddx < -0.5d0) then 
          xc = xc + 1.0d0
       end if
       yc = x(2)
       ddy = yc - dom%cu%center(2)
       if (ddy > 0.5d0) then 
          yc = yc -1.0d0 
       else if (ddy < -0.5d0) then 
          yc = yc + 1.0d0
       end if
       zc = x(3)
       ddz = zc - dom%cu%center(3)
       if (ddz > 0.5d0) then 
          zc = zc -1.0d0 
       else if (ddz < -0.5d0) then 
          zc = zc + 1.0d0
       end if
       ddx = min((dom%cu%center(1)+dom%cu%size*0.5d0-xc), (xc-dom%cu%center(1)+dom%cu%size*0.5d0))
       ddy = min((dom%cu%center(2)+dom%cu%size*0.5d0-yc), (yc-dom%cu%center(2)+dom%cu%size*0.5d0))
       ddz = min((dom%cu%center(3)+dom%cu%size*0.5d0-zc), (zc-dom%cu%center(3)+dom%cu%size*0.5d0))
       if((ddx>=0.0d0).and.(ddy>=0.0d0).and.(ddz>=0.0d0))then
          ! inside domain
          domain_distance_to_border = min(ddx,ddy,ddz)
       else
          ! outside domain
          domain_distance_to_border = -sqrt((min(0.0d0,ddx))**2 + (min(0.0d0,ddy))**2 + (min(0.0d0,ddz))**2)  
       endif

    case('slab')
       zc = x(3)
       ddz = zc - dom%sl%zc
       if (ddz > 0.5d0) then 
          zc = zc -1.0d0 
       else if (ddz < -0.5d0) then 
          zc = zc + 1.0d0
       end if
       domain_distance_to_border = min((dom%sl%zc+dom%sl%thickness*0.5d0-zc), (zc-dom%sl%zc+dom%sl%thickness*0.5d0)) 

    end select

    return
  end function domain_distance_to_border



  function domain_distance_to_border_along_k(x,k,dom)
    ! return the distance of point x to the closest border of domain dom along propagation vector k.
    ! inputs:
    ! - x : position vector, in simulation-box units
    ! - k : direction vector, normalized
    ! - dom : a domain
    
    real(kind=8),dimension(3),intent(in) :: x, k
    type(domain),intent(in)              :: dom
    real(kind=8)                         :: domain_distance_to_border_along_k
    ! variables for the spherical case
    real(kind=8) :: b, c, delta, dx, dy, dz, ddx
    ! variables for the shell case
    real(kind=8) :: t1,t2,tin,tout
    ! variables for the cube case
    real(kind=8),dimension(3) :: x_dom,xc
    integer(kind=4) :: i
    ! variables for the slab case
    real(kind=8) :: zc
    
    ! point x should be inside domain dom
    if(.not.(domain_contains_point(x,dom)))then
       print*,'ERROR: point x outside domain dom in domain_distance_to_border_along_k...'
       stop
    endif
    
    select case(trim(dom%type))
       
    case('sphere')
       ! correct position for periodic boundaries
       dx = x(1)-dom%sp%center(1)
       if (dx > 0.5d0) then 
          dx = dx -1.0d0 
       else if (dx < -0.5d0) then 
          dx = dx + 1.0d0
       end if
       dy = x(2)-dom%sp%center(2)
       if (dy > 0.5d0) then 
          dy = dy -1.0d0 
       else if (dy < -0.5d0) then 
          dy = dy + 1.0d0
       end if
       dz = x(3)-dom%sp%center(3)
       if (dz > 0.5d0) then 
          dz = dz -1.0d0 
       else if (dz < -0.5d0) then 
          dz = dz + 1.0d0
       end if
       !dx = x(1) - dom%sp%center(1)
       !dy = x(2) - dom%sp%center(2)
       !dz = x(3) - dom%sp%center(3)
       b = 2.d0 * ( k(1)*dx + k(2)*dy +  k(3)*dz )
       c = dx*dx + dy*dy + dz*dz - dom%sp%radius*dom%sp%radius
       delta = b*b - 4.0d0*c
       if (delta <= 0.0d0) then
          print*,'ERROR: pb with domain_distance_to_border_along_k' ! This can only mean that x is out of the domain ... 
          stop
       end if
       domain_distance_to_border_along_k = (-b + sqrt(delta))*0.5d0 ! select the only positive solution

    case('shell')
       ! correct position for periodic boundaries
       dx = x(1) - dom%sh%center(1)
       if (dx > 0.5d0) then 
          dx = dx -1.0d0 
       else if (dx < -0.5d0) then 
          dx = dx + 1.0d0
       end if
       dy = x(2) - dom%sh%center(2)
       if (dy > 0.5d0) then 
          dy = dy -1.0d0 
       else if (dy < -0.5d0) then 
          dy = dy + 1.0d0
       end if
       dz = x(3) - dom%sh%center(3)
       if (dz > 0.5d0) then 
          dz = dz -1.0d0 
       else if (dz < -0.5d0) then 
          dz = dz + 1.0d0
       end if
       !dx = x(1) - dom%sh%center(1)
       !dy = x(2) - dom%sh%center(2)
       !dz = x(3) - dom%sh%center(3)
       ! inner shell intersection
       b = 2.0d0 * ( k(1)*dx + k(2)*dy +  k(3)*dz )
       c = dx*dx + dy*dy + dz*dz - dom%sh%r_inbound*dom%sh%r_inbound
       delta = b*b - 4.0d0*c
       tin   = 2.0d0 ! larger than box size 
       if (delta >= 0.0d0) then
          t1 =  (-b + sqrt(delta))*0.5d0
          t2 = (-b - sqrt(delta))*0.5d0
          ! select smallest positive solution
          if (t1 > 0 .and. t1 < tin) tin = t1
          if (t2 > 0 .and. t2 < tin) tin = t2
       end if
       ! outer shell intersection
       b = 2.0d0 * ( k(1)*dx + k(2)*dy +  k(3)*dz )
       c = dx*dx + dy*dy + dz*dz - dom%sh%r_outbound*dom%sh%r_outbound
       delta = b*b - 4.0d0*c
       tout  = 2.0d0 ! larger than box size 
       if (delta >= 0.0d0) then
          t1 =  (-b + sqrt(delta))*0.5d0
          t2 = (-b - sqrt(delta))*0.5d0
          ! select smallest positive solution
          if (t1 > 0 .and. t1 < tout) tout = t1
          if (t2 > 0 .and. t2 < tout) tout = t2
       end if
       ! Select smallest distance
       domain_distance_to_border_along_k = min(tin,tout)
       
    case('slab')
       zc = x(3)
       dz = zc - dom%sl%zc
       if (dz > 0.5d0) then 
          zc = zc - 1.0d0 
       else if (dz < -0.5d0) then 
          zc = zc + 1.0d0
       end if
       if(k(3)<0.0d0)then
          domain_distance_to_border_along_k = (dom%sl%zc-dom%sl%thickness*0.5d0-zc)/k(3)
       else
          domain_distance_to_border_along_k = (dom%sl%zc+dom%sl%thickness*0.5d0-zc)/k(3)
       endif
       if (domain_distance_to_border_along_k < 0.0d0)then
          print *,'ERROR: pb with distance to border along k, slab case, point outside domain...'
          stop
       endif

    case('cube')  
       ! correct position for periodic boundaries
       do i = 1,3
          xc(i) = x(i)
          ddx = xc(i) - dom%cu%center(i)
          if (ddx > 0.5d0) then 
             xc(i) = xc(i) -1.0d0 
          else if (ddx < -0.5d0) then 
             xc(i) = xc(i) + 1.0d0
          end if
       end do
       ! get position relative to the domain
       x_dom = (xc - dom%cu%center)/dom%cu%size + 0.5d0
       if((x_dom(1) < 0.0d0).or.(x_dom(1)>1.0d0).or.&
            (x_dom(2) < 0.0d0).or.(x_dom(2)>1.0d0).or.&
            (x_dom(3) < 0.0d0).or.(x_dom(3)>1.0d0))then
          ! x is outside the domain
          print *,'ERROR: pb with distance to border along k, cube case, point outside domain...'
          stop
       endif
       domain_distance_to_border_along_k = path(x_dom,k) * dom%cu%size
       
    case default
       print *,'ERROR: domain type not defined ',trim(dom%type)
       stop

    end select
       
    return
    
  end function domain_distance_to_border_along_k



  subroutine domain_get_bounding_box(dom,xmin,xmax,ymin,ymax,zmin,zmax)
    implicit none
    type(domain),intent(in)     :: dom
    real(kind=8),intent(inout)  :: xmin,xmax,ymin,ymax,zmin,zmax

    select case(dom%type)
    case('sphere')
       xmax = dom%sp%center(1) + dom%sp%radius
       xmin = dom%sp%center(1) - dom%sp%radius
       ymax = dom%sp%center(2) + dom%sp%radius
       ymin = dom%sp%center(2) - dom%sp%radius
       zmax = dom%sp%center(3) + dom%sp%radius
       zmin = dom%sp%center(3) - dom%sp%radius
    case('shell')
       xmax = dom%sh%center(1) + dom%sh%r_outbound
       xmin = dom%sh%center(1) - dom%sh%r_outbound
       ymax = dom%sh%center(2) + dom%sh%r_outbound
       ymin = dom%sh%center(2) - dom%sh%r_outbound
       zmax = dom%sh%center(3) + dom%sh%r_outbound
       zmin = dom%sh%center(3) - dom%sh%r_outbound
    case('cube')
       xmax = dom%cu%center(1) + dom%cu%size*0.5d0
       xmin = dom%cu%center(1) - dom%cu%size*0.5d0
       ymax = dom%cu%center(2) + dom%cu%size*0.5d0
       ymin = dom%cu%center(2) - dom%cu%size*0.5d0
       zmax = dom%cu%center(3) + dom%cu%size*0.5d0
       zmin = dom%cu%center(3) - dom%cu%size*0.5d0
    case('slab')
       xmax = 1.0d0
       xmin = 0.0d0
       ymax = 1.0d0
       ymin = 0.0d0
       zmax = dom%sl%zc + dom%sl%thickness*0.5d0
       zmin = dom%sl%zc - dom%sl%thickness*0.5d0
    end select
    
    return
  end subroutine domain_get_bounding_box



end module module_domain

