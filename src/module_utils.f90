module module_utils
  
  ! general-purpose functions : 
  ! 
  ! - isotropic_direction
  ! - anisotropic_direction_HIcore
  ! - anisotropic_direction_Rayleigh
  ! - anisotropic_direction_Dust
  ! - locatedb
  ! - path
  ! - print_rascas_header
  
  use module_constants, only : twopi
  use module_random
  
  public

contains

  
  subroutine isotropic_direction(k,iran)

    ! ---------------------------------------------------------------------------------
    ! return k vector pointing to a random direction (uniform on the sphere)
    ! ---------------------------------------------------------------------------------
    ! INPUTS:
    ! - iran : state of random number generator
    ! OUTPUTS :
    ! - k    : normalized direction vector
    ! - iran : updated state of random number generator
    ! ---------------------------------------------------------------------------------

    implicit none

    real(kind=8),intent(out)      :: k(3)
    integer(kind=4),intent(inout) :: iran
    real(kind=8)                  :: cos_theta,sin_theta,phi,knorm
    
    phi   = twopi*ran3(iran)
    cos_theta = 1.0d0 - 2.0d0 * ran3(iran)  ! in [-1,1]
    sin_theta = sqrt(1.0d0 - cos_theta**2) ! in [0,1]
    k(1) = sin_theta * cos(phi)   !x
    k(2) = sin_theta * sin(phi)   !y
    k(3) = cos_theta              !z
    ! force normalisation at numerical precision
    knorm = sqrt(k(1)*k(1)+k(2)*k(2)+k(3)*k(3))
    k     = k / knorm

  end subroutine isotropic_direction

  

  subroutine anisotropic_direction_HIcore(kin,kout,mu,bu,iran)

    ! ---------------------------------------------------------------------------------
    ! sends back new direction vector kout as a function of incident direction kin, for a phase function
    ! described by P(mu) = 11/24 + 3/24 * mu**2 (with mu = cos(theta) = kin.kout). This is a good
    ! description, e.g., of core scatterings on HI atoms.
    ! ---------------------------------------------------------------------------------
    ! INPUTS:
    ! - kin  : normalized direction vector of incident photon
    ! - iran : state of random number generator
    ! OUTPUTS:
    ! - kout : normalized direction vector of scattered photon
    ! - mu   : dot-product between kin and kout (i.e. cos(theta))
    ! - bu   : sin(theta) (i.e. sqrt(1-mu**2))
    ! - iran : updated state of random number generator
    ! 
    ! Notes on the method : 
    ! ---------------------
    ! To draw theta values, we compute the cumulative probability from above :
    ! P(< mu) = 1/2 + 11/24 * mu + 1/24 * mu**3
    ! It is analytically integrable and invertible. The solution of this cubic polynomial
    ! is a function of the form: mu = (A+B)**(1.d0/3.d0) - (A-B)**(1.d0/3.d0)
    ! with B = (2*x-1)*6 and A = sqrt(B**2+11**3/27)
    ! -> to get a value of mu (hence theta), we draw x in [0,1] and compute mu from the above function
    ! ---------------------------------------------------------------------------------

    implicit none

    real(kind=8),intent(in)       :: kin(3)
    real(kind=8),intent(out)      :: kout(3)
    real(kind=8),intent(out)      :: mu,bu
    integer(kind=4),intent(inout) :: iran
    real(kind=8)                  :: phi,x,cti,sti,cpi,spi,ct1,st1,cp1,sp1,knorm,A,B
    
    phi = twopi * ran3(iran)
    x   = ran3(iran)
    B   = (2.d0*x - 1.d0)*6.d0
    A   = sqrt(B*B+11.d0**3/27.d0)
    mu  = (A+B)**(1.d0/3.d0) - (A-B)**(1.d0/3.d0)
    ! angular description of kin in external frame (box coordinates)
    cti = kin(3)
    sti = sqrt(1.0d0 - cti**2)  ! sin(theta) is positive for theta in [0,pi]. 
    if (sti > 0) then 
       cpi = kin(1)/sti
       spi = kin(2)/sti
    else
       cpi = 1.0d0
       spi = 0.0d0
    end if
    ! angular description of kout (relative to k)
    ct1 = mu
    st1 = sqrt(1.0d0 - ct1*ct1)
    bu  = st1
    cp1 = cos(phi)
    sp1 = sin(phi)
    ! vector kout (such that indeed knew . k = ct1) in external frame (box coords.)
    kout(1) = cti*cpi*st1*cp1 + sti*cpi*ct1 - spi*st1*sp1
    kout(2) = cti*spi*st1*cp1 + sti*spi*ct1 + cpi*st1*sp1
    kout(3) = -sti*st1*cp1 + cti*ct1
    ! force normalisation at numerical precision
    knorm = sqrt(kout(1)*kout(1)+kout(2)*kout(2)+kout(3)*kout(3))
    kout  = kout / knorm
    
  end subroutine anisotropic_direction_HIcore

  
  subroutine anisotropic_direction_Rayleigh(kin,kout,mu,bu,iran)

    ! ---------------------------------------------------------------------------------
    ! Sends back new direction vector kout as a function of incident direction kin, for a phase function
    ! described by P(mu) = 3/8 * (1 + mu**2) (with mu = cos(theta) = kin.kout). This is a good
    ! description, e.g., of wing scatterings on HI atoms.
    ! It is actually the phase function of Rayleigh scattering.
    ! ---------------------------------------------------------------------------------
    ! INPUTS:
    ! - kin  : normalized direction vector of incident photon
    ! - iran : state of random number generator
    ! OUTPUTS:
    ! - kout : normalized direction vector of scattered photon
    ! - mu   : dot-product between kin and kout (i.e. cos(theta))
    ! - bu   : sin(theta) (i.e. sqrt(1-mu**2))
    ! - iran : updated state of random number generator
    ! 
    ! Notes on the method : 
    ! ---------------------
    ! To draw theta values, we compute the cumulative probability from above :
    ! P(< mu) = 1/2 + 3/8 * mu + 1/8 * mu**3
    ! It is analytically integrable and invertible. The solution of this cubic polynomial
    ! is a function of the form: mu = (A+B)**(1.d0/3.d0) - (A-B)**(1.d0/3.d0)
    ! with B = 4*x-2 and A = sqrt(B**2+1)
    ! -> to get a value of mu (hence theta), we draw x in [0,1] and compute mu from the above fit
    ! ---------------------------------------------------------------------------------

    implicit none
    
    real(kind=8),intent(in)       :: kin(3)
    real(kind=8),intent(out)      :: kout(3)
    real(kind=8),intent(out)      :: mu,bu
    integer(kind=4),intent(inout) :: iran
    real(kind=8)                  :: phi,x,cti,sti,cpi,spi,ct1,st1,cp1,sp1,knorm,A,B

    phi = twopi * ran3(iran)
    x   = ran3(iran)
    B   = 4.d0*x - 2.d0
    A   = sqrt(B*B+1.d0)
    mu  = (A+B)**(1.d0/3.d0) - (A-B)**(1.d0/3.d0)
    ! angular description of kin in external frame (box coordinates)
    cti = kin(3)
    sti = sqrt(1.0d0 - cti**2)  ! sin(theta) is positive for theta in [0,pi]. 
    if (sti > 0) then 
       cpi = kin(1)/sti
       spi = kin(2)/sti
    else
       cpi = 1.0d0
       spi = 0.0d0
    end if
    ! angular description of kout (relative to k)
    ct1 = mu
    st1 = sqrt(1.0d0 - ct1*ct1)
    bu  = st1
    cp1 = cos(phi)
    sp1 = sin(phi)
    ! vector kout (such that indeed knew . k = ct1) in external frame (box coords.)
    kout(1) = cti*cpi*st1*cp1 + sti*cpi*ct1 - spi*st1*sp1
    kout(2) = cti*spi*st1*cp1 + sti*spi*ct1 + cpi*st1*sp1
    kout(3) = -sti*st1*cp1 + cti*ct1
    ! force normalisation at numerical precision
    knorm = sqrt(kout(1)*kout(1)+kout(2)*kout(2)+kout(3)*kout(3))
    kout  = kout / knorm

  end subroutine anisotropic_direction_Rayleigh


  subroutine anisotropic_direction_Dust(kin,kout,mu,iran,g_dust)

    ! -------------------------------------------------------------------------------------------------
    ! Returns new direction vector kout as a function of incident direction kin, for a phase function
    ! given by Henyey-Greenstein.
    ! INPUTS:
    ! - kin    : normalized direction vector of incident photon
    ! - iran   : state of random number generator
    ! - g_dust : g parameter of the Henyey-Greenstein phase function for dust scattering
    ! OUTPUTS:
    ! - kout : normalized direction vector of scattered photon
    ! - mu   : dot-product between kin and kout (i.e. cos(theta))
    ! - iran : updated state of random number generator
    ! -------------------------------------------------------------------------------------------------

    implicit none

    real(kind=8),intent(in)       :: kin(3)
    real(kind=8),intent(out)      :: kout(3)
    real(kind=8),intent(out)      :: mu
    integer(kind=4),intent(inout) :: iran
    real(kind=8)                  :: phi,cti,sti,cpi,spi,ct1,st1,cp1,sp1,x,knorm
    real(kind=8),intent(in)       :: g_dust


    ! determine scattering angle (in atom's frame)
    ! use White 79 approximation for the "reciprocal" of cumulative Henyey-Greenstein phase fct:
    x  = ran3(iran) 
    mu = (1.0d0+g_dust*g_dust-((1.0d0-g_dust*g_dust)/(1.0d0-g_dust+2.0d0*g_dust*x))**2)/(2.0d0*g_dust)

    ! angular description of kin in external frame (box coordinates)
    ! ---------------------------------------------------------------------------------
    ! kx = sin(theta) * cos(phi)
    ! ky = sin(theta) * sin(phi)
    ! kz = cos(theta)
    ! cti, sti, cpi, spi correspond to kin
    ! ---------------------------------------------------------------------------------
    cti = kin(3)              
    sti = sqrt(1.0d0 - cti**2)  ! sin(theta) is positive for theta in [0,pi]. 
    if (sti > 0) then 
       cpi = kin(1)/sti       
       spi = kin(2)/sti       
    else
       cpi = 1.0d0
       spi = 0.0d0
    end if
    
    ! angular description of kout (relative to kin)
    ct1 = mu
    st1 = sqrt(1.0d0 - ct1*ct1)
    phi = twopi * ran3(iran)
    cp1 = cos(phi)
    sp1 = sin(phi)
    
    ! vector kout (such that indeed kout . kin = ct1) in external frame (box coords.)
    kout(1) = cti*cpi*st1*cp1 + sti*cpi*ct1 - spi*st1*sp1
    kout(2) = cti*spi*st1*cp1 + sti*spi*ct1 + cpi*st1*sp1
    kout(3) = -sti*st1*cp1 + cti*ct1
    ! force normalisation at numerical precision
    knorm = sqrt(kout(1)*kout(1)+kout(2)*kout(2)+kout(3)*kout(3))
    kout  = kout / knorm

  end subroutine anisotropic_direction_Dust


  subroutine locatedb(xx,n,x,j)
    
    ! subroutine which locates the position j of a value x in an array xx of n elements
    ! NB : here xx is double precision
    
    implicit none
    
    integer(kind=4),intent(in)           :: n
    integer(kind=4),intent(out)          :: j
    integer(kind=4)                      :: jl,ju,jm
    real(kind=8),intent(in)              :: x
    real(kind=8),dimension(n),intent(in) :: xx
    
    
    jl = 0
    ju = n+1
    do while (ju-jl > 1) 
       jm = (ju+jl)/2
       if ((xx(n) > xx(1)) .eqv. (x > xx(jm))) then
          jl = jm
       else
          ju = jm
       endif
    enddo
    j = jl
    
    return
    
  end subroutine locatedb


  function path(pos,dir)

    ! compute distance to border of a cube (of side 1), from position
    ! pos (in cube units, i.e. [0,1]) and in direction dir (normalized). 
    
    implicit none

    real(kind=8),intent(in) :: pos(3)   ! position of photon in cell units
    real(kind=8),intent(in) :: dir(3)   ! propagation direction of photon
    integer(kind=4)         :: i
    real(kind=8)            :: dx(3)
    real(kind=8)            :: path     ! distance from pos to exit point

    do i = 1,3
       if(dir(i) < 0.0d0) then
          dx(i) = -pos(i) / dir(i)
       else if (dir(i) > 0.0d0) then
          dx(i) = (1.0d0 - pos(i)) / dir(i)
       else ! dir(i) == 0
          dx(i) = 10.0d0  ! larger than maximum extent of cell (sqrt(3)) in cell units
       end if
    end do
    path = minval(dx)

    return
    
  end function path


  subroutine print_rascas_header
    write(*,'(a)') '                                             '
    write(*,'(a)') '    _____  _____  _____  _____  _____  _____ '
    write(*,'(a)') '   / ___ \/ ___ \/ ____\/   __\/ ___ \/ ____\'
    write(*,'(a)') '   | \_/ || \_/ |\ \___ |  /   | \_/ |\ \___ '
    write(*,'(a)') '   |    _/|  _  | \___ \|  |   |  _  | \___ \'
    write(*,'(a)') '   | |\ \ | | | |_____/ |  \___| | | |____/ |'
    write(*,'(a)') '   |_| \_\|_| |_|\_____/\_____/|_| |_|\_____/'
    write(*,'(a)') '                                             '
    write(*,'(a)') '                                             '
  end subroutine print_rascas_header


end module module_utils
