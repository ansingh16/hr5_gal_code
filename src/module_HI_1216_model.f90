module module_HI_1216_model

  use module_constants
  use module_utils, only : isotropic_direction, anisotropic_direction_HIcore, anisotropic_direction_Rayleigh
  use module_uparallel
  use module_random
  use module_voigt

  implicit none

  private

  ! definition of atomic values
  real(kind=8),parameter,public :: lambda_0=1215.67d0                           ![A] Lya wavelength
  real(kind=8),parameter        :: gamma=6.265d8                                ! Einstein coeff = damping constant for Voigt Function(gamma_alpha)  
  real(kind=8),parameter        :: f12=0.416d0                                  ! oscillator strength for Ly-alpha
  ! useful pre-computed quantities
  real(kind=8),parameter,public :: lambda_0_cm = lambda_0 / cmtoA               ! cm
  real(kind=8),parameter,public :: nu_0 = clight / lambda_0_cm                  ! Hz
  real(kind=8),parameter        :: sigmaH_factor = sqrtpi*e_ch**2*f12/me/clight ! H cross-section factor-> multiply by Voigt(x,a)/delta_nu_doppler to get sigma.
  real(kind=8),parameter,public :: gamma_over_fourpi = gamma / fourpi

  ! --------------------------------------------------------------------------
  ! user-defined parameters - read from section [HI] in the parameter file 
  ! --------------------------------------------------------------------------
  logical                       :: recoil       = .true.      ! if set to true, recoil effect is computed [default is true]
  logical                       :: isotropic    = .false.     ! if set to true, scattering events will be isotropic [default is false]
  !--CORESKIP--
  logical                       :: HI_core_skip = .false.     ! if true, skip scatterings in the core of the line (as in Smith+15).
  real(kind=8)                  :: xcritmax     = -1d10        ! core-skipping will truncate at min(xcrit, xcritmax) -> set to a low (but positive) value to activate. 
  !--PIKSEROC--
  ! --------------------------------------------------------------------------

  public :: get_tau_HI_1216, scatter_HI_1216, read_HI_1216_params, print_HI_1216_params
  !--CORESKIP--
  public :: HI_core_skip 
  !--PIKSEROC-- 
  
contains

  function get_tau_HI_1216(nhi, vth, distance_to_border_cm, nu_cell)

    ! --------------------------------------------------------------------------
    ! compute optical depth of Hydrogen over a given distance
    ! --------------------------------------------------------------------------
    ! INPUTS:
    ! - nhi      : number density of neutral HI atoms                      [ cm^-3 ]
    ! - vth      : thermal (+ small-scale turbulence) velocity of HI atoms [ cm / s ]
    ! - distance_to_border_cm : distance over which we compute tau        [ cm ]
    ! - nu_cell  : photon's frequency in the frame of the cell            [ Hz ]
    ! OUTPUT :
    ! - get_tau_HI_1216 : optical depth of Hydrogen's Lya line over distance_to_border_cm
    ! --------------------------------------------------------------------------
    
    real(kind=8),intent(in) :: nhi,vth,distance_to_border_cm,nu_cell
    real(kind=8)            :: delta_nu_doppler,x_cell,sigmaH,a,h_cell, get_tau_HI_1216

    ! compute Doppler width and a-parameter, for H 
    delta_nu_doppler = vth / lambda_0_cm 
    a = gamma_over_fourpi / delta_nu_doppler
 
    ! Cross section of H 
    x_cell = (nu_cell - nu_0)/delta_nu_doppler
    h_cell = voigt_function(x_cell,a)
    sigmaH = sigmaH_factor / delta_nu_doppler * h_cell

    get_tau_HI_1216 = sigmaH * nhi * distance_to_border_cm

    return

  end function get_tau_HI_1216


  
  !--CORESKIP--
  subroutine scatter_HI_1216(vcell,vth,nu_cell,k,nu_ext,iran,xcrit)
  !subroutine scatter_HI_1216(vcell,vth,nu_cell,k,nu_ext,iran)
  !--PIKSEROC--
    ! ---------------------------------------------------------------------------------
    ! perform scattering event on a Hydrogen atom with an anisotropic phase function
    ! ---------------------------------------------------------------------------------
    ! INPUTS :
    ! - vcell    : bulk velocity of the gas (i.e. cell velocity)       [ cm / s ] 
    ! - vth      : thermal (+turbulent) velocity dispersion of H atoms [ cm / s ] 
    ! - nu_cell  : frequency of incoming photon in cell's rest-frame   [ Hz ] 
    ! - k        : propagaction vector (normalized) 
    ! - nu_ext   : frequency of incoming photon, in external frame     [ Hz ]
    ! - iran     : random number generator seed
    ! OUTPUTS :
    ! - nu_cell  : updated frequency in cell's frame   [ Hz ]
    ! - nu_ext   : updated frequency in external frame [ Hz ]
    ! - k        : updated propagation direction
    ! _ iran     : updated value of seed
    ! ---------------------------------------------------------------------------------
    !
    ! Notes on the phase function :
    ! -----------------------------
    ! - for core photons (|x| < 0.2) we use P(mu) = 11/24 + 3/24 * mu**2
    ! - for wing photons (|x| > 0.2) we use P(mu) = 3/8 * (1 + mu**2) [this is Rayleigh]
    ! where mu = cos(theta), (and theta in [0,pi]).
    ! ---------------------------------------------------------------------------------

    real(kind=8),intent(inout)              :: nu_cell, nu_ext
    real(kind=8),dimension(3),intent(inout) :: k
    real(kind=8),dimension(3),intent(in)    :: vcell
    real(kind=8),intent(in)                 :: vth
    !--CORESKIP--
    real(kind=8),intent(in)                 :: xcrit
    real(kind=8)                            :: xc
    !--PIKSEROC--
    integer(kind=4),intent(inout)           :: iran
    real(kind=8)                            :: delta_nu_doppler, a, x_cell, upar, ruper
    real(kind=8)                            :: r2, uper, nu_atom, mu, bu, scalar
    real(kind=8)                            :: x_atom
    real(kind=8),dimension(3)               :: knew

    !--CORESKIP--  sanity check ... 
    if (.not. HI_core_skip .and. xcrit .ne. 0.0d0) then
       print*,'ERROR: core skipping is off but xcrit is not zero ... '
       stop
    end if
    if (HI_core_skip)  then
       xc = min(xcrit,xcritmax)
    else
       xc=0.0d0
    endif
    !--PIKSEROC-- 
    
    ! define x_cell & a
    delta_nu_doppler = vth / lambda_0_cm 
    a = gamma_over_fourpi / delta_nu_doppler
    x_cell = (nu_cell - nu_0) / delta_nu_doppler

    ! 1/ component parallel to photon's propagation
    ! -> get velocity of interacting atom parallel to propagation
    upar = get_uparallel(x_cell,a,iran)
    upar = upar * vth    ! upar is an x -> convert to a velocity 

    ! 2/ component perpendicular to photon's propagation
    ruper  = ran3(iran)
    r2     = ran3(iran)
    !--CORESKIP--
    uper   = sqrt(xc**2-log(ruper))*cos(twopi*r2)
    !uper   = sqrt(-log(ruper))*cos(twopi*r2)
    !--PIKSEROC--
    uper   = uper * vth  ! from x to velocity

    ! 3/ incoming frequency in atom's frame
    nu_atom = nu_cell - nu_ext * upar/clight

    ! 4/ determine direction of scattered photon
    if (isotropic) then
       call isotropic_direction(knew,iran)
       mu = k(1)*knew(1) + k(2)*knew(2) + k(3)*knew(3)
       bu = sqrt(1.0d0 - mu*mu)
    else
       x_atom  = (nu_atom -nu_0) / delta_nu_doppler
       if (abs(x_atom) < 0.2d0) then ! core scattering 
          call anisotropic_direction_HIcore(k,knew,mu,bu,iran)
       else ! wing scattering 
          call anisotropic_direction_Rayleigh(k,knew,mu,bu,iran)
       end if
    end if

    ! 5/ recoil effect 
    if (recoil) then 
       nu_atom = nu_atom / (1.0d0 + ((planck*nu_atom)/(mp*clight*clight))*(1.0d0-mu))
    end if
    
    ! 6/ compute atom freq. in external frame, after scattering
    scalar = knew(1) * vcell(1) + knew(2) * vcell(2) + knew(3)* vcell(3)
    nu_ext = nu_atom * (1.0d0 + scalar/clight + (upar*mu + bu*uper)/clight)
    nu_cell = (1.0d0 - scalar/clight) * nu_ext 
    k = knew

  end subroutine scatter_HI_1216

  

  subroutine read_HI_1216_params(pfile)
    
    ! ---------------------------------------------------------------------------------
    ! subroutine which reads parameters of current module in the parameter file pfile
    !
    ! default parameter values are set at declaration (head of module)
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
       if (line(1:4) == '[HI]') then
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
          case ('recoil')
             read(value,*) recoil
          case ('isotropic')
             read(value,*) isotropic
          !--CORESKIP--
          case ('HI_core_skip') 
             read(value,*) HI_core_skip
          case ('xcritmax')
             read(value,*) xcritmax
          !--PIKSEROC--
          end select
       end do
    end if
    close(10)

    !--CORESKIP--  sanity check ... 
    if (HI_core_skip .and. xcritmax <= 0.0d0) then
       print*,'ERROR: core skipping is on but xcritmax is not set... '
       stop
    end if
    !--PIKSEROC-- 
    
    call read_uparallel_params(pfile)
    call read_voigt_params(pfile)

    return

  end subroutine read_HI_1216_params


  
  subroutine print_HI_1216_params(unit)
    
    ! ---------------------------------------------------------------------------------
    ! write parameter values to std output or to an open file if argument unit is
    ! present.
    ! ---------------------------------------------------------------------------------

    integer(kind=4),optional,intent(in) :: unit

    if (present(unit)) then 
       write(unit,'(a,a,a)') '[HI]'
       write(unit,'(a,L1)') '  recoil    = ',recoil
       write(unit,'(a,L1)') '  isotropic = ',isotropic
       !--CORESKIP--
       write(unit,'(a,L1)')     '  HI_core_skip = ',HI_core_skip
       write(unit,'(a,ES10.3)') '  xcritmax     = ',xcritmax
       !--PIKSEROC--
       write(unit,'(a)') ''
       call print_uparallel_params(unit)
       call print_voigt_params(unit)
    else
       write(*,'(a,a,a)') '[HI]'
       write(*,'(a,L1)') '  recoil    = ',recoil
       write(*,'(a,L1)') '  isotropic = ',isotropic
       !--CORESKIP--
       write(*,'(a,L1)')     '  HI_core_skip = ',HI_core_skip
       write(*,'(a,ES10.3)') '  xcritmax     = ',xcritmax
       !--PIKSEROC--
       write(*,'(a)') ''
       call print_uparallel_params()
       call print_voigt_params()
    end if
    
    return
    
  end subroutine print_HI_1216_params


end module module_HI_1216_model
