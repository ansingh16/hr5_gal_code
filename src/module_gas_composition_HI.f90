module module_gas_composition

  ! Pure HI gas. 
  ! - The HI content is derived from RAMSES directly

  use module_HI_1216_model
  use module_random
  use module_constants

  implicit none

  private

  character(100),parameter :: moduleName = 'module_gas_composition_HI.f90'
  
  type, public :: gas
     ! fluid
     real(kind=8) :: v(3)      ! gas velocity [cm/s]
     ! Hydrogen 
     real(kind=8) :: nHI       ! HI numerical density [HI/cm3]
     real(kind=8) :: temp  ! Doppler width [cm/s]
     real(kind=8) :: metal 
  end type gas
  real(kind=8),public :: box_size_cm   ! size of simulation box in cm. 

  ! --------------------------------------------------------------------------
  ! user-defined parameters - read from section [gas_composition] of the parameter file
  ! --------------------------------------------------------------------------
  ! possibility to overwrite ramses values with an ad-hoc model 
  logical                  :: gas_overwrite       = .false. ! if true, define cell values from following parameters 
  real(kind=8)             :: fix_nhi             = 0.0d0   ! ad-hoc HI density (H/cm3)
  real(kind=8)             :: fix_vth             = 1.0d5   ! ad-hoc thermal velocity (cm/s)
  real(kind=8)             :: fix_vel             = 0.0d0   ! ad-hoc cell velocity (cm/s) -> NEED BETTER PARAMETERIZATION for more than static... 
  real(kind=8)             :: fix_box_size_cm     = 1.0d8   ! ad-hoc box size in cm. 
  real(kind=8)             :: fix_ndust           = 0.0d0   ! ad-hoc dust number density (/cm3)

  ! --------------------------------------------------------------------------

  ! public functions:
  public :: gas_from_ramses_leaves,get_gas_velocity,dump_gas
  public :: read_gas,gas_destructor,read_gas_composition_params,print_gas_composition_params

contains
  

  subroutine gas_from_ramses_leaves(repository,snapnum,nleaf,nvar,ramses_var, g)

    ! define gas contents from ramses raw data

    use module_ramses

    character(2000),intent(in)                     :: repository 
    integer(kind=4),intent(in)                     :: snapnum
    integer(kind=4),intent(in)                     :: nleaf,nvar
    real(kind=8),intent(in),dimension(nvar,nleaf)  :: ramses_var
    type(gas),dimension(:),allocatable,intent(out) :: g
    integer(kind=4)                                :: ileaf
    real(kind=8),dimension(:),allocatable          :: T, nhi, metallicity
    real(kind=8),dimension(:,:),allocatable        :: v

    ! allocate gas-element array
    allocate(g(nleaf))

    if (gas_overwrite) then
       call overwrite_gas(g)
    else
       box_size_cm = ramses_get_box_size_cm(repository,snapnum)
       ! compute velocities in cm / s
       ! write(*,*) '-- module_gas_composition_HI : extracting velocities from ramses '
       allocate(v(3,nleaf))
       call ramses_get_velocity_cgs(repository,snapnum,nleaf,nvar,ramses_var,v)
       do ileaf = 1,nleaf
          g(ileaf)%v = v(:,ileaf)
       end do
       deallocate(v)
       ! get nHI and temperature from ramses
       ! write(*,*) '-- module_gas_composition_HI : extracting nHI and T from ramses '
       allocate(T(nleaf),nhi(nleaf))
       call ramses_get_T_nhi_cgs(repository,snapnum,nleaf,nvar,ramses_var,T,nhi)
       g(:)%nHI = nhi(:)
       ! compute thermal velocity 
       ! ++++++ TURBULENT VELOCITY >>>>> parameter to add and use here
       g(:)%temp = T 


       ! get ndust (pseudo dust density from Laursen, Sommer-Larsen, Andersen 2009)
       ! write(*,*) '-- module_gas_composition_HI_D_dust : extracting ndust from ramses '
       allocate(metallicity(nleaf))
       call ramses_get_metallicity(nleaf,nvar,ramses_var,metallicity)
       
       do ileaf = 1,nleaf
          g(ileaf)%metal = metallicity(ileaf) 
       end do

       deallocate(metallicity,T,nhi)

    end if

    return

  end subroutine gas_from_ramses_leaves

  

  subroutine overwrite_gas(g)
    ! overwrite ramses values with an ad-hoc model

    type(gas),dimension(:),intent(inout) :: g

    box_size_cm   = fix_box_size_cm
    
    g(:)%v(1)     = fix_vel
    g(:)%v(2)     = fix_vel
    g(:)%v(3)     = fix_vel
    g(:)%nHI      = fix_nhi
    g(:)%temp     = fix_vth
    g(:)%metal    = fix_ndust
  end subroutine overwrite_gas



  function get_gas_velocity(cell_gas)
    type(gas),intent(in)      :: cell_gas
    real(kind=8),dimension(3) :: get_gas_velocity
    get_gas_velocity(:) = cell_gas%v(:)
    return
  end function get_gas_velocity



  

  subroutine dump_gas(unit,g)
    type(gas),dimension(:),intent(in) :: g
    integer(kind=4),intent(in)        :: unit
    integer(kind=4)                   :: i,nleaf
    nleaf = size(g)
    write(unit) (g(i)%v(:), i=1,nleaf)
    write(unit) (g(i)%nHI, i=1,nleaf)
    write(unit) (g(i)%temp, i=1,nleaf)
    write(unit) (g(i)%metal, i=1,nleaf)

    ! write(*,*) g(1)
    
    write(unit) box_size_cm
  end subroutine dump_gas




  subroutine read_gas(unit,n,g)
    integer(kind=4),intent(in)                     :: unit,n
    type(gas),dimension(:),allocatable,intent(out) :: g
    integer(kind=4)                                :: i
    allocate(g(1:n))
    if (gas_overwrite) then
       call overwrite_gas(g)
    else
       read(unit) (g(i)%v(:),i=1,n)
       read(unit) (g(i)%nHI,i=1,n)
       read(unit) (g(i)%temp,i=1,n)
       read(unit) (g(i)%metal,i=1,n)
       read(unit) box_size_cm
    end if
  end subroutine read_gas

  

  subroutine gas_destructor(g)
    type(gas),dimension(:),allocatable,intent(inout) :: g
    deallocate(g)
  end subroutine gas_destructor


  subroutine read_gas_composition_params(pfile)

    ! ---------------------------------------------------------------------------------
    ! subroutine which reads parameters of current module in the parameter file pfile
    ! default parameter values are set at declaration (head of module)
    !
    ! ALSO read parameter form used modules (HI_model)
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
       if (line(1:17) == '[gas_composition]') then
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
          case ('gas_overwrite')
             read(value,*) gas_overwrite
          case ('fix_nhi')
             read(value,*) fix_nhi
          case ('fix_vth')
             read(value,*) fix_vth
          case ('fix_vel')
             read(value,*) fix_vel
          case ('fix_box_size_cm')
             read(value,*) fix_box_size_cm
          end select
       end do
    end if
    close(10)

    call read_HI_1216_params(pfile)

    return

  end subroutine read_gas_composition_params



  subroutine print_gas_composition_params(unit)

    ! ---------------------------------------------------------------------------------
    ! write parameter values to std output or to an open file if argument unit is
    ! present.
    ! ---------------------------------------------------------------------------------

    integer(kind=4),optional,intent(in) :: unit

    if (present(unit)) then 
       write(unit,'(a,a,a)') '[gas_composition]'
       write(unit,'(a,a)')     '# code compiled with: ',trim(moduleName)
       write(unit,'(a)')       '# overwrite parameters'
       write(unit,'(a,L1)')    '  gas_overwrite       = ',gas_overwrite
       if(gas_overwrite)then
          write(unit,'(a,ES10.3)') '  fix_nhi            = ',fix_nhi
          write(unit,'(a,ES10.3)') '  fix_vth            = ',fix_vth
          write(unit,'(a,ES10.3)') '  fix_vel            = ',fix_vel
          write(unit,'(a,ES10.3)') '  fix_box_size_cm    = ',fix_box_size_cm
       endif
       write(unit,'(a)')             ' '
       call print_HI_1216_params(unit)
    else
       write(*,'(a,a,a)') '[gas_composition]'
       write(*,'(a,a)')     '# code compiled with: ',trim(moduleName)
       write(*,'(a)')       '# overwrite parameters'
       write(*,'(a,L1)')    '  gas_overwrite       = ',gas_overwrite
       if(gas_overwrite)then
          write(*,'(a,ES10.3)') '  fix_nhi            = ',fix_nhi
          write(*,'(a,ES10.3)') '  fix_vth            = ',fix_vth
          write(*,'(a,ES10.3)') '  fix_vel            = ',fix_vel
          write(*,'(a,ES10.3)') '  fix_box_size_cm    = ',fix_box_size_cm
       endif
       write(*,'(a)')             ' '
       call print_HI_1216_params
    end if

    return

  end subroutine print_gas_composition_params



end module module_gas_composition
