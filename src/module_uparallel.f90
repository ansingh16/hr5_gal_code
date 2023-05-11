module module_uparallel

  use module_constants, only:pi
  use module_random

  private

  ! --------------------------------------------------------------------------
  ! user-defined parameters - read from section [uparallel] in the parameter file 
  ! --------------------------------------------------------------------------
  character(20)   :: method = 'RASCAS'   ! may be 'Smith', 'Semelin', or 'RASCAS'
  real(kind=8)    :: xForGaussian = 8.0  ! above this value, use a Gaussian to draw u_parallel
  ! --------------------------------------------------------------------------

  logical         :: isRead=.False., isPrinted=.False. ! to avoid multiple reads and prints when called from different modules
  integer(kind=4) :: methodKey
  
  public :: get_uparallel,read_uparallel_params,print_uparallel_params

contains
  
  function get_uparallel(y,a,iran)
    
    ! --------------------------------------------------------------------------------
    ! Draw a realisation from the probability function p(u)=exp(-u*u)/((x-u)**2+a**2)`
    ! Method from Zheng & Miralda-Escude (2002), modified to incorporate variants by
    ! Semelin+02 (who kindly provided the first version of this routine), Smith+15, and
    ! out own RASCAS version (Michel-Dansac+18). 
    !
    ! PARAMETERS 
    ! - y    : this is the frequency of the photon in doppler units (i.e. y is x)
    ! - a    : this is a
    ! - iran : current seed for random number generator 
    ! --------------------------------------------------------------------------------
    
    implicit none
    
    real(kind=8),intent(in)       :: y,a
    integer(kind=4),intent(inout) :: iran
    real(kind=8)                  :: get_uparallel
    real(kind=8)                  :: u,coeff,u_0, signe,x,r,p,theta1,xcw,la,la2
    logical                       :: success
    real(kind=8),parameter        :: pi2=pi/2.
    integer(kind=4)               :: ctloc

    if(a < 1.d-6) then
       print*,'Error using module_uparallel.f90: get_uparallel: a too small'
       STOP
    endif

    signe=sign(dble(1),y)
    x=y*signe
    
    if( x < xForGaussian ) then
       ! Define a value of u_0 (the optimisation parameter from Zheng & Miralda-Escude 2002). 

       ! select case(trim(method))
       ! case ('Semelin')
       !    if (x>3) then 
       !       u_0=1.85-log(a)/6.73+log(log(x))  ! Eq. 17 of Semelin+07
       !    else
       !       u_0 = 0.0d0
       !    end if
       ! case ('Smith')
       !    xcw = 6.9184721d0 + 81.766279d0 / (log10(a) - 14.651253d0)  ! Eq. 21 of Smith+15
       !    if (x < xcw) then
       !       u_0 = x - 1.0d0 / (x + exp(1.d0-x*x)/a)   ! Eq. 31 of Smith+15
       !    else
       !       u_0 = xcw - 1.d0/xcw + 0.15*(x-xcw)       ! Eq. 32 of Smith+15
       !    end if
       ! case ('RASCAS')
       !    ! empirical 2D polynomial for u0 (JB-2017)
       !    if (x<0.6) then 
       !       u_0 = 0
       !    else
       !       la = log10(a)
       !       la2 = la*la
       !       u_0 =2.648963+2.014446*la+0.351479*la2 + x*(-4.058673-3.675859*la-0.640003*la2 &
       !            + x*(3.017395+2.117133*la+0.370294*la2 + x*(-0.869789-0.565886*la-0.096312*la2 &
       !            + x*(0.110987+0.070103*la+0.011557*la2 + x*(-0.005200-0.003240*la-0.000519*la2)))))
       !    end if
       ! case default
       !    print*,'ERROR: method not known in module_uparallel.f90:PROB_FUNC : ',trim(method)
       !    stop
       ! end select

       select case(methodKey)
       case (1)
          if (x>3) then 
             u_0=1.85-log(a)/6.73+log(log(x))  ! Eq. 17 of Semelin+07
          else
             u_0 = 0.0d0
          end if
       case (2)
          xcw = 6.9184721d0 + 81.766279d0 / (log10(a) - 14.651253d0)  ! Eq. 21 of Smith+15
          if (x < xcw) then
             u_0 = x - 1.0d0 / (x + exp(1.d0-x*x)/a)   ! Eq. 31 of Smith+15
          else
             u_0 = xcw - 1.d0/xcw + 0.15*(x-xcw)       ! Eq. 32 of Smith+15
          end if
       case (3)
          ! empirical 2D polynomial for u0 (JB-2017)
          if (x<0.6) then 
             u_0 = 0
          else
             la = log10(a)
             la2 = la*la
             u_0 =2.648963+2.014446*la+0.351479*la2 + x*(-4.058673-3.675859*la-0.640003*la2 &
                  + x*(3.017395+2.117133*la+0.370294*la2 + x*(-0.869789-0.565886*la-0.096312*la2 &
                  + x*(0.110987+0.070103*la+0.011557*la2 + x*(-0.005200-0.003240*la-0.000519*la2)))))
          end if
       case default
          print*,'ERROR: method not known in module_uparallel.f90: get_uparallel: ',trim(method)
          stop
       end select
       
       ! Perform the rejection method given u_0
       ctloc=0
       success=.false.
       coeff=exp(-u_0*u_0)
       theta1=atan((u_0-x)/a)
       p=(theta1+pi2)/((1.-coeff)*theta1+(1+coeff)*pi2)
       do while(.not.success)
          r = ran3(iran)
          if(r<p) then
             r = ran3(iran)
             r=r*(theta1+pi2)-pi2
             u=a*tan(r)+x
             r = ran3(iran)
             if(r < exp(-u*u) ) then
                success=.true.
             endif
          else
             r = ran3(iran)
             r=r*(-theta1+pi2)+theta1
             u=a*tan(r)+x
             r = ran3(iran)
             if(r < exp(-u*u+u_0*u_0) ) then
                success=.true.
             endif
          endif
          ctloc=ctloc+1
          if(mod(ctloc,10000000)==0) print*,'rrr upar...',ctloc,r,u,u_0,a,x
       enddo

    else  ! For large x's, simply draw from a Gaussian

       u = upar_from_gaussian(x,iran)

    end if

    get_uparallel=u*signe

    return
    
  end function get_uparallel

  

  function upar_from_gaussian(xin,iran)
    
    implicit none
    
    real(kind=8),intent(in)       :: xin
    integer(kind=4),intent(inout) :: iran
    real(kind=8)                  :: r1,r2
    real(kind=8)                  :: u_mean
    real(kind=8)                  :: upar_from_gaussian

    r1 = ran3(iran)
    r2 = ran3(iran)
    
    u_mean = 1.0d0 / xin
    upar_from_gaussian = sqrt(-2.*log(r1)) * cos(2.0d0*pi*r2)
    upar_from_gaussian = upar_from_gaussian * 1.0d0 / sqrt(2.0d0) + u_mean ! width is 1 since u = x = v / vth
    
    return
    
  end function upar_from_gaussian


  subroutine read_uparallel_params(pfile)
    
    ! ---------------------------------------------------------------------------------
    ! subroutine which reads parameters of current module in the parameter file pfile
    !
    ! default parameter values are set at declaration (head of module)
    ! ---------------------------------------------------------------------------------

    character(*),intent(in) :: pfile
    character(1000)         :: line,name,value
    integer(kind=4)         :: err,i
    logical                 :: section_present
    
    if(.not.(isRead)) then
       section_present = .false.
       open(unit=10,file=trim(pfile),status='old',form='formatted')
       ! search for section start
       do
          read (10,'(a)',iostat=err) line
          if(err/=0) exit
          if (line(1:11) == '[uparallel]') then
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
             case ('method')
                write(method,'(a)') trim(value)
             case ('xForGaussian')
                read(value,*) xForGaussian
             end select
          end do
       end if
       close(10)
    endif

    select case(trim(method))
    case('Semelin')
       methodKey = 1
    case('Smith')
       methodKey = 2
    case('RASCAS')
       methodKey = 3
    case default
       print*,'ERROR: method not known in module_uparallel.f90: read_uparallel_params: ',trim(method)
       stop
    end select
    
    isRead = .True.
    
    return
    
  end subroutine read_uparallel_params


  
  subroutine print_uparallel_params(unit)
    
    ! ---------------------------------------------------------------------------------
    ! write parameter values to std output or to an open file if argument unit is
    ! present.
    ! ---------------------------------------------------------------------------------

    integer(kind=4),optional,intent(in) :: unit

    if(.not.(isPrinted)) then
       if (present(unit)) then
          !write(unit,'(a)') ''
          write(unit,'(a,a,a)')    '[uparallel]'
          write(unit,'(a,a)')      '  method       = ',method
          write(unit,'(a,ES10.3)') '  xForGaussian = ',xForGaussian
          write(unit,'(a)') ''
       else
          !write(*,*) ''
          write(*,'(a,a,a)')    '[uparallel]'
          write(*,'(a,a)')      '  method       = ',method
          write(*,'(a,ES10.3)') '  xForGaussian = ',xForGaussian
          write(*,*) ''
       end if
    end if

    isPrinted = .True.
    
    return
    
  end subroutine print_uparallel_params

  
end module module_uparallel



