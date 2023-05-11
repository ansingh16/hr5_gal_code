module module_voigt
  
  use module_constants, only : pi, sqrtpi
  
  private

  ! --------------------------------------------------------------------------
  ! user-defined parameters - read from section [voigt] in the parameter file 
  ! --------------------------------------------------------------------------
  character(20)   :: approximation = 'COLT'  ! could be 'COLT', 'Tasitsiomi' or 'Humlicek_w4'
  ! --------------------------------------------------------------------------

  integer(kind=4) :: approxKey
  logical         :: isRead=.False., isPrinted=.False. ! to avoid multiple reads and prints when called from different modules
  
  public :: voigt_function, read_voigt_params, print_voigt_params

contains

  function voigt_function(x,a)

    ! returns H(x,a) = a/pi * integral(exp(-z**2) dz / (a**2+(z-x)**2))
    
    implicit none
    
    real(kind=8),intent(in) :: x,a
    real(kind=8)            :: voigt_function

    !========================================
    ! if(approximation=='COLT')then
    !    voigt_function = colt_approx(x,a)
    ! else 
    !    if(approximation=='Tasitsiomi')then
    !       voigt_function = tasitsiomi_approx(x,a)
    !    else
    !       voigt_function = humlicek_w4(x,a)
    !    end if
    ! end if
    !========================================
    ! select case(trim(approximation))
    ! case('Tasitsiomi')
    !    voigt_function = tasitsiomi_approx(x,a)
    ! case('COLT')
    !    voigt_function = colt_approx(x,a)
    ! case('Humlicek_w4')
    !    voigt_function = humlicek_w4(x,a)
    ! end select
    !========================================
    select case(approxKey)
    case(1)
       voigt_function = tasitsiomi_approx(x,a)
    case(2)
       voigt_function = colt_approx(x,a)
    case(3)
       voigt_function = humlicek_w4(x,a)
    end select
    !========================================
    
    return

  end function voigt_function


  function tasitsiomi_approx(x,a)
    ! from Tasitsiomi 2006
    implicit none
    real(kind=8),intent(in) :: x,a
    real(kind=8)            :: tasitsiomi_approx
    real(kind=8)            :: q,z,x2
    x2 = x**2
    z  = (x2 - 0.855d0) / (x2 + 3.42d0)
    if (z > 0.0d0) then 
       q = z * (1.0d0 + 21.0d0/x2) * a / pi / (x2 + 1.0d0)
       q = q * (((5.674d0*z - 9.207d0)*z + 4.421d0)*z + 0.1117)
    else
       q = 0.0d0 
    end if
    tasitsiomi_approx = sqrtpi*q + exp(-x2)
    return
  end function tasitsiomi_approx


  function colt_approx(x,a)
    ! from Smith+15, appendix A1
    implicit none
    real(kind=8), intent(in) :: x,a
    real(kind=8)             :: A0,A1,A2,A3,A4,A5,A6
    real(kind=8)             :: B0,B1,B2,B3,B4,B5,B6,B7,B8
    real(kind=8)             :: z, colt_approx
    !--- coefficients for the rational function (Table A1, Smith+15) ---
    A0 = 15.75328153963877d0
    A1 = 286.9341762324778d0
    A2 = 19.05706700907019d0
    A3 = 28.22644017233441d0
    A4 = 9.526399802414186d0
    A5 = 35.29217026286130d0
    A6 = 0.8681020834678775d0
    B0 = 0.0003300469163682737d0
    B1 = 0.5403095364583999d0
    B2 = 2.676724102580895d0
    B3 = 12.82026082606220d0
    B4 = 3.21166435627278d0
    B5 = 32.032981933420d0
    B6 = 9.0328158696d0
    B7 = 23.7489999060d0
    B8 = 1.82106170570d0
    !------------------------------
    z = x*x
    if(z <= 3.0d0)then
       colt_approx = exp(-z) * (1.0d0-a*(A0+A1/(z-A2+A3/(z-A4+A5/(z-A6)))))
    else
       if(z < 25.0d0) then
          colt_approx = exp(-z) + a * (B0+B1/(z-B2+B3/(z+B4+B5/(z-B6+B7/(z-B8)))))
       else
          colt_approx = a/sqrtpi/(z-1.5d0-1.5d0/(z-3.5d0-5.0d0/(z-5.5d0)))
       end if
    end if
    return
  end function colt_approx


  function humlicek_w4(x,a)
    ! From Humlicek82
    real(kind=8), intent(in) :: x,a
    complex(kind=8)          :: t, u, w4
    real(kind=8)             :: s, absx, humlicek_w4, oneoversqrtpi
    
    t = cmplx(a,-x,8)
    u = t*t
    absx = abs(x)
    s = absx+abs(a)
    oneoversqrtpi = 1.d0/sqrtpi
    
    if(s >= 15.d0)then
       ! region 1
       w4 = t * oneoversqrtpi / (0.5d0 + u)
    else
       if(s >= 5.5)then
          ! region 2
          w4 = t*(1.410474d0+u*oneoversqrtpi)/(0.75d0+u*(3.d0+u))
       else
          if(abs(a) >= (0.195*absx-0.176)) then
             ! region 3
             w4 = (16.4955d0 + t*(20.20933d0 + t*(11.96482d0 + &
             t*(3.778987d0 + t*oneoversqrtpi)))) / &
            (16.4955d0 + t*(38.82363d0 + t*(39.27121d0 + & 
            t*(21.69274d0 + t*(6.699398d0+t)))))
          else
             ! region 4
             w4 = exp(u)-t*(36183.30536-u*(3321.990492-u*(1540.786893-u*(219.0312964- & 
             u*(35.76682780-u*(1.320521697-u*oneoversqrtpi)))))) / &
                  (32066.59372-u*(24322.84021-u*(9022.227659-u*(2186.181081-u*(364.2190727- &
                  u*(61.57036588-u*(1.841438936-u)))))))
          end if
       end if
    end if
    humlicek_w4 = real(w4)
    return
  end function humlicek_w4

  
  subroutine read_voigt_params(pfile)    
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
          if (line(1:7) == '[voigt]') then
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
             case ('approximation')
                write(approximation,'(a)') trim(value)
             end select
          end do
       end if
       close(10)
    endif
    isRead = .True.

    ! encode approximation choice into an integer
    select case(trim(approximation))
    case('Tasitsiomi')
       approxKey = 1
    case('COLT')
       approxKey = 2
    case('Humlicek_w4')
       approxKey = 3
    case default
       print*,'ERROR: approximation not known in module_voigt.f90: read_voigt_params: ',trim(approximation)
       stop
    end select

    return
  end subroutine read_voigt_params


  subroutine print_voigt_params(unit)
    ! ---------------------------------------------------------------------------------
    ! write parameter values to std output or to an open file if argument unit is
    ! present.
    ! ---------------------------------------------------------------------------------
    integer(kind=4),optional,intent(in) :: unit
    if(.not.(isPrinted)) then
       if (present(unit)) then
          !write(unit,'(a)') ''
          write(unit,'(a,a,a)')    '[voigt]'
          write(unit,'(a,a)')      '  approximation = ',approximation
          !write(unit,'(a)') ''
       else
          !write(*,*) ''
          write(*,'(a,a,a)')    '[voigt]'
          write(*,'(a,a)')      '  approximation = ',approximation
          !write(*,*) ''
       end if
    end if
    isPrinted = .True.
    return
  end subroutine print_voigt_params
  
end module module_voigt
