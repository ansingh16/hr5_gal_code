module module_random

  implicit none

  ! Some variables useful for a "MPI ran3"
  integer(kind=4)           :: ix_ran=-1, iy_ran=-1
  real(kind=8)              :: am_ran
  integer(kind=4),parameter :: k4b=selected_int_kind(9)
  integer(kind=4),parameter :: ia_ran=16807,im_ran=2147483647,iq_ran=127773,ir_ran=2836

contains

!  function ran3(idum)
!    ! no MPI
!
!    implicit none
!
!    integer(kind=4),intent(inout) :: idum
!    real(kind=8)                  :: ran3
!    integer(kind=4),parameter     :: ia=16807,im=2147483647,iq=127773,ir=2836
!    real(kind=8),save             :: am
!    integer(kind=4),save          :: ix=-1, iy=-1, k
!
!    if(idum <=0 .or. iy < 0) then
!       am=nearest(1.0,-1.0)/im
!       iy=ior(ieor(888889999,abs(idum)),1)
!       ix=ieor(777755555,abs(idum))
!       idum=abs(idum)+1
!    end if
!
!    ix=ieor(ix,ishft(ix,13))
!    ix=ieor(ix,ishft(ix,-17))
!    ix=ieor(ix,ishft(ix,5))
!
!    k=iy/iq
!    iy=ia*(iy-k*iq)-ir*k
!    if(iy < 0) iy=iy+im
!
!    ran3=am*ior(iand(im,ieor(ix,iy)),1)
!
!  end function ran3


  function ran3(idum)
    ! A version of ran3 supposed to be used in the MPI version

    implicit none
  
    integer(kind=4),intent(inout) :: idum
    real(kind=8)                  :: ran3
    integer(kind=4),save          :: k
    
    if(idum <=0 .or. iy_ran < 0) then
       am_ran=nearest(1.0,-1.0)/im_ran
       iy_ran=ior(ieor(888889999,abs(idum)),1)
       ix_ran=ieor(777755555,abs(idum))
       idum=abs(idum)+1
    endif
  
    ix_ran=ieor(ix_ran,ishft(ix_ran,13))
    ix_ran=ieor(ix_ran,ishft(ix_ran,-17))
    ix_ran=ieor(ix_ran,ishft(ix_ran,5))
    
    k=iy_ran/iq_ran
    iy_ran=ia_ran*(iy_ran-k*iq_ran)-ir_ran*k
    if(iy_ran < 0) iy_ran=iy_ran+im_ran
    
    ran3=am_ran*ior(iand(im_ran,ieor(ix_ran,iy_ran)),1)
    
  end function ran3

end module module_random
