module module_select

  use module_gas_composition

  implicit none

  interface select_from_domain
     module procedure xleaf_select, leaflevel_select, gas_select
  end interface select_from_domain

contains

  subroutine xleaf_select(arr_in, ind_sel, arr_out)

    real(kind=8),dimension(:,:),intent(in)              :: arr_in
    integer(kind=4),dimension(:),intent(in)             :: ind_sel
    real(kind=8),dimension(:,:),allocatable,intent(out) :: arr_out
    integer(kind=4)                                     ::i,n,ii

    n = size(ind_sel)
    
    allocate(arr_out(1:n,1:3))
    do i=1,n
       ii=ind_sel(i)
       arr_out(i,1:3)=arr_in(ii,1:3)
    enddo

  end subroutine xleaf_select


  subroutine leaflevel_select(arr_in, ind_sel, arr_out)

    integer(kind=4),dimension(:),intent(in)              :: arr_in
    integer(kind=4),dimension(:),intent(in)              :: ind_sel
    integer(kind=4),dimension(:),allocatable,intent(out) :: arr_out
    integer(kind=4)                                      ::i,n,ii

    n = size(ind_sel)
    
    allocate(arr_out(1:n))
    do i=1,n
       ii=ind_sel(i)
       arr_out(i)=arr_in(ii)
    enddo

  end subroutine leaflevel_select

  
  subroutine gas_select(arr_in, ind_sel, arr_out)

    type(gas),intent(in),dimension(:)              :: arr_in
    integer(kind=4),dimension(:),intent(in)        :: ind_sel
    type(gas),intent(out),dimension(:),allocatable :: arr_out
    integer(kind=4)                                ::i,n,ii

    n = size(ind_sel)
    
    allocate(arr_out(1:n))
    do i=1,n
       ii=ind_sel(i)
       arr_out(i) = arr_in(ii)
    enddo

  end subroutine gas_select


end module module_select
    
