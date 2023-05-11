

  subroutine friedman(O_mat_0,O_vac_0,O_k_0,alpha,axp_min,axp_out,hexp_out,tau_out,t_out,ntable,age_tot,dtau)
    
    integer,intent(in)::ntable
    real(kind=8),intent(in)::O_mat_0, O_vac_0, O_k_0
    real(kind=8),intent(in)::alpha,axp_min
    real(kind=8),intent(out)::age_tot,dtau
    real(kind=8),dimension(0:ntable),intent(out)::axp_out,hexp_out,tau_out,t_out
    real(kind=8)::axp_tau, axp_t
    real(kind=8)::axp_tau_pre, axp_t_pre
    real(kind=8)::dadtau, dadt
    real(kind=8)::dt
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
       axp_tau_pre = axp_tau - dadtau(axp_tau,O_mat_0,O_vac_0,O_k_0)*dtau/2d0
       axp_tau = axp_tau - dadtau(axp_tau_pre,O_mat_0,O_vac_0,O_k_0)*dtau
       tau = tau - dtau

       dt = alpha * axp_t / dadt(axp_t,O_mat_0,O_vac_0,O_k_0)
       axp_t_pre = axp_t - dadt(axp_t,O_mat_0,O_vac_0,O_k_0)*dt/2d0
       axp_t = axp_t - dadt(axp_t_pre,O_mat_0,O_vac_0,O_k_0)*dt
       t = t - dt

    end do

    age_tot=-t
    write(*,666)-t
666 format(' Age of the Universe (in unit of 1/H0)=',1pe10.3)

    nskip=nstep/ntable

    axp_t = 1d0
    t = 0d0
    axp_tau = 1d0
    tau = 0d0
    nstep = 0
    nout=0
    t_out(nout)=t
    tau_out(nout)=tau
    axp_out(nout)=axp_tau
    hexp_out(nout)=dadtau(axp_tau,O_mat_0,O_vac_0,O_k_0)/axp_tau

    do while ( (axp_tau .ge. axp_min) .or. (axp_t .ge. axp_min) )

       nstep = nstep + 1
       dtau = alpha * axp_tau / dadtau(axp_tau,O_mat_0,O_vac_0,O_k_0)
       axp_tau_pre = axp_tau - dadtau(axp_tau,O_mat_0,O_vac_0,O_k_0)*dtau/2d0
       axp_tau = axp_tau - dadtau(axp_tau_pre,O_mat_0,O_vac_0,O_k_0)*dtau
       tau = tau - dtau

       dt = alpha * axp_t / dadt(axp_t,O_mat_0,O_vac_0,O_k_0)
       axp_t_pre = axp_t - dadt(axp_t,O_mat_0,O_vac_0,O_k_0)*dt/2d0
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

  end subroutine friedman



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