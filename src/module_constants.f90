module module_constants
  
  public

  ! from the NIST Constant Index
  real(kind=8),parameter :: clight  = 2.99792458d+10          ![cm/s] light speed
  real(kind=8),parameter :: cmtoA   = 1.0d8                   ! from cm to A
  real(kind=8),parameter :: evtoerg = 1.6021766d-12           ! from eV to erg
  real(kind=8),parameter :: kb      = 1.38064852d-16          ![erg/K] Boltzman constant
  real(kind=8),parameter :: mp      = 1.672621898d-24         ![g] proton mass
  real(kind=8),parameter :: mpc     = 3.08d24                 ![cm] from Mpc to cm
  real(kind=8),parameter :: me      = 9.10938356d-28          ![g] electron 
  real(kind=8),parameter :: pi      = 3.1415926535898d0       ! pi
  real(kind=8),parameter :: twopi   = 2.0d0 * pi              ! 2 x pi
  real(kind=8),parameter :: fourpi  = 4.0d0 * pi              ! 4 x pi 
  real(kind=8),parameter :: sqrtpi  = 1.77245387556703d0      ! sqrt(pi)
  real(kind=8),parameter :: grtoev  = 1.782661d-33*clight**2  ! from gr to eV
  real(kind=8),parameter :: e_ch    = 4.80320451d-10          ![esu] electron charge
  real(kind=8),parameter :: planck  = 6.626070040d-27         ![erg s] Planck's constant
  real(kind=8),parameter :: sqrt_H2Deut_mass_ratio = 0.7071067811865d0   ! == sqrt(mp/mdeut) = 1/sqrt(2).
  real(kind=8),parameter :: XH = 0.76d0
  real(kind=8),parameter :: amu = 1.66054d-24                 ![g] atomic mass unit
  real(kind=8),parameter :: msun = 1.989d33                   ![g] solar mass 
  real(kind=8),parameter :: mSi = 28.085d0 * amu              ![g] mass of Silicon
  real(kind=8),parameter :: mMg = 24.305d0 * amu              ![g] mass of Magnesium
  real(kind=8),parameter :: mFe = 55.845d0 * amu              ![g] mass of Fer

end module module_constants

