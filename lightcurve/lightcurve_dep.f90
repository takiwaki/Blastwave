program main

  implicit none

  real(8), parameter :: day_in_second = 8.64d4
  real(8), parameter :: Msun = 1.9891d33 ! Solar mass in g
  real(8), parameter :: clight = 2.99792458d10 ! speed of light in cm/s
  real(8) :: pi
  
  real(8) :: Mej,Ek,Mni,kappa,v_ej,t_diff
  real(8) :: Lsn,t,h,k1,k2,k3,k4

  open(1,file='lightcurve_dep.dat',status='new')
  ! status = 'new' : forbid to overwrite an existing file
  ! status = 'unknown' : allow to overwrite an existing file
    
  ! ejecta mass in Msun
  Mej = 1.4d0

  ! kinetic energy in 1e51 erg
  Ek  = 1.3d0

  ! 56Ni mass in Msun
  Mni = 0.6d0
  
  ! ejecta opacity in cm2/g
  kappa = 0.1d0

  write(1,'("# light curve with Mej =",x,f6.3,x,"Msun, Ek =",x,f6.3,x,"x 10^51 erg, and Mni =",x,f6.3,x,"Msun")') Mej, Ek, Mni
  write(1,'(a)') '# time (day), luminosity (erg/s)'
  
  Mej = Mej*Msun
  Ek  = Ek*1d51

  pi = 2d0*asin(1d0)
  
  v_ej = sqrt(2d0*Ek/Mej)
  t_diff = sqrt(3d0*Mej*kappa/(4d0*pi*v_ej*clight))
  
  Lsn = 0d0
  t = 1d3
  h = 1d3
  
  do while (t < 100d0*day_in_second) ! light-curve calculation up to 100 days
     
     k1 = h*f(t,Lsn)
     k2 = h*f(t+h/2d0,Lsn+k1/2d0)
     k3 = h*f(t+h/2d0,Lsn+k2/2d0)
     k4 = h*f(t+h,Lsn+k3)

     Lsn = Lsn + (k1+2d0*(k2+k3)+k4)/6d0

     write(1,'(f9.5,x,1pe13.6)') t/day_in_second,Lsn

     t = t+h
     
  end do

  close(1)
  
  stop
  
contains

  function f(x,y)

    implicit none

    real(8), intent(in)  :: x,y
    real(8) :: f

    real(8), parameter :: tau_56Ni = 8.8d0   ! 56Ni decay time in days
    real(8), parameter :: tau_56Co = 111.3d0 ! 56Co decay time in days
    real(8) :: Qdot,kappa_gamma,tau_gamma,f_dep
    
    Qdot = (6.45d43*exp(-x/day_in_second/tau_56Ni) &
         +1.45d43*exp(-x/day_in_second/tau_56Co))*Mni

    kappa_gamma = 0.03d0
    tau_gamma = 3d0*kappa_gamma*Mej/(4d0*pi*v_ej**2d0*x**2d0)
    f_dep = 1d0 - exp(-tau_gamma)
    
    f = x/t_diff**2d0*(Qdot*f_dep-y)

  end function f
  
end program main
