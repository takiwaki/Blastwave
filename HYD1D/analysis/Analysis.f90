module unitsmod
  implicit none
  real(8),parameter::    pc  = 3.085677581d18   ! parsec in [cm]
  real(8),parameter::    mu  = 1.660539066d-24  ! g
  real(8),parameter:: Msolar = 1.989e33         ! g
  real(8),parameter::   kbol = 1.380649d-23     ! J/K
  real(8),parameter::   year = 365.0d0*24*60*60 ! sec  

  real(8),parameter:: mp=1.67262192369d-24  !! proton mass [g]
  real(8),parameter:: kB=1.380649d-16  !! Boltzmann constant [erg/K]
  real(8),parameter:: erg_to_keV= 6.242d8  !! erg => keV
end module unitsmod
      

module fieldmod
    implicit none
    integer:: incr
    real(8):: time,dt
    integer:: in,jn,kn
    integer:: izone,jzone,kzone
    integer:: igs,jgs,kgs
    integer:: is,js,ks,ie,je,ke
    real(8),dimension(:),allocatable:: x1b,x2b,dvl1a
    real(8),dimension(:),allocatable:: x1a,x2a
    real(8),dimension(:,:,:),allocatable:: d,v1,v2,v3,p,ei,gp
    real(8),dimension(:,:,:),allocatable:: Tem,edot
    real(8):: dx
    real(8):: gam,rho0,Eexp

    real(8):: rshock,Msw,kTshock,Vshock,Lbol

end module fieldmod

program data_analysis
  use fieldmod
  implicit none
  integer:: fbeg, fend
  logical:: flag
  integer,parameter:: unitcon=100

  INQUIRE(FILE ="control.dat",EXIST = flag)
  if(flag) then
     open (unitcon,file="control.dat" &
     &        ,status='old',form='formatted')
     read (unitcon,*) fbeg,fend
     close(unitcon)
  endif

  FILENUMBER: do incr  = fbeg,fend
     write(6,*) "file index",incr
     call ReadData
     call EstimateEmissivity
     call FindShockRadiusAnswer
     call Visualize1D
     call TimeProfleAnswer
  enddo FILENUMBER

  stop
end program data_analysis

subroutine ReadData
  use fieldmod
  implicit none   
  character(20),parameter::dirname="../bindata/"
  character(40)::filename
  integer,parameter::unitinp=13
  integer,parameter::unitbin=14
  character(8)::dummy
  logical flag
  logical,save:: is_inited
  data is_inited / .false. /

  write(filename,'(a3,i5.5,a4)')"unf",incr,".dat"
  filename = trim(dirname)//filename

  INQUIRE(FILE =filename,EXIST = flag)
  if(.not. flag) then
     write(6,*) "FILE:",filename
     stop 'Cannot Open  data'
  endif
  open(unitinp,file=filename,status='old',form='formatted')
  read(unitinp,*) dummy,time,dt
  read(unitinp,*) dummy,izone,igs
  close(unitinp)
  in=izone+2*igs
  jn=1
  kn=1

  is=1+igs
  js=1
  ks=1
  ie=in-igs
  je=1
  ke=1

  if(.not. is_inited)then
     allocate( x1b(in),x1a(in),dvl1a(in))
     allocate( d(in,jn,kn))
     allocate(v1(in,jn,kn))
     allocate(v2(in,jn,kn))
     allocate(v3(in,jn,kn))
     allocate( p(in,jn,kn))
     allocate(ei(in,jn,kn))
     is_inited = .true.
  endif

  write(filename,'(a3,i5.5,a4)')"bin",incr,".dat"
  filename = trim(dirname)//filename
  open(unitbin,file=filename,status='old',form='unformatted',access="stream")
  read(unitbin)x1b(:),x1a(:),dvl1a(:)
!  read(unitbin)x2b(:),x2a(:)
  read(unitbin)  d(:,:,:)
  read(unitbin) v1(:,:,:)
  read(unitbin) v2(:,:,:)
  read(unitbin) v3(:,:,:)
  read(unitbin)  p(:,:,:)
  read(unitbin) ei(:,:,:)
  close(unitbin)
  
  return
end subroutine ReadData
subroutine EstimateEmissivity
  use unitsmod
  use fieldmod
  implicit none
  integer::i,j,k
  logical:: is_inited
  data is_inited / .false. /
  if(.not. is_inited) then
     allocate(Tem(in,jn,kn))
     allocate(edot(in,jn,kn))
     is_inited = .true.
  endif
  k=ks
  j=js
  do i =is,ie
     ! p = n k T => T = p/(n)/k since kbol[J/K] kbol*1.0d5 [erg/K] 
     Tem(i,j,k) = p(i,j,k)/(d(i,j,k)/mp) / (kbol*1.0d5) ! [K]
     edot(i,j,k) = 1.4d-27 * (d(i,j,k)/mp)**2 *sqrt(Tem(i,j,k)) ! erg/s/cm^3
  enddo
end subroutine EstimateEmissivity
  
subroutine FindShockRadius
  use unitsmod
  use fieldmod
  implicit none
  integer::i,j,k
  
  k = ks
  j = js 
  do i=is,ie
     ! find pressure max and set rshock, note x1b(i) is the radius
     rshock = 0.0d0
     ! use p = n T and set T_shock 
     kTshock = 0.0d0! T [keV]
     Vshock = 0.0d0 ! v [km/s]
  enddo
  !print *, "rshock=",rshock/pc,"[pc]"
  
end subroutine FindShockRadius


subroutine Visualize1D
  use unitsmod
  use fieldmod
  implicit none
  integer::i,j,k

  character(20),parameter::dirname="output/"
  character(40)::filename
  integer,parameter::unit1D=123

  logical,save:: is_inited
  data is_inited / .false. /

  if(.not. is_inited)then
     call makedirs(dirname)
     is_inited = .true.
  endif

  write(filename,'(a6,i5.5,a4)')"onepro",incr,".dat"
  filename = trim(dirname)//filename
  open(unit1D,file=filename,status='replace',form='formatted')

  write(unit1D,'(1a,a7,1(1x,E12.3))') "#","  time=",time
!                                     1234567890123   1234567890123   1234567890123   1234567890123
  write(unit1D,'(1a,4(1x,a13))') "#","1:r[cm]      ","2:den[1/cm^3]","3:p[erg/cm3] ","4:vel[cm/s]  "
  k=ks
  j=js
  do i=is,ie
     write(unit1D,'(1x,SP,4(1x,E13.3))') x1b(i),d(i,j,k)/mu,p(i,j,k),v1(i,j,k)
  enddo
  close(unit1D)

  return
end subroutine Visualize1D

subroutine TimeProfle
  use unitsmod
  use fieldmod
  implicit none
  integer::i,j,k

  character(20),parameter::dirname="output/"
  character(40)::filename
  integer::unittpr
  real(8)::Etot,pi

  logical,save:: is_inited
  data is_inited / .false. /

  if(.not. is_inited)then
     call makedirs(dirname)
  endif

  pi = acos(-1.0d0)
  Msw = 0.0d0
  Etot=0.0d0
  k=ks
  j=js
  do i=is,ie
     ! add Msw if possible
     Etot = Etot + (0.5d0*d(i,j,k)*v1(i,j,k)**2+ei(i,j,k))*dvl1a(i)*4.0d0*pi
  enddo
  !print *, "Msw=",Msw/Msolar,"[M_s]"
  write(filename,'(a3,i5.5,a4)')"tpr",incr,".dat"
  filename = trim(dirname)//filename
  open(newunit=unittpr,file=filename,status='replace',form='formatted')
  if(.not. is_inited) write(unittpr,'(1a,1x,A)') "#"," time[year] rshock[pc] Msw[Ms] kTshock[keV] Etot[erg]"

  write(unittpr,'(1x,5(1x,E13.3))') time/year,rshock/pc,Msw/Msolar,kTshock,Etot
  close(unittpr)
  
  is_inited = .true.

  return
end subroutine  TimeProfle

subroutine makedirs(outdir)
  implicit none
  character(len=*), intent(in) :: outdir
  character(len=256) command
  write(command, *) 'if [ ! -d ', trim(outdir), ' ]; then mkdir -p ', trim(outdir), '; fi'
  write(*, *) trim(command)
  call system(command)
end subroutine makedirs

!---------------------------
! in the following
!---------------------------


subroutine FindShockRadiusAnswer
  use unitsmod
  use fieldmod
  implicit none
  integer::i,j,k
  real(8):: pmax
  
  rshock = 0.0d0
  pmax = 0.0d0
  kTshock = 0.0d0
  Vshock  = 0.0d0
  k = ks
  j = js 
  do i=is,ie
     if(pmax < p(i,j,k)) then
        pmax = p(i,j,k)
        rshock = x1b(i)
        ! p = n T 
        kTshock = (kbol*1.0d5)*Tem(i,j,k)*erg_to_keV ! T [keV]
        Vshock = v1(i,j,k)/1.0e5 ! cm/s => km/s
     endif
  enddo
  !print *, "rshock=",rshock/pc,"[pc]"
  !print *, "kTshock=",kTshock,"[keV]"
  !print *, "Vshock=",Vshock,"[km/s]"
  
end subroutine FindShockRadiusAnswer

subroutine TimeProfleAnswer
  use unitsmod
  use fieldmod
  implicit none
  integer::i,j,k

  character(20),parameter::dirname="./"
  character(40)::filename
  integer,save::unittpr
  real(8)::Etot,pi

  logical,save:: is_inited
  data is_inited / .false. /

  if(.not. is_inited)then
     call makedirs(dirname)
  endif

  pi = acos(-1.0d0)
  Msw  = 0.0d0
  Etot = 0.0d0
  Lbol = 0.0d0
  k=ks
  j=js
  do i=is,ie
     if(x1b(i) <= rshock ) Msw  = Msw  + d(i,j,k)*dvl1a(i)*4.0d0*pi
     Lbol = Lbol + edot(i,j,k) * dvl1a(i)*4.0d0*pi
     Etot = Etot + (0.5d0*d(i,j,k)*v1(i,j,k)**2+ei(i,j,k))*dvl1a(i)*4.0d0*pi
  enddo
  !print *, "Msw=",Msw/Msolar,"[M_s]"
  write(filename,'(A)')"t-prof.dat"
  filename = trim(dirname)//filename
  if(.not. is_inited) print *,"time evolution is written in", filename
  if(.not. is_inited) open(newunit=unittpr,file=filename,status='replace',form='formatted')
  if(.not. is_inited) write(unittpr,'(1a,1x,A)') "#"," time[year] rshock[pc] Msw[Ms] Tshock[keV] Vshock[km/s] Lbol[erg/s] Etot[erg]"
  write(unittpr,'(1x,7(1x,E13.3))') time/year,rshock/pc,Msw/Msolar,kTshock,Vshock,Lbol,Etot
  
  !close(unittpr)
  
  is_inited = .true.

  return
end subroutine  TimeProfleAnswer
