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
    real(8),dimension(:),allocatable:: x1a,x2a,dvl2a
    real(8),dimension(:,:,:),allocatable:: d,v1,v2,v3,p,ei,gp
    real(8),dimension(:,:,:),allocatable:: Tem,edot
    real(8):: dx
    real(8):: gam,rho0,Eexp
    real(8),dimension(:,:),allocatable:: rshock_ray
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
     call FindShockRadius
     call Visualize1D
     call Visualize2D
     call TimeProfle
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
  read(unitinp,*) dummy,jzone,jgs
  close(unitinp)
  in=izone+2*igs
  jn=jzone+2*jgs
  kn=1
!  write(6,*)igs,jgs
  is=1+igs
  js=1+jgs
  ks=1
  ie=in-igs
  je=jn-jgs
  ke=1

  if(.not. is_inited)then
     allocate( x1b(in),x1a(in),dvl1a(in))
     allocate( x2b(jn),x2a(jn),dvl2a(jn))
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
  read(unitbin)x2b(:),x2a(:),dvl2a(:)
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
  do j=js,je
  do i =is,ie
     ! p = n k T => T = p/(n)/k since kbol[J/K] kbol*1.0d5 [erg/K] 
     Tem(i,j,k) = p(i,j,k)/(d(i,j,k)/mp) / (kbol*1.0d5) ! [K]
     edot(i,j,k) = 1.4d-27 * (d(i,j,k)/mp)**2 *sqrt(Tem(i,j,k)) ! erg/s/cm^3
  enddo
  enddo
end subroutine EstimateEmissivity

subroutine FindShockRadius
  use unitsmod
  use fieldmod
  implicit none
  integer::i,j,k
  real(8),dimension(:,:),allocatable,save:: pmax_ray,kTshock_ray,Vshock_ray
  logical::is_inited
  data is_inited /.false./

  if(.not. is_inited)then
     allocate(pmax_ray(jn,kn))
     allocate(rshock_ray ,mold=pmax_ray)
     allocate(kTshock_ray,mold=pmax_ray)
     allocate(Vshock_ray ,mold=pmax_ray)
  endif
  
  rshock_ray(:,:) = 0.0d0
  pmax_ray(:,:) = 0.0d0
  kTshock_ray(:,:) = 0.0d0
  Vshock_ray(:,:)  = 0.0d0
  k = ks
  do j=js,je 
     do i=is,ie
        if(pmax_ray(j,k) < p(i,j,k)) then
           pmax_ray(j,k) = p(i,j,k)
           rshock_ray = x1b(i)
           ! p = n T 
           kTshock_ray(j,k) = (kbol*1.0d5)*Tem(i,j,k)*erg_to_keV ! T [keV]
           Vshock_ray(j,k)  = v1(i,j,k)/1.0e5 ! cm/s => km/s
        endif
     enddo
  enddo
  rshock = 0.0d0
  do j=js,je
     if(rshock < rshock_ray(j,k) )then
        rshock  =  rshock_ray(j,k)
        kTshock = kTshock_ray(j,k)
        Vshock  =  Vshock_ray(j,k)
     endif
  enddo
  !print *, "rshock=",rshock/pc,"[pc]"
  !print *, "kTshock=",kTshock,"[keV]"
  !print *, "Vshock=",Vshock,"[km/s]"
  is_inited = .true.
end subroutine FindShockRadius

subroutine Visualize2D
  use unitsmod
  use fieldmod
  implicit none
  integer::i,j,k

  character(20),parameter::dirname="output/"
  character(40)::filename
  integer,parameter::unit2D=432

  real(8),dimension(:,:),allocatable,save::d2d,p2d,v12d,ed2d

  logical,save:: is_inited
  data is_inited / .false. /

  if(.not. is_inited)then
     call makedirs(dirname)
     is_inited = .true.
     allocate( d2d(in,jn))
     allocate( p2d,mold=d2d)
     allocate(v12d,mold=d2d)
     allocate(ed2d,mold=d2d)
  endif

  k = ks
! boundary 
  do i=is,ie
      d(i,js-1,k) =   d(i,js,k)
      p(i,js-1,k) =   p(i,js,k)
     v1(i,js-1,k) =  v1(i,js,k)
   edot(i,js-1,k) =edot(i,js,k)

      d(i,je+1,k) =   d(i,je,k)
      p(i,je+1,k) =   p(i,je,k)
     v1(i,je+1,k) =  v1(i,je,k)
   edot(i,je+1,k) =edot(i,je,k)
  enddo

  do j=js,je+1
  do i=is,ie
       d2d(i,j) =  0.5d0*(   d(i,j,k)+   d(i,j-1,k))
       p2d(i,j) =  0.5d0*(   p(i,j,k)+   p(i,j-1,k))
      v12d(i,j) =  0.5d0*(  v1(i,j,k)+  v1(i,j-1,k))
      ed2d(i,j) =  0.5d0*(edot(i,j,k)+edot(i,j-1,k))
  enddo
  enddo


  write(filename,'(a6,i5.5,a4)')"twopro",incr,".dat"
  filename = trim(dirname)//filename
  open(unit2D,file=filename,status='replace',form='formatted')

  write(unit2D,'(1a,a6,1(1x,E12.3))') "#"," time=",time
!                                    12345678    1234567890123   1234567890123   123456789012
  write(unit2D,'(1a,2(1x,a7,i0))') "#"," Nrad= ",ie-is+1," Nthe= ",je-js+2

  write(unit2D,'(1a,6(1x,a13))') "#","1:r[cm] ","2:theta[rad] ","3:den[1/cm^3] ","4:p[erg/cm3] ","5:vel[cm/s] ","6:edot[cgs]"

  do j=js,je+1
  do i=is,ie
     write(unit2D,'(1x,SP,6(1x,E13.3))') x1b(i),x2a(j),d2d(i,j)/mu,p2d(i,j),v12d(i,j),ed2d(i,j)
  enddo
     write(unit2D,*)
  enddo

  close(unit2D)

  return
end subroutine Visualize2D

subroutine Visualize1D
  use unitsmod
  use fieldmod
  implicit none
  integer::i,j,k

  character(20),parameter::dirname="output/"
  character(40)::filename
  integer,parameter::unit1D=123

  real(8),dimension(:),allocatable,save::d1d,p1d,v11d

  logical,save:: is_inited
  data is_inited / .false. /

  if(.not. is_inited)then
     call makedirs(dirname)
     is_inited = .true.
     allocate( d1d(in))
     allocate( p1d(in))
     allocate(v11d(in))
  endif


  d1d(:) = 0.0d0
  p1d(:) = 0.0d0
  v11d(:) = 0.0d0
  k=ks
  do i=is,ie
  do j=js,je
      d1d(i) =  d1d(i) +  d(i,j,k)*dvl2a(j)
      p1d(i) =  p1d(i) +  p(i,j,k)*dvl2a(j)
     v11d(i) = v11d(i) + v1(i,j,k)*dvl2a(j)
  enddo
      d1d(i) =  d1d(i)/sum(dvl2a(:))
      p1d(i) =  p1d(i)/sum(dvl2a(:))
     v11d(i) = v11d(i)/sum(dvl2a(:))
  enddo

  write(filename,'(a6,i5.5,a4)')"onepro",incr,".dat"
  filename = trim(dirname)//filename
  open(unit1D,file=filename,status='replace',form='formatted')

  write(unit1D,'(1a,a6,1(1x,E12.3))') "#"," time=",time
!                                     1234567890123   1234567890123   1234567890123   1234567890123
  write(unit1D,'(1a,4(1x,a13))') "#","1:r[cm]      ","2:den[1/cm^3]","3:p[erg/cm3] ","4:vel[cm/s]  "

  do i=is,ie
     write(unit1D,'(1x,SP,4(1x,E13.3))') x1b(i),d1d(i)/mu,p1d(i),v11d(i)
  enddo
  close(unit1D)

  return
end subroutine Visualize1D

subroutine TimeProfle
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

  Etot=0.0d0
  k=ks
  do j=js,je
  do i=is,ie
     if(x1b(i) <= rshock_ray(j,k) ) Msw  = Msw  + d(i,j,k)*dvl1a(i)*4.0d0*pi
     Lbol = Lbol + edot(i,j,k) * dvl1a(i)*4.0d0*pi
     Etot = Etot + (0.5d0*d(i,j,k)*(v1(i,j,k)**2+v2(i,j,k)**2)+ei(i,j,k))*dvl1a(i)*dvl2a(j)*2.0d0*pi
  enddo
  enddo

  if(.not. is_inited) then
     write(filename,'(A)')"t-prof.dat"
     filename = trim(dirname)//filename
     print *,"time evolution is written in", filename
     open(newunit=unittpr,file=filename,status='replace',form='formatted')
     write(unittpr,'(1a,1x,A)') "#"," time[year] rshock[pc] Msw[Ms] Tshock[keV] Vshock[km/s] Lbol[erg/s] Etot[erg]"
  endif
!  write(unittot,'(1a,4(1x,E12.3))') "#",time/year
!                                    12345678   1234567890123     1234567890123   123456789012
!  write(unittot,'(1a,4(1x,a13))') "#","1:r[pc] ","2:den[1/cm^3] ","3:p[erg/cm3] ","4:vel[km/s] "

  write(unittpr,'(1x,7(1x,E13.3))') time/year,rshock/pc,Msw/Msolar,kTshock,Vshock,Lbol,Etot
  ! close(unittot)

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


