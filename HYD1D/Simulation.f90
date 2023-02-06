
      module unitsmod
      implicit none
      real(8),parameter::    pc  = 3.085677581d18   ! parsec in [cm]
      real(8),parameter::    mu  = 1.660539066d-24  ! g
      real(8),parameter:: Msolar = 1.989e33         ! g
      real(8),parameter::   kbol = 1.380649d-23     ! J/K
      real(8),parameter::   year = 365.0d0*24*60*60 ! sec

      end module unitsmod

      module commons
      use unitsmod
      implicit none
      integer::nhy
      integer,parameter::nhymax=100000
      real(8)::time,dt
      real(8),parameter:: Coul=0.25d0
      data time / 0.0d0 /
      real(8),parameter:: timemax=1.0d5*year
      real(8),parameter:: dtout=timemax/500

      integer,parameter::izones=200
      integer,parameter::jzones=1
      integer,parameter::mgn=2
      integer,parameter::in=izones+2*mgn+1 &
     &                  ,jn=1 &
     &                  ,kn=1
      integer,parameter::is=mgn+1 &
     &                  ,js=1 &
     &                  ,ks=1
      integer,parameter::ie=izones+mgn &
     &                  ,je=1 &
     &                  ,ke=1


      real(8),parameter:: x1min=0.0d0,x1max=1.0d2*pc
      real(8),dimension(in)::x1a,x1b,dvl1a
      real(8),dimension(jn)::x2a,x2b
      real(8),dimension(kn)::x3a,x3b

      real(8),dimension(in,jn,kn)::d,et,mv1,mv2,mv3
      real(8),dimension(in,jn,kn)::p,ei,v1,v2,v3,cs
      real(8),dimension(in,jn,kn)::gp,gp1a,gp2a


      end module commons
     
      module eosmod
      implicit none
! adiabatic
      real(8),parameter::gam=4.0d0/3.0d0 !! adiabatic index
! isothermal
!      real(8)::csiso  !! isothemal sound speed
end module eosmod

      module fluxmod
      use commons, only : in,jn,kn
      implicit none
      integer,parameter::nden=1,nve1=2,nve2=3,nve3=4,nene=5,npre=6,ncsp=7
      integer,parameter::nhyd=7
      real(8),dimension(nhyd,in,jn,kn):: svc

      integer,parameter::mudn=1,muvu=2,muvv=3,muvw=4,muet=5  &
     &                  ,mfdn=6,mfvu=7,mfvv=8,mfvw=9,mfet=10 &
     &                  ,mcsp=11,mvel=12,mpre=13
      integer,parameter:: mflx=5,madd=3

      integer,parameter:: mden=1,mrv1=2,mrv2=3,mrv3=4,meto=5 &
     &                          ,mrvu=muvu,mrvv=muvv,mrvw=muvw
      real(8),dimension(mflx,in,jn,kn):: nflux1,nflux2,nflux3
      real(8),dimension(in,jn,kn):: grvsrc1,grvsrc2,grvsrc3
      real(8),dimension(in,jn,kn):: adsrca

      end module fluxmod

      program main
      use commons
      implicit none
      integer,parameter:: nhyspan=1000
      write(6,*) "setup grids and fields"
      call GenerateGrid
      call GenerateProblem
      call ConsvVariable
      write(6,*) "entering main loop"
! main loop
                                  write(6,*)"step","time [yr]","dt [yr]"
      mloop: do nhy=1,nhymax
         call TimestepControl
!         if(mod(nhy,nhyspan) .eq. 0 ) write(6,*)nhy,time/year,dt/year
         write(6,*)nhy,time/year,dt/year
         call BoundaryCondition
         call StateVector
         call NumericalFlux1
         call UpdateConsv
         call PrimVariable
         time=time+dt
         call Output
         if(time > timemax) exit mloop
      enddo mloop

      write(6,*) "program has been finished"
      end program main

      subroutine GenerateGrid
      use commons
      implicit none
      real(8)::dx,dy
      integer::i,j,k
      dx=(x1max-x1min)/izones
      do i=1,in
         x1a(i) = dx*(i-(mgn+1))+x1min
      enddo
      do i=1,in-1
         x1b(i) = 0.5d0*(x1a(i+1)+x1a(i))
      enddo
      do i=is,ie
         dvl1a(i) = (x1a(i+1)**3-x1a(i)**3)/3.0d0
      enddo

      return
      end subroutine GenerateGrid

      subroutine GenerateProblem
      use commons
      use eosmod
      implicit none
      integer::i,j,k
      real(8),parameter:: neu = 3.0d0
      real(8):: rho1,rho2
      real(8):: ein0
      real(8):: pre1,pre2
      real(8):: vel1,vel2
      real(8):: dr
      real(8):: pi
      real(8):: frac,eexp

      pi =acos(-1.0d0)
      dr = 8.0d0*(x1a(is+1)-x1a(is)) ! 8 mesh
      write(6,*) "shell length [pc]",dr/pc

! circum steller  medium
      rho2 = 1.0d0*mu ! Intersteller medium 1 [1/cm^3]
      pre2 = rho2* kbol *1.0d4 ! 10^4 [K]
      vel2 = 0.0d0

! blast wave
      frac = 0.8d0
      rho1 = (10.0d0*Msolar)/(4.0*pi/3.0d0*dr**3)
      eexp = frac*(1.0d51)
      pre1 = eexp/(4.0*pi/3.0d0*dr**3)*(gam-1.0d0)  
      vel1 = sqrt((1.0d0-frac)*eexp/(4.0*pi/3.0d0*dr**3)/rho1)

      write(6,*) "Eex= ",frac   ,"[10^51 erg]"
      write(6,*) "rho= ",rho1/mu,"[1/cm^3]"
      write(6,*) "vel= ",vel1   ,"[cm/s]"
      write(6,*) "pre= ",pre1   ,"[erg/cm^3]"
     
      d(:,:,:) = rho2
  
      do k=ks,ke
      do j=js,je
      do i=is,ie
         if(x1b(i) < dr)then
             d(i,j,k) = max(rho1*(x1b(i)/dr)**(neu/(gam-1)),rho2)
             p(i,j,k) = pre1
            v1(i,j,k) = vel1*max(x1b(i)/dr,0.0d0)
         else
             d(i,j,k) = rho2
             p(i,j,k) = pre2
            v1(i,j,k) = vel2
         endif
      enddo
      enddo
      enddo

      do k=ks,ke
      do j=js,je
      do i=is,ie
          ei(i,j,k) = p(i,j,k)/(gam-1.0d0)
          cs(i,j,k) = sqrt(gam*p(i,j,k)/d(i,j,k))
      enddo
      enddo
      enddo
      

      call BoundaryCondition

      return
      end subroutine GenerateProblem


      subroutine BoundaryCondition
      use commons
      implicit none
      integer::i,j,k

!reflection
      do k=ks,ke
      do j=js,je
      do i=1,mgn
           d(is-i,j,k) =   d(is+i-1,j,k)
          ei(is-i,j,k) =  ei(is+i-1,j,k)
          v1(is-i,j,k) =  v1(is+i-1,j,k)
          v2(is-i,j,k) =  v2(is+i-1,j,k)
          v3(is-i,j,k) =  v3(is+i-1,j,k)
          gp(is-i,j,k) =  gp(is+i-1,j,k)
      enddo
      enddo
      enddo

      do k=ks,ke
      do j=js,je
      do i=1,mgn
           d(ie+i,j,k) =   d(ie-i+1,j,k)
          ei(ie+i,j,k) =  ei(ie-i+1,j,k)
          v1(ie+i,j,k) =  v1(ie-i+1,j,k)
          v2(ie+i,j,k) =  v2(ie-i+1,j,k)
          v3(ie+i,j,k) =  v3(ie-i+1,j,k)
          gp(ie+i,j,k) =  gp(ie-i+1,j,k)
      enddo
      enddo
      enddo

      return
      end subroutine BoundaryCondition

      subroutine ConsvVariable
      use commons
      implicit none
      integer::i,j,k
      do k=ks,ke
      do j=js,je
      do i=is,ie
          et(i,j,k) = 0.5d0*d(i,j,k)*(   &
     &                    +v1(i,j,k)**2  &
     &                    +v2(i,j,k)**2  &
     &                    +v3(i,j,k)**2) &
     &                    +ei(i,j,k)
          mv1(i,j,k) =d(i,j,k)*v1(i,j,k)
          mv2(i,j,k) =d(i,j,k)*v2(i,j,k)
          mv3(i,j,k) =d(i,j,k)*v3(i,j,k)
      enddo
      enddo
      enddo
      
      return
      end subroutine Consvvariable

      subroutine PrimVariable
      use commons
      use eosmod
      implicit none
      integer::i,j,k
      do k=ks,ke
      do j=js,je
      do i=is,ie
          v1(i,j,k) = mv1(i,j,k)/d(i,j,k)
          v2(i,j,k) = mv2(i,j,k)/d(i,j,k)
          v3(i,j,k) = mv3(i,j,k)/d(i,j,k)

          ei(i,j,k) =  et(i,j,k)          &
     &          -0.5d0*d(i,j,k)*(         &
     &                    +v1(i,j,k)**2   &
     &                    +v2(i,j,k)**2   &
     &                    +v3(i,j,k)**2)

! adiabatic
           p(i,j,k) =  ei(i,j,k)*(gam-1.0d0)
          cs(i,j,k) =  sqrt(gam*p(i,j,k)/d(i,j,k))
! isotermal
!           p(i,j,k) =  d(i,j,k)*csiso**2
!          cs(i,j,k) =  csiso
      enddo
      enddo
      enddo

      return
      end subroutine PrimVariable

      subroutine TimestepControl
      use commons
      implicit none
      real(8)::dtl1
      real(8)::dtl2
      real(8)::dtl3
      real(8)::dtlocal
      real(8)::dtmin
      integer::i,j,k
      dtmin=1.0d90
      do k=ks,ke
      do j=js,je
      do i=is,ie
         dtl1 =(x1a(i+1)-x1a(i))/(abs(v1(i,j,k)) +cs(i,j,k))
!         dtl2 =(x2a(j+1)-x2a(j))/(abs(v2(i,j,k)) +cs(i,j,k))
!         dtl3 =(x1a(i+1)-x1a(i))/(abs(v1(i,j,k)) +cs(i,j,k))
!         dtlocal = min (dtl1,dtl2)
         dtlocal = dtl1
         if(dtlocal .lt. dtmin) dtmin = dtlocal
      enddo
      enddo
      enddo

      dt = Coul * dtmin
!      write(6,*)"dt",dt
      return
      end subroutine TimestepControl

      subroutine StateVector
      use commons
      use fluxmod
      use eosmod
      implicit none
      integer::i,j,k

!      do j=1,jn-1

      do k=ks,ke
      do j=js,je
      do i=1,in-1
         svc(nden,i,j,k) =  d(i,j,k)
         svc(nve1,i,j,k) = v1(i,j,k)
         svc(nve2,i,j,k) = v2(i,j,k)
         svc(nve3,i,j,k) = v3(i,j,k)
! adiabatic
         svc(nene,i,j,k) = ei(i,j,k)/d(i,j,k)
         svc(npre,i,j,k) = ei(i,j,k)*(gam-1.0d0)
         svc(ncsp,i,j,k) = sqrt(gam*(gam-1.0d0)*ei(i,j,k)/d(i,j,k))
! isotermal
!         svc(nene,i,j,k) = csiso**2
!         svc(npre,i,j,k) = d(i,j,k)*csiso**2
!         svc(ncsp,i,j,k) = csiso

 ! for output
         p(i,j,k) = svc(npre,i,j,k) 

!      integer,parameter:: mden=1,mrv1=2,mrv2=3,mrv3=4,meto=5 &
! (P_tt+P_pp) 
         adsrca(i,j,k) = &
     &   2.0d0* p(i,j,k)     &
     & + d(i,j,k)*(v2(i,j,k)**2+v3(i,j,k)**2)
 
      enddo
      enddo
      enddo



      return
      end subroutine StateVector

      subroutine minmod(a,b,d)
      use fluxmod, only : nhyd
      implicit none
      real(8),dimension(nhyd),intent(in)::a,b
      real(8),dimension(nhyd),intent(out)::d
      integer:: n

      do n=1,nhyd
         d(n) = sign(1.0d0,a(n))*max(0.0d0,min(abs(a(n))                &
     &                                        ,sign(1.0d0,a(n))*b(n)))
      enddo

      return
      end subroutine minmod


      subroutine vanLeer(dvp,dvm,cf,cb,dv)
      use fluxmod, only : nhyd
      implicit none
      real(8),dimension(nhyd),intent(in)::dvp,dvm
      real(8),intent(in)::cf,cb
      real(8),dimension(nhyd),intent(out)::dv
      integer:: n

      do n=1,nhyd
         if(dvp(n)*dvm(n) .gt. 0.0d0)then
            dv(n) =  dvp(n)*dvm(n) &
     &              * (cf*dvm(n)+cb*dvp(n))&
     &              / ( dvm(n)**2 +(cf+cb-2.0d0)*dvp(n)*dvm(n) +dvp(n)**2)
         else
            dv(n) = 0.0d0
         endif
      enddo

      return
      end subroutine vanLeer

      subroutine MClimiter(a,b,c,d)
      use fluxmod, only : nhyd
      implicit none
      real(8),dimension(nhyd),intent(in)::a,b,c
      real(8),dimension(nhyd),intent(out)::d
      integer:: n

      do n=1,nhyd
         d(n) = sign(1.0d0,a(n))*max(0.0d0,min(abs(a(n))         &
     &                                  ,sign(1.0d0,a(n))*b(n)   &
     &                                  ,sign(1.0d0,a(n))*c(n))) 
      enddo

      return
      end subroutine MClimiter

      subroutine NumericalFlux1
      use commons, only: is,ie,in,js,je,jn,ks,ke,kn,mgn,x1a,x1b
      use fluxmod
      implicit none
      integer::i,j,k
      real(8),dimension(nhyd):: dsvp,dsvm,dsvc,dsv
      real(8),dimension(nhyd,in,jn,kn):: leftpr,rigtpr
      real(8),dimension(2*mflx+madd,in,jn,kn):: leftco,rigtco
      real(8),dimension(2*mflx+madd):: leftst,rigtst
      real(8),dimension(mflx):: nflux
      real(8),dimension(in),save:: x1c,ctl,ctr
      real(8),dimension(in),save:: bck,frd,cf,cb
      real(8):: xii,cflo,cblo
!
! Mignone 2014 Ref. [1]
!
      logical,save:: is_inited
      data is_inited / .false. /
      do i=is-mgn,ie+mgn
         xii = (x1a(i+1)-x1a(i))/x1b(i  )
! Eq. (C10) of Ref. [1]
         x1c(i) = x1b(i  )*(1.0d0 + xii**2 * 2.0d0/(12.0d0 + xii**2))
! Eq. (C2), (C3) of Ref. [1]
         ctl(i) = 1.0d0+xii*(4.0d0/(12.0d0 + xii**2))
         ctr(i) = 1.0d0-xii*(4.0d0/(12.0d0 + xii**2))
      enddo

      do i=is-mgn+1,ie+mgn-1
! Eq. (C8), (C9) of Ref. [1]
         bck(i) = (x1a(i+1)-x1a(i))/(x1c(i  )- x1c(i-1))
         frd(i) = (x1a(i+1)-x1a(i))/(x1c(i+1)- x1c(i  ))

! Eq. (C12), (C13) of Ref. [1]
         cf(i) = (x1c(i+1)- x1c(i  ))/( x1a(i+1) -x1c(i  ))
         cb(i) = (x1c(i  )- x1c(i-1))/( x1c(i  ) -x1a(i  ))
     enddo

      k=ks
      do j=js,je
      do i=is-1,ie+1
         dsvp(:) = (svc(:,i+1,j,k) -svc(:,i,j,k)                 )
         dsvm(:) = (                svc(:,i,j,k) - svc(:,i-1,j,k))

         dsvp(:) = dsvp(:) * frd(i)
         dsvm(:) = dsvm(:) * bck(i)
         cflo = cf(i)
         cblo = cb(i)
         call vanLeer(dsvp,dsvm,cflo,cblo,dsv)
!         call minmod(dsvp,dsvm,dsv)
         leftpr(:,i+1,j,k) = svc(:,i,j,k) + 0.5d0*dsv(:)*ctr(i)
         rigtpr(:,i  ,j,k) = svc(:,i,j,k) - 0.5d0*dsv(:)*ctl(i)
      enddo
      enddo

      do j=js,je
      do i=is,ie+1
         leftco(mudn,i,j,k)=leftpr(nden,i,j,k) ! rho
         leftco(muvu,i,j,k)=leftpr(nve1,i,j,k)*leftpr(nden,i,j,k)  ! rho v_x
         leftco(muvv,i,j,k)=leftpr(nve2,i,j,k)*leftpr(nden,i,j,k)  ! rho v_y
         leftco(muvw,i,j,k)=leftpr(nve3,i,j,k)*leftpr(nden,i,j,k)  ! rho v_z
         leftco(muet,i,j,k)=leftpr(nene,i,j,k)*leftpr(nden,i,j,k) &! e_i+ rho v^2/2
     &               +0.5d0*leftpr(nden,i,j,k)*(    &
     &                     +leftpr(nve1,i,j,k)**2   &
     &                     +leftpr(nve2,i,j,k)**2   &
     &                     +leftpr(nve3,i,j,k)**2)

         leftco(mfdn,i,j,k)=leftpr(nden,i,j,k)                   *leftpr(nve1,i,j,k)
         leftco(mfvu,i,j,k)=leftpr(nden,i,j,k)*leftpr(nve1,i,j,k)*leftpr(nve1,i,j,k) &
     &                     +leftpr(npre,i,j,k)
         leftco(mfvv,i,j,k)=leftpr(nden,i,j,k)*leftpr(nve2,i,j,k)*leftpr(nve1,i,j,k)
         leftco(mfvw,i,j,k)=leftpr(nden,i,j,k)*leftpr(nve3,i,j,k)*leftpr(nve1,i,j,k)
         leftco(mfet,i,j,k)=(leftpr(nene,i,j,k)*leftpr(nden,i,j,k)  &
     &               +0.5d0*leftpr(nden,i,j,k)*(   &
     &                     +leftpr(nve1,i,j,k)**2  &
     &                     +leftpr(nve2,i,j,k)**2  &
     &                     +leftpr(nve3,i,j,k)**2) &
     &                     +leftpr(npre,i,j,k)     &
     &                       )                                  *leftpr(nve1,i,j,k) 

         leftco(mcsp,i,j,k)= leftpr(ncsp,i,j,k)
         leftco(mvel,i,j,k)= leftpr(nve1,i,j,k)
         leftco(mpre,i,j,k)= leftpr(npre,i,j,k)


         rigtco(mudn,i,j,k)=rigtpr(nden,i,j,k)
         rigtco(muvu,i,j,k)=rigtpr(nve1,i,j,k)*rigtpr(nden,i,j,k)
         rigtco(muvv,i,j,k)=rigtpr(nve2,i,j,k)*rigtpr(nden,i,j,k)
         rigtco(muvw,i,j,k)=rigtpr(nve3,i,j,k)*rigtpr(nden,i,j,k)
         rigtco(muet,i,j,k)=rigtpr(nene,i,j,k)*rigtpr(nden,i,j,k) &
     &               +0.5d0*rigtpr(nden,i,j,k)*(  &
     &                     +rigtpr(nve1,i,j,k)**2 &
     &                     +rigtpr(nve2,i,j,k)**2 &
     &                     +rigtpr(nve3,i,j,k)**2)

         rigtco(mfdn,i,j,k)=rigtpr(nden,i,j,k)                   *rigtpr(nve1,i,j,k)
         rigtco(mfvu,i,j,k)=rigtpr(nden,i,j,k)*rigtpr(nve1,i,j,k)*rigtpr(nve1,i,j,k) &
     &                     +rigtpr(npre,i,j,k)
         rigtco(mfvv,i,j,k)=rigtpr(nden,i,j,k)*rigtpr(nve2,i,j,k)*rigtpr(nve1,i,j,k)
         rigtco(mfvw,i,j,k)=rigtpr(nden,i,j,k)*rigtpr(nve3,i,j,k)*rigtpr(nve1,i,j,k)
         rigtco(mfet,i,j,k)=(rigtpr(nene,i,j,k)*rigtpr(nden,i,j,k) &
     &               +0.5d0*rigtpr(nden,i,j,k)*(   &
     &                     +rigtpr(nve1,i,j,k)**2  &
     &                     +rigtpr(nve2,i,j,k)**2  &
     &                     +rigtpr(nve3,i,j,k)**2) &
     &                     +rigtpr(npre,i,j,k)     &
     &                      )                                    *rigtpr(nve1,i,j,k)

         rigtco(mcsp,i,j,k)= rigtpr(ncsp,i,j,k)
         rigtco(mvel,i,j,k)= rigtpr(nve1,i,j,k)
         rigtco(mpre,i,j,k)= rigtpr(npre,i,j,k)

      enddo
      enddo

      do j=js,je
      do i=is,ie+1
         leftst(:)=leftco(:,i,j,k)
         rigtst(:)=rigtco(:,i,j,k)
!         call HLLE(leftst,rigtst,nflux)
         call HLLC(leftst,rigtst,nflux)
         nflux1(mden,i,j,k)=nflux(mden)
         nflux1(mrv1,i,j,k)=nflux(mrvu)
         nflux1(mrv2,i,j,k)=nflux(mrvv)
         nflux1(mrv3,i,j,k)=nflux(mrvw)
         nflux1(meto,i,j,k)=nflux(meto)
      enddo
      enddo

      return
      end subroutine Numericalflux1

      subroutine HLLE(leftst,rigtst,nflux)
      use fluxmod
      implicit none
      real(8),dimension(2*mflx+madd),intent(in)::leftst,rigtst
      real(8),dimension(mflx),intent(out)::nflux
      real(8),dimension(mflx)::ul,ur,fl,fr
      real(8)::csl,csr
      real(8):: vl, vr
      real(8):: sl, sr

      ul(1:mflx) = leftst(1:mflx)
      fl(1:mflx) = leftst(mflx+1:2*mflx)
      ur(1:mflx) = rigtst(1:mflx)
      fr(1:mflx) = rigtst(mflx+1:2*mflx)
      csl=leftst(mcsp)
      csr=rigtst(mcsp)
       vl=leftst(mvel)
       vr=rigtst(mvel)

       sl = min(vl,vr) - max(csl,csr)
       sl = min(0.0d0,sl)
       sr = max(vl,vr) + max(csl,csr)
       sr = max(0.0d0,sr)

       nflux(:) = (sr*fl(:)-sl*fr(:) +sl*sr*(ur(:)-ul(:)))/(sr-sl)

      return
      end subroutine HLLE

      subroutine HLLC(leftst,rigtst,nflux)
!=====================================================================
!
! HLLC Scheme
!
! Purpose
! Calculation of Numerical Flux by HLLC method
!
! Reference
!  Toro EF, Spruce M, Speares W. (1992,1994)
!
! Input
! Output
!=====================================================================
      use fluxmod, only: mflx,madd                 &
     &                 , mudn,muvu,muvv,muvw,muet  &
     &                 , mfdn,mfvu,mfvv,mfvw,mfet  &
     &                 , mcsp,mvel,mpre            &
     &                 , mden,mrvu,mrvv,mrvw,meto

      implicit none
      real(8),dimension(2*mflx+madd),intent(in)::leftst,rigtst
      real(8),dimension(mflx),intent(out)::nflux

!----- U -----
! qql :: left state
! qqr :: right state
      real(8) :: rol,vxl,vyl,vzl,ptl,eel
      real(8) :: ror,vxr,vyr,vzr,ptr,eer
      real(8) :: rxl,ryl,rzl
      real(8) :: rxr,ryr,rzr
      real(8) :: ptst

!----- U* ----
! qqlst ::  left state
! qqrst :: right state
      real(8) :: rolst,vxlst,vylst,vzlst,eelst
      real(8) :: rorst,vxrst,vyrst,vzrst,eerst
      real(8) :: rxlst,rylst,rzlst
      real(8) :: rxrst,ryrst,rzrst

!----- flux ---
! fqql ::  left physical flux
! fqqr :: right physical flux
      real(8) :: frol,frxl,fryl,frzl,feel
      real(8) :: fror,frxr,fryr,frzr,feer

!----- wave speed ---
! sl ::  left-going fastest signal velocity
! sr :: right-going fastest signal velocity
! sm :: contact discontinuity velocity
! slst ::  left-going alfven velocity
! srst :: right-going alfven velocity
      real(8) :: sm,sl,sr

! cfl :: left-state Fast wave velocity
! cfr :: right-sate Fast wave velocity
      real(8) :: cfl,cfr

!--------------------
! temporary variables
      real(8) :: sdl,sdr,sdml,sdmr,isdml,isdmr,rosdl,rosdr
      real(8) :: temp
  
! no if
      real(8) :: sign1,maxs1,mins1
      real(8) :: msl,msr

!----- Step 0. ----------------------------------------------------------|

!---- Left state
        
        rol = leftst(mudn)
        eel = leftst(muet)
        rxl = leftst(muvu)
        ryl = leftst(muvv)
        rzl = leftst(muvw)
        vxl = leftst(muvu)/leftst(mudn)
        vyl = leftst(muvv)/leftst(mudn)
        vzl = leftst(muvw)/leftst(mudn)
        ptl = leftst(mpre)

!---- Right state
        
        ror = rigtst(mudn)
        eer = rigtst(muet)
        rxr = rigtst(muvu)
        ryr = rigtst(muvv)
        rzr = rigtst(muvw)
        vxr = rigtst(muvu)/rigtst(mudn)
        vyr = rigtst(muvv)/rigtst(mudn)
        vzr = rigtst(muvw)/rigtst(mudn)
        ptr = rigtst(mpre)
!----- Step 1. ----------------------------------------------------------|
! Compute wave left & right wave speed
!
         
        cfl = leftst(mcsp)
        cfr = rigtst(mcsp)

        sl = min(vxl,vxr)-max(cfl,cfr) ! note sl is negative
        sr = max(vxl,vxr)+max(cfl,cfr)
!----- Step 2. ----------------------------------------------------------|
! compute L/R fluxs
!
! Left value
        frol = leftst(mfdn)
        feel = leftst(mfet)
        frxl = leftst(mfvu)
        fryl = leftst(mfvv)
        frzl = leftst(mfvw)

! Right value
! Left value
        fror = rigtst(mfdn)
        feer = rigtst(mfet)
        frxr = rigtst(mfvu)
        fryr = rigtst(mfvv) 
        frzr = rigtst(mfvw)

!----- Step 4. ----------------------------------------------------------|
! compute middle and alfven wave
!
        sdl = sl - vxl
        sdr = sr - vxr
        rosdl = rol*sdl
        rosdr = ror*sdr

        temp = 1.0d0/(rosdr - rosdl)
! Eq. 45
        sm = (rosdr*vxr - rosdl*vxl - ptr + ptl)*temp
           
        sdml = sl - sm; isdml = 1.0d0/sdml
        sdmr = sr - sm; isdmr = 1.0d0/sdmr
        
!----- Step 5. ----------------------------------------------------------|
! compute intermediate states
!
! Eq. 49
        ptst = (rosdr*ptl-rosdl*ptr+rosdl*rosdr*(vxr-vxl))*temp

!----- Step 5A. ----------------------------------------------------------|
! compute Ul*
!

        rolst = rol*sdl   *isdml
        vxlst = sm
        rxlst = rolst*vxlst
           
        vylst = vyl
        rylst = rolst*vylst
        vzlst = vzl
        rzlst = rolst*vzlst

        eelst =(sdl*eel - ptl*vxl + ptst*sm  )*isdml

!----- Step 5B. ----------------------------------------------------------|
! compute Ur*
!

        rorst   = rosdr   *isdmr
        vxrst = sm
        rxrst = rorst*vxrst
        vyrst = vyr
        ryrst = rorst*vyrst
        vzrst = vzr
        rzrst = rorst*vzrst
           
        eerst = (sdr*eer - ptr*vxr  + ptst*sm  )*isdmr
              
!----- Step 6. ----------------------------------------------------------|
! compute flux
        sign1 = sign(1.0d0,sm)    ! 1 for sm>0, -1 for sm<0
        maxs1 =  max(0.0d0,sign1) ! 1 sm>0, 0 for sm<0
        mins1 = -min(0.0d0,sign1) ! 0 sm>0,-1 for sm<0

        msl   = min(sl  ,0.0d0)   ! 0 for sl > 0, sl for sl < 0
        msr   = max(sr  ,0.0d0)   ! S_R > 0

        nflux(mden) = (frol+msl*(rolst-rol))*maxs1 &
     &               +(fror+msr*(rorst-ror))*mins1
        nflux(meto) = (feel+msl*(eelst-eel))*maxs1 &
     &               +(feer+msr*(eerst-eer))*mins1
        nflux(mrvu) = (frxl+msl*(rxlst-rxl))*maxs1 &
     &               +(frxr+msr*(rxrst-rxr))*mins1
        nflux(mrvv) = (fryl+msl*(rylst-ryl))*maxs1 &
     &               +(fryr+msr*(ryrst-ryr))*mins1
        nflux(mrvw) = (frzl+msl*(rzlst-rzl))*maxs1 &
     &               +(frzr+msr*(rzrst-rzr))*mins1

      return
      end subroutine HLLC

      subroutine GravForce
      use commons
      use fluxmod
      implicit none
      integer :: i,j,k,n

      do k=ks,ke
      do j=js,je
      do i=is,ie+1
         gp1a(i  ,j,k) = gp(i,j,k) &
     & - 0.5d0*(gp(i  ,j,k)-gp(i-1,j,k))

         gp1a(i+1,j,k) = gp(i,j,k) &
     & + 0.5d0*(gp(i+1,j,k)-gp(i  ,j,k))

       grvsrc1(i,j,k) = (gp1a(i+1,j,k)-gp1a(i,j,k))/(x1a(i+1)-x1a(i))*d(i,j,k)

      enddo
      enddo
      enddo

      do k=ks,ke
      do i=is,ie
      do j=js,je+1
         gp2a(i  ,j,k) = gp(i,j,k) &
     & - 0.5d0*(gp(i  ,j,k)-gp(i,j-1,k))

         gp2a(i,j+1,k) = gp(i,j,k) &
     & + 0.5d0*(gp(i,j+1,k)-gp(i  ,j,k))

       grvsrc2(i,j,k) = (gp2a(i,j+1,k)-gp2a(i,j,k))/(x2a(j+1)-x2a(j))*d(i,j,k)

      enddo
      enddo
      enddo

       grvsrc3(:,:,:) = 0.0d0

      return
      end subroutine  GravForce

      subroutine UpdateConsv
      use commons
      use fluxmod
      implicit none
      integer::i,j,k
      real(8),dimension(in),save::as1,dv1i
      logical,save:: is_inited
      data is_inited / .false. /

      if (.not. is_inited ) then
         do i=is-1,ie+1 
            as1(i)  =  x1a(i)**2
         enddo
         do i=is,ie
            dv1i(i) = 1.0d0/dvl1a(i) 
         enddo
         
         is_inited = .true.
      endif


      do k=ks,ke
      do j=js,je
      do i=is,ie
         
         d(i,j,k) = d(i,j,k)                       &
     & +dt*(                                       &
     & +(- nflux1(mden,i+1,j,k)*as1(i+1)           &
     &   + nflux1(mden,i  ,j,k)*as1(i  ))*dv1i(i)  &
     &      )

         mv1(i,j,k) = mv1(i,j,k)                   &
     & +dt*(                                       &
     &      +  grvsrc1(i,j,k)                      &
     & +(- (nflux1(mrv1,i+1,j,k)-0.5d0*adsrca(i,j,k))*as1(i+1)           &
     &   + (nflux1(mrv1,i  ,j,k)-0.5d0*adsrca(i,j,k))*as1(i  ))*dv1i(i)  &
     &      )

         mv2(i,j,k) = mv2(i,j,k)                   &    
     & +dt*(                                       &
     &      +  grvsrc2(i,j,k)                      &
     & +(- nflux1(mrv2,i+1,j,k)*as1(i+1)           &
     &   + nflux1(mrv2,i  ,j,k)*as1(i  ))*dv1i(i)  &
     &      )

         mv3(i,j,k) = mv3(i,j,k)                   & 
     & +dt*(                                       &
     &      +  grvsrc3(i,j,k)                      &
     & +(- nflux1(mrv3,i+1,j,k)*as1(i+1)           &
     &   + nflux1(mrv3,i  ,j,k)*as1(i  ))*dv1i(i)  &
     &      )

          et(i,j,k) = et(i,j,k)                    &
     & +dt*(                                       &
     &  (- nflux1(meto,i+1,j,k)*as1(i+1)           &
     &   + nflux1(meto,i  ,j,k)*as1(i  ))*dv1i(i)  &
     &      )
      enddo
      enddo
      enddo

      return
      end subroutine UpdateConsv

      subroutine Output
      use commons
      implicit none
      integer::i,j,k
      character(20),parameter::dirname="bindata/"
      character(40)::filename
      real(8),save::tout
      data tout / 0.0d0 /
      integer::nout
      data nout / 1 /
      integer,parameter:: unitout=17
      integer,parameter:: unitbin=13
      integer,parameter:: gs=0
      integer,parameter:: nvar=6
      real(8)::x1out(is-gs:ie+gs,3)
!      real(8)::x2out(js-gs:je+gs,2)
      real(8)::hydout(is-gs:ie+gs,js-gs:je+gs,ks,nvar)

      logical, save:: is_inited
      data is_inited /.false./

      if (.not. is_inited) then
         call makedirs("bindata")
         is_inited =.true.
      endif

      if(time .lt. tout+dtout) return

      write(filename,'(a3,i5.5,a4)')"unf",nout,".dat"
      filename = trim(dirname)//filename

      open(unitout,file=filename,status='replace',form='formatted')
      write(unitout,*) "# ",time,dt
      write(unitout,*) "# ",izones,gs
      close(unitout)

      x1out(is-gs:ie+gs,1) = x1b(is-gs:ie+gs)
      x1out(is-gs:ie+gs,2) = x1a(is-gs:ie+gs)
      x1out(is-gs:ie+gs,3) = dvl1a(is-gs:ie+gs)

!      x2out(js-gs:je+gs,1) = x2b(js-gs:je+gs)
!      x2out(js-gs:je+gs,2) = x2a(js-gs:je+gs)

      hydout(is-gs:ie+gs,js-gs:je+gs,ks,1) =  d(is-gs:ie+gs,js-gs:je+gs,ks)
      hydout(is-gs:ie+gs,js-gs:je+gs,ks,2) = v1(is-gs:ie+gs,js-gs:je+gs,ks)
      hydout(is-gs:ie+gs,js-gs:je+gs,ks,3) = v2(is-gs:ie+gs,js-gs:je+gs,ks)
      hydout(is-gs:ie+gs,js-gs:je+gs,ks,4) = v3(is-gs:ie+gs,js-gs:je+gs,ks)
      hydout(is-gs:ie+gs,js-gs:je+gs,ks,5) =  p(is-gs:ie+gs,js-gs:je+gs,ks)
      hydout(is-gs:ie+gs,js-gs:je+gs,ks,6) = ei(is-gs:ie+gs,js-gs:je+gs,ks)
!      hydout(is-gs:ie+gs,js-gs:je+gs,ks,6) = gp(is-gs:ie+gs,js-gs:je+gs,ks)

      write(filename,'(a3,i5.5,a4)')"bin",nout,".dat"
      filename = trim(dirname)//filename
      open(unitbin,file=filename,status='replace',form='binary') 
      write(unitbin) x1out(:,:)
!      write(unitbin) x2out(:,:)
      write(unitbin) hydout(:,:,:,:)
      close(unitbin)

      write(6,*) "output:",nout,time,dt

      nout=nout+1
      tout=time

      return
      end subroutine Output

      subroutine makedirs(outdir)
      implicit none
      character(len=*), intent(in) :: outdir
      character(len=256) command
      write(command, *) 'if [ ! -d ', trim(outdir), ' ]; then mkdir -p ', trim(outdir), '; fi'
      write(*, *) trim(command)
      call system(command)
      end subroutine makedirs
