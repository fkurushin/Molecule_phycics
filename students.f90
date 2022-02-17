module constants
  ! DOI: 10.1103/PhysRevB.71.035211
  implicit real*8 (a-z)
  real*8, parameter:: D0=3.24d0
  real*8, parameter:: r0=2.222d0
  real*8, parameter:: S=1.57d0
  real*8, parameter:: betta=1.4760d0
  real*8, parameter:: gamma=0.09253d0
  real*8, parameter:: c=1.13681d0
  real*8, parameter:: d=0.63397d0
  real*8, parameter:: h=0.335d0
  real*8, parameter:: R=2.90d0
  real*8, parameter:: DD=0.15d0
end module constants

module coords
  integer nat ! number of atoms
  integer, parameter:: Nmax=10000
  real*8 x(Nmax),y(Nmax), z(Nmax) ! coordinates
  integer atnumber(Nmax) ! atomic number
  real*8 fx(Nmax),fy(Nmax),fz(Nmax)
end module coords

subroutine set_coords
use coords
  nat=5
  x(1)=0.; y(1)=0.; z(1)=0.;

  x(2)=2.; y(2)=0.; z(2)=0.;
  x(3)=4.; y(3)=0.5; z(3)=0.;
  x(4)=4.7; y(4)=1.; z(4)=0.5;
  x(5)=5.; y(5)=3.; z(5)=1.;  
end subroutine set_coords

subroutine out_coords
use coords
  open(1,file='coords.xyz',access='append')
  write(1,*) nat
  write(1,*) 'coords'
  do i=1,nat
    write(1,*) 'Si  ',x(i),y(i),z(i)
  enddo
  close(1)
end subroutine out_coords

subroutine Verlet
use coords
implicit real*8 (a-z)
real*8 time,dtime,vfx(nat),vfy(nat),vfz(nat),vx(nat),vy(nat),vz(nat),massa,kT
integer i,k
  call set_rand(1885,1986)
  temp =300.d0; kT=temp/11602.d0; dtemp =580.d0;  
  time=0.d0; dtime=1.0; massa=2925.6d0 ! time in fs
  do i=1,nat
    vx(i)=gaus(0.d0, kT/massa)
    vy(i)=gaus(0.d0, kT/massa)
    vz(i)=gaus(0.d0, kT/massa)
  enddo  
  call a_force
  do k=1,time
    do i=1,nat
        x(i)=x(i)+vx(i)*dtime+0.5d0*fx(i)*dtime*dtime/massa
        y(i)=y(i)+vy(i)*dtime+0.5d0*fy(i)*dtime*dtime/massa
        z(i)=z(i)+vz(i)*dtime+0.5d0*fz(i)*dtime*dtime/massa 
        vfx(i)=fx(i); vfy(i)=fy(i); vfz(i)=fz(i)
     enddo
     call a_force
     do i=1,nat
         vx(i)=vx(i)+0.5d0*dtime*(vfx(i)+fx(i))/massa
         vy(i)=vy(i)+0.5d0*dtime*(vfy(i)+fy(i))/massa
         vz(i)=vz(i)+0.5d0*dtime*(vfz(i)+fz(i))/massa
     enddo
end subroutine Verlet

program main
use coords
real*8 e; integer i
  call set_coords
  open(1,file='coords.xyz'); write(1,*) ''; close(1)
  call out_coords
  call energy(e)
  print*, e; read*
  call relax
end program main

subroutine a_force
use coords
use constants
real*8 e,fcij,fcik,rij,rik,Vr,Va,b,cosQ,Gijk,Gijk1
! *********************
real*8 fcij1x,fcik1x,Vr1x,Va1x,b1x(nat),cosQ1ix,cosQ1jx,cosQ1kx
real*8 fcij1y,fcik1y,Vr1y,Va1y,b1y(nat),cosQ1iy,cosQ1jy,cosQ1ky
real*8 fcij1z,fcik1z,Vr1z,Va1z,b1z(nat),cosQ1iz,cosQ1jz,cosQ1kz
! *********************
integer i,j,k
  do i=1,nat
    fx(i)=0.d0; fy(i)=0.d0; fz(i)=0.d0
  enddo
  do i=1,nat
  do j=1,nat
  if (i.ne.j) then
       rij=dsqrt((x(i)-x(j))**2.d0+(y(i)-y(j))**2.d0+(z(i)-z(j))**2.d0)
       if (rij.lt.(R+DD)) then
          fcij=1.d0; fcij1x=0.d0; fcij1y=0.d0; fcij1z=0.d0
          if (rij.gt.(R-DD)) then
              fcij=0.5d0-0.5d0*dsin(0.5d0*3.14159265d0*(rij-R)/DD)
              fcij1x=-0.5d0*dcos(0.5d0*3.14159265d0*(rij-R)/DD)*0.5d0*3.14159265d0*(x(i)-x(j))/(rij*DD)
              fcij1y=-0.5d0*dcos(0.5d0*3.14159265d0*(rij-R)/DD)*0.5d0*3.14159265d0*(y(i)-y(j))/(rij*DD)
              fcij1z=-0.5d0*dcos(0.5d0*3.14159265d0*(rij-R)/DD)*0.5d0*3.14159265d0*(z(i)-z(j))/(rij*DD)
          endif
          Vr=D0*dexp(-betta*dsqrt(2.d0*S)*(rij-r0))/(S-1.d0)
          Vr1x=Vr*(-betta)*dsqrt(2.d0*S)*(x(i)-x(j))/rij
          Vr1y=Vr*(-betta)*dsqrt(2.d0*S)*(y(i)-y(j))/rij
          Vr1z=Vr*(-betta)*dsqrt(2.d0*S)*(z(i)-z(j))/rij
          Va=S*D0*dexp(-betta*dsqrt(2.d0/S)*(rij-r0))/(S-1.d0)
          Va1x=Va*(-betta)*dsqrt(2.d0/S)*(x(i)-x(j))/rij
          Va1y=Va*(-betta)*dsqrt(2.d0/S)*(y(i)-y(j))/rij
          Va1z=Va*(-betta)*dsqrt(2.d0/S)*(z(i)-z(j))/rij
          b=1.d0
          do k=1,nat
          if (((k-i)*(k-j)).ne.0) then
             rik=dsqrt((x(i)-x(k))**2.d0+(y(i)-y(k))**2.d0+(z(i)-z(k))**2.d0)
             if (rik.lt.(R+DD)) then
                fcik=1.d0
                if (rik.gt.(R-DD)) then
                   fcik=0.5d0-0.5d0*dsin(0.5d0*3.14159265d0*(rik-R)/DD)
                 endif
                 cosQ=(x(j)-x(i))*(x(k)-x(i))+(y(j)-y(i))*(y(k)-y(i))+(z(j)-z(i))*(z(k)-z(i))
                 cosQ=cosQ/(rij*rik)
                 b=b+fcik*gamma*(1.d0+c*c/(d*d)-c*c/(d*d+(h+cosQ)*(h+cosQ)))
             endif
          endif            
          enddo
          b=1.d0/dsqrt(b)
!*****************************************
          b1x(1:nat)=0.d0; b1y(1:nat)=0.d0; b1z(1:nat)=0.d0
          do k=1,nat
          if (((k-i)*(k-j)).ne.0) then
             rik=dsqrt((x(i)-x(k))**2.d0+(y(i)-y(k))**2.d0+(z(i)-z(k))**2.d0)
             if (rik.lt.(R+DD)) then
                fcik=1.d0; fcik1x=0.d0; fcik1y=0.d0; fcik1z=0.d0
                if (rik.gt.(R-DD)) then
                   fcik=0.5d0-0.5d0*dsin(0.5d0*3.14159265d0*(rik-R)/DD)
                   fcik1x=0.5d0*dcos(0.5d0*3.14159265d0*(rik-R)/DD)*0.5d0*3.14159265d0*(x(i)-x(k))/(rik*DD)
                   fcik1y=0.5d0*dcos(0.5d0*3.14159265d0*(rik-R)/DD)*0.5d0*3.14159265d0*(y(i)-y(k))/(rik*DD)
                   fcik1z=0.5d0*dcos(0.5d0*3.14159265d0*(rik-R)/DD)*0.5d0*3.14159265d0*(z(i)-z(k))/(rik*DD)
                 endif
                 cosQ=(x(j)-x(i))*(x(k)-x(i))+(y(j)-y(i))*(y(k)-y(i))+(z(j)-z(i))*(z(k)-z(i))
                 cosQ=cosQ/(rij*rik)
                 cosQ1ix=(2.d0*x(i)-x(j)-x(k))-((x(i)-x(k))*(rij/rik)+(x(i)-x(j))*(rik/rij))*cosQ; cosQ1ix=cosQ1ix/(rij*rik)
                 cosQ1jx=(x(k)-x(i))*rij-(x(j)-x(i))*cosQ*rik; cosQ1jx=cosQ1jx/(rik*rij*rij)
                 cosQ1kx=(x(j)-x(i))*rik-(x(k)-x(i))*cosQ*rij; cosQ1kx=cosQ1kx/(rij*rik*rik)
                 cosQ1iy=(2.d0*y(i)-y(j)-y(k))-((y(i)-y(k))*(rij/rik)+(y(i)-y(j))*(rik/rij))*cosQ; cosQ1iy=cosQ1iy/(rij*rik)
                 cosQ1jy=(y(k)-y(i))*rij-(y(j)-y(i))*cosQ*rik; cosQ1jy=cosQ1jy/(rik*rij*rij)
                 cosQ1ky=(y(j)-y(i))*rik-(y(k)-y(i))*cosQ*rij; cosQ1ky=cosQ1ky/(rij*rik*rik)
                 cosQ1iz=(2.d0*z(i)-z(j)-z(k))-((z(i)-z(k))*(rij/rik)+(z(i)-z(j))*(rik/rij))*cosQ; cosQ1iz=cosQ1iz/(rij*rik)
                 cosQ1jz=(z(k)-z(i))*rij-(z(j)-z(i))*cosQ*rik; cosQ1jz=cosQ1jz/(rik*rij*rij)
                 cosQ1kz=(z(j)-z(i))*rik-(z(k)-z(i))*cosQ*rij; cosQ1kz=cosQ1kz/(rij*rik*rik)
                 Gijk=gamma*(1.d0+c*c/(d*d)-c*c/(d*d+(h+cosQ)*(h+cosQ)))
                 Gijk1=-gamma*c*c*2.d0*(h+cosQ)/((d*d+(h+cosQ)*(h+cosQ))**2.d0)
                 b1x(i)=b1x(i)+0.5d0*(b*b*b)*(fcik1x*Gijk+fcik*Gijk1*cosQ1ix)
                 b1x(j)=b1x(j)+0.5d0*(b*b*b)*fcik*Gijk1*cosQ1jx
                 b1x(k)=b1x(k)+0.5d0*(b*b*b)*(-fcik1x*Gijk+fcik*Gijk1*cosQ1kx)
                 b1y(i)=b1y(i)+0.5d0*(b*b*b)*(fcik1y*Gijk+fcik*Gijk1*cosQ1iy)
                 b1y(j)=b1y(j)+0.5d0*(b*b*b)*fcik*Gijk1*cosQ1jy
                 b1y(k)=b1y(k)+0.5d0*(b*b*b)*(-fcik1y*Gijk+fcik*Gijk1*cosQ1ky)
                 b1z(i)=b1z(i)+0.5d0*(b*b*b)*(fcik1z*Gijk+fcik*Gijk1*cosQ1iz)
                 b1z(j)=b1z(j)+0.5d0*(b*b*b)*fcik*Gijk1*cosQ1jz
                 b1z(k)=b1z(k)+0.5d0*(b*b*b)*(-fcik1z*Gijk+fcik*Gijk1*cosQ1kz)
              endif
          endif           
          enddo     
          fx(i)=fx(i)-0.5d0*fcij*(Vr1x-b*Va1x)-0.5d0*(Vr-b*Va)*fcij1x; 
          fx(j)=fx(j)+0.5d0*fcij*(Vr1x-b*Va1x)+0.5d0*(Vr-b*Va)*fcij1x;
          fy(i)=fy(i)-0.5d0*fcij*(Vr1y-b*Va1y)-0.5d0*(Vr-b*Va)*fcij1y; 
          fy(j)=fy(j)+0.5d0*fcij*(Vr1y-b*Va1y)+0.5d0*(Vr-b*Va)*fcij1y;
          fz(i)=fz(i)-0.5d0*fcij*(Vr1z-b*Va1z)-0.5d0*(Vr-b*Va)*fcij1z; 
          fz(j)=fz(j)+0.5d0*fcij*(Vr1z-b*Va1z)+0.5d0*(Vr-b*Va)*fcij1z;
          do i1=1,nat
            fx(i1)=fx(i1)+0.5d0*fcij*b1x(i1)*Va
            fy(i1)=fy(i1)+0.5d0*fcij*b1y(i1)*Va
            fz(i1)=fz(i1)+0.5d0*fcij*b1z(i1)*Va
          enddo     
!***************************************
       endif
  endif
  enddo
  enddo
end subroutine a_force

subroutine n_force
use coords
real*8 d,e0,e1
integer i
  d=1.d-4
  call energy(e0)
  print*, 'NUMERICAL FORCES:'
  do i=1,nat
     x(i)=x(i)+d; call energy(e1); x(i)=x(i)-d; fx(i)=-(e1-e0)/d
     y(i)=y(i)+d; call energy(e1); y(i)=y(i)-d; fy(i)=-(e1-e0)/d
     z(i)=z(i)+d; call energy(e1); z(i)=z(i)-d; fz(i)=-(e1-e0)/d
     ! write(*,'(I5,3F10.4)') i,fx(i),fy(i),fz(i)
  enddo
end subroutine n_force

subroutine relax
use coords
real*8 e,step
  step=0.005d0
  do
    call a_force
    do i=1,nat
      x(i)=x(i)+step*fx(i); y(i)=y(i)+step*fy(i); z(i)=z(i)+step*fz(i)
    enddo
    call energy(e); print*, e;
    call out_coords
  enddo
end subroutine relax

subroutine energy(e)
use constants
use coords
real*8 e,fcij,fcik,rij,rik,Vr,Va,b,cosQ
integer i,j,k
  e=0.d0
  do i=1,nat
  do j=1,nat
  if (i.ne.j) then
       rij=dsqrt((x(i)-x(j))**2.d0+(y(i)-y(j))**2.d0+(z(i)-z(j))**2.d0)
       if (rij.lt.(R+DD)) then
          fcij=1.d0
          if (rij.gt.(R-DD)) then
              fcij=0.5d0-0.5d0*dsin(0.5d0*3.14159265d0*(rij-R)/DD)
          endif
          Vr=D0*dexp(-betta*dsqrt(2.d0*S)*(rij-r0))/(S-1.d0)
          Va=S*D0*dexp(-betta*dsqrt(2.d0/S)*(rij-r0))/(S-1.d0)
          b=1.d0
          do k=1,nat
          if (((k-i)*(k-j)).ne.0) then
             rik=dsqrt((x(i)-x(k))**2.d0+(y(i)-y(k))**2.d0+(z(i)-z(k))**2.d0)
             if (rik.lt.(R+DD)) then
                fcik=1.d0
                if (rik.gt.(R-DD)) then
                   fcik=0.5d0-0.5d0*dsin(0.5d0*3.14159265d0*(rik-R)/DD)
                 endif
                 cosQ=(x(j)-x(i))*(x(k)-x(i))+(y(j)-y(i))*(y(k)-y(i))+(z(j)-z(i))*(z(k)-z(i))
                 cosQ=cosQ/(rij*rik)
                 b=b+fcik*gamma*(1.d0+c*c/(d*d)-c*c/(d*d+(h+cosQ)*(h+cosQ)))
              endif
          endif            
          enddo
          b=1.d0/dsqrt(b)
          e=e+0.5d0*fcij*(Vr-b*Va)
       endif
  endif   
  enddo
  enddo
end subroutine energy

real*8 function gaus(centre,disp) 
  implicit real*8(a-z)
  real*8 centre,disp,ro,x,pi,r  
  pi=3.14159265359d0
  r=rndm()
  ro=dsqrt(-2.d0*dlog(r))
  r=rndm()
  x=ro*dcos(2.d0*pi*r)
  x=centre+x*dsqrt(2.d0*disp)
  gaus=x
end function gaus
 
subroutine set_rand(IJ,KL)
   integer ij,kl,i97,j97,i7,k7,j7,l7,ii7,jj7,m7
   double precision c, cd,cm,r1,r2,s,tt,ugen
   common/UGEN/UGEN(97),C,CD,CM
   common/UGEN1/r1,r2
   common/UGEN2/I97,J97
   if ( IJ<0 .OR. IJ>31328 .OR. KL<0 .OR. KL>30081) then
      print '(A)',' The first seed must be between 0 and 31328'
      print '(A)',' The second seed must be between 0 and 30081'
      stop
    endif
    I7=MOD(IJ/177,177)+2 
    J7=MOD(IJ,177)+2
    K7=MOD(KL/169,178)+1
    L7=MOD(KL,169)
    do ii7=1,97
       S=0.D0
       Tt=0.5D0
       do jj7=1,24
          M7=MOD(MOD(I7*J7,179)*K7,179)
          I7=J7
          J7=K7 
          K7=M7 
          L7=MOD(53*L7+1,169)
          if (MOD(L7*M7,64)>32) S=S+Tt
          Tt=0.5D0*Tt
        enddo
        UGEN(ii7)=S
     enddo
     C=362436.D0/16777216.D0
     CD=7654321.D0/16777216.D0
     CM=16777213.D0/16777216.D0
     I97=97
     J97=33
     r1=5.d-15
     r2=1.d-14
     return
end subroutine set_rand
      
real*8 function rndm()
   double precision c, cd,cm,r1,r2,ugen,rval
   common/UGEN/UGEN(97), C, CD, CM
   common/UGEN1/r1,r2
   common/UGEN2/I97, J97
   integer i97,j97
   RVAL=UGEN(I97)-UGEN(J97)
   if (RVAL<0.D0) RVAL=RVAL+1.0d0
   UGEN(I97)=RVAL 
   I97=I97-1
   if (I97==0) I97=97
   J97=J97-1; if (J97 == 0) J97=97
   C=C-CD; if (C < 0.D0 ) C=C+CM
   RVAL=RVAL-C
   if ( RVAL .LT. 0.D0 ) RVAL = RVAL + 1.0d0
   rndm=max(RVAL-r1,r2)
   return
end function rndm