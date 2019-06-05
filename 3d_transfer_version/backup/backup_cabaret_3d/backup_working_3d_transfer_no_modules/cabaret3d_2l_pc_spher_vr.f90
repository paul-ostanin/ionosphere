! Solve mass continuity equation on the sphere
! 3d spherical layer case, C-grid,
! uniform grid in longitude,latitude and nonuniform in height

program cabaret2d

! library to estimate cpu time
use ifport

implicit none

character(len=*),parameter::maska='ion_upb_trd2'
integer, parameter::&
&                   sp=3,&
&                   ic=20,& 
&                   sf=2**(sp-1),&
&                   Nx=145,&!25*sf+1,&
&                   Ny=91,&!12*sf+1,&
&                   Nz=81,&!4*sf+1,&!2
&                   Nout=24,&!200*sf
&                   Nt=Nout*sf*ic,&
&                   Nx1=Nx-1,&
&                   Ny1=Ny-1,&
&                   Nz1=Nz-1,&
&                   Nnorms=5,&
&                   Nfc=Nt+1
real, parameter::pi=3.1415926,&
&                tend=1.*86400,&
&                h_atm=1.0e5,&
&                a=6.37e6,&
&                dx=2*pi/Nx1,&
&                dy=pi/Ny1,&
&                dz=h_atm/Nz1,&
&                dt=tend/Nt,&
&                u0=2*pi*a/tend,&
&                N_vert_kol=1.,&
&                Omega=2*pi*N_vert_kol/tend,&
&                hraz=0.1*h_atm,&
&                w0=hraz*Omega,&
&                sigma_tot=0.25,&
&                sigma_r=1.5,&
&                sigma_h=0.5,&!&0.05*h_atm,&!0.25,&
&                r0=1./3*a,&
&                rs=r0/3.,&
&                lambda_c=-80.*pi/180.,&
&                phi_c=0.*pi/180.,&
&                h0=0.5*h_atm,&!0.1*h_atm,&
&                hb=h_atm,&
&                z_c=2.5*h_atm,&
&                eps=1e-5,&
&                alpha=0.*pi/180,&
&                beta=-0./tend, &
&                phicr=0.*pi/180
real lambda(1:Nx),lambda_m(1:Nx+1),phi(1:Ny1),phi_m(1:Ny),cphi_m(1:Ny),cphi(1:Ny1),z_m(1:Nz),z(1:Nz1),&
&    t(1:Nt),ue(1:Nx1,1:Ny1,1:Nz1),c(1:Nx,1:Ny,1:Nz,1:3,1:2),sigc(1:Nx,1:Ny,1:Nz,1:3),v(1:Nx1,1:Ny1,1:Nz1,3),&
&    u(1:Nx1,1:Ny1,1:Nz1,1:2),fu(1:Nx1,1:Ny1,1:Nz1),ux(1:Nx,1:Ny1,1:Nz1,1:2),uy(1:Nx,1:Ny,1:Nz1,1:2),uz(1:Nx,1:Ny1,1:Nz,1:2),&
&    uef(1:Nx1,1:Ny1,1:Nz1),ua(1:Nx,1:Ny1,1:Nz1),kurant(1:3),dphi(1:Ny1),dlev(1:Nz1),dlam(1:Nx),lambda_0(1:2),b(1:2)
real*4 res(1:Nnorms),ta(1:2)
real*4 uout(1:Nx1,1:Ny1,1:Nz1),uinp(1:Nz1)

integer i,j,k,n,m,i1,p,q,nrec,ip1,ip2,jq1,jq2,jp1,itr,kr1,kr2,r
real max_val,min_val,mass(3),mom2(2),mass_rel,mom2_rel,time_cpu,up1,up2,uxmax,uxmin,uymax,uymin,uzmax,uzmin,h,rad,rad_tot,s1,s2,cgrid,mass_rel1

character(2) spc
character(1) phic
character(3) lamc

write(spc,'(i2)')ic
write(phic,'(i1)')nint(phi_c*180/pi)
write(lamc,'(i3)')nint(lambda_c*180/pi)

print*,maska


open (20,file='res/res_'//maska//'.dat',status='replace',access='direct',form='unformatted',recl=Nx1*Ny1*Nz1)
open (21,file='res/norms_'//maska//'.txt',status='replace',access='sequential',form='formatted')

! initialize grid

! read vertical levels
open (22,file='ion_wind/zet.dat',status='old',access='sequential',form='formatted')
do i=1,Nz1
read(22,'(i12,f11.2)')k,z(i)
end do
print*,z


do i=1,Nx+1
   lambda_m(i) = (i-1)*dx
end do

do i=1,Nx   
   lambda(i) = 0.5*(lambda_m(i)+lambda_m(i+1))
   dlam(i) = lambda_m(i+1)  -  lambda_m(i)
end do
 cgrid=2.
do j=1,Ny
   phi_m(j) = -pi/2+(j-1)*dy
end do
do j=1,Ny1
   phi(j) = 0.5*(phi_m(j)+phi_m(j+1))
   dphi(j)=phi_m(j+1)-phi_m(j)
end do
z_m(1)=(3*z(1)-z(2))/2.

do k=2,Nz
   z_m(k) = 2*z(k-1)-z_m(k-1)
end do

do k=1,Nz1
   dlev(k)=z_m(k+1)-z_m(k)
end do

do n=1,Nt
   t(n)=(n-1)*dt
end do

 cphi=cos(phi)
 cphi_m=cos(phi_m)

! initialize tracer

!do k=1,Nz1
!  do j=1,Ny1
!     do i=1,Nx1
!!         rad=a*acos(sin(phi_c)*sin(phi(j))+cos(phi_c)*cos(phi(j))*cos(lambda(i)-lambda_c))!/sigma_r
!         rad=a*acos(sin(phi_c)*sin(phi(j))+cos(phi_c)*cos(phi(j)))!/sigma_r
!         h=abs(z(k)-z_c)
!!         rad_tot=sqrt(rad**2+h**2)
!!        ue(i,j,k) = exp(-(rad/sigma_r)**2)*(1+sign(1.,r0-rad))*exp(-(h/sigma_h)**2)*(1+sign(1.,h0-h))/4.
!!         ue(i,j,k) = (1+sign(1.,r0-rad_tot))/2.*exp(-(rad_tot/sigma_tot)**2)!
!!         williamson test
!!         ue(i,j,k) = (1+sign(1.,r0-rad))/2.*(1+cos(pi*rad/r0))/2*(1+sign(1.,h0-h))/2.
!          ue(i,j,k) = (1+sign(1.,r0-rad))/2.*(1+sign(1.,h0-h))/2.
!    end do
!  end do
!end do

open(9,file='ion_wind/res_Nz_80_Nphi_90_tau_5sec-1.txt',status='old',access='sequential',form='formatted')

do j=1,Ny1
    read(9,'(80f10.7)')uinp
!    print*,'j=',j,uinp
    do k=1,Nz1
       ue(:,j,k)=uinp(Nz1-k+1)
    end do    
end do
do i=1,Nx1
   ue(i,:,:)=ue(i,:,:)*sin(lambda(i)-lambda_c)**2
end do

u(1:Nx1,1:Ny1,1:Nz1,1)=ue(1:Nx1,1:Ny1,1:Nz1)
u(1:Nx1,1:Ny1,1:Nz1,2)=u(1:Nx1,1:Ny1,1:Nz1,1)

! read wind

open(10,file='ion_wind/upb.std',status='old',access='direct',form='unformatted',recl=Nx1*Ny1*Nz1)
open(11,file='ion_wind/vpb.std',status='old',access='direct',form='unformatted',recl=Nx1*Ny1*Nz1)
open(12,file='ion_wind/wpb.std',status='old',access='direct',form='unformatted',recl=Nx1*Ny1*Nz1)

read(10,rec=1)uout
v(:,:,:,1)=uout
read(11,rec=1)uout
v(:,:,:,2)=uout
read(12,rec=1)uout
v(:,:,:,3)=uout

!initialization of tracer variables and velocity components
n=0

! define velocities
 call wind_ini(Nx,Ny,Nz,Nx1,Ny1,Nz1,v,c(:,:,:,:,1),lambda,lambda_m,phi,z,z_m,u0,w0,alpha,Omega,eps,a,dt,pi,n)            
! c(:,:,:,:,2)=c(:,:,:,:,1)

!   to define tracer at cell boundaries (flux variables)

call flux_ini(Nx,Ny,Nz,Nx1,Ny1,Nz1,u(:,:,:,1),ux(:,:,:,1),uy(:,:,:,1),uz(:,:,:,1),lambda,phi,z,z_m,dlam,dphi,dlev,cphi,pi)

uef(1:Nx1,1:Ny1,1:Nz1)=u(1:Nx1,1:Ny1,1:Nz1,1)

! define Courant numbers

kurant(1)=dt*maxval(maxval(maxval(abs(c(1:Nx1,1:Ny1,1:Nz1,1,1)),1),2)/(a*cphi(1:Ny1)*dx),1)
kurant(2)=dt*maxval(maxval(maxval(abs(c(1:Nx1,1:Ny1,1:Nz1,2,1)),1),2)/(a*dy),1)
kurant(3)=dt*maxval(maxval(maxval(abs(c(1:Nx1,1:Ny1,1:Nz1,3,1)),1),1)/(-dlev),1)

write(*,'(7a7)'),'dt(s)','dlam','dphi','dz(km)','cx','cy','cz'
write(*,'(7f7.3)'),dt,dx*180/pi,dy*180/pi,maxval(-dlev)/1.e3,kurant

write(*,'(3a10)'),'lambda0','phi0','z0(km)'
write(*,'(3f10.1)'),lambda_c*180/pi,phi_c*180/pi,z_c/1.e3

! print*,dt,dx,dy,dz

! write output to file

nrec=1
uout(1:Nx1,1:Ny1,1:Nz1)=u(1:Nx1,1:Ny1,1:Nz1,1)
write(20,rec=nrec)uout

! to estimate cpu time

time_cpu=dtime(ta)

! main loop on time
do n=1,Nt

 call wind_ini(Nx,Ny,Nz,Nx1,Ny1,Nz1,v,c(:,:,:,:,2),lambda,lambda_m,phi,z,z_m,u0,w0,alpha,Omega,eps,a,dt,pi,n)

 call forcing_calc(Nx,Ny,Nz,Nx1,Ny1,Nz1,c,u,fu,dlam,dphi,dlev,cphi,cphi_m,a,dt)

 call cab_adv(Nx,Ny,Nz,Nx1,Ny1,Nz1,c,u,fu,ux,uy,uz,dlam,dphi,dlev,phi_m,cphi,cphi_m,lambda,lambda_m,z_m,a,pi,dt)

  u(1:Nx1,1:Ny1,1:Nz1,1)=u(1:Nx1,1:Ny1,1:Nz1,2)
  ux(1:Nx1,1:Ny1,1:Nz1,1)=ux(1:Nx1,1:Ny1,1:Nz1,2)
  uy(1:Nx1,1:Ny,1:Nz1,1)=uy(1:Nx1,1:Ny,1:Nz1,2)
  uz(1:Nx1,1:Ny1,1:Nz,1)=uz(1:Nx1,1:Ny1,1:Nz,2)
  c(1:Nx1,1:Ny,1:Nz,1:3,1)=c(1:Nx1,1:Ny,1:Nz,1:3,2)

! write output to file

if (mod(n,Nt/Nout)==0) then
nrec=nrec+1
uout(1:Nx1,1:Ny1,1:Nz1)=u(1:Nx1,1:Ny1,1:Nz1,1)
write(20,rec=nrec)uout
print*,'Nout=',n*Nout/Nt
end if

end do

time_cpu=dtime(ta)

! estimation norms of error field

max_val=maxval(u(1:Nx1,1:Ny1,1:Nz1,1))/maxval(ue(1:Nx1,1:Ny1,1:Nz1))
min_val=minval(u(1:Nx1,1:Ny1,1:Nz1,1))/maxval(ue(1:Nx1,1:Ny1,1:Nz1))
call mass_calc(Nx,Ny,Nz,Nx1,Ny1,Nz1,ue,dlam,dphi,dlev,cphi,mass(1))
call mass_calc(Nx,Ny,Nz,Nx1,Ny1,Nz1,u(:,:,:,1),dlam,dphi,dlev,cphi,mass(2))
call mass_calc(Nx,Ny,Nz,Nx1,Ny1,Nz1,ue**2,dlam,dphi,dlev,cphi,mom2(1))
call mass_calc(Nx,Ny,Nz,Nx1,Ny1,Nz1,u(:,:,:,1)**2,dlam,dphi,dlev,cphi,mom2(2))

mass_rel=mass(2)/mass(1)
mom2_rel=mom2(2)/mom2(1)

res(1)=mass_rel
res(2)=mom2_rel
res(3)=min_val
res(4)=max_val
res(5)=time_cpu/60.


write(*,'(5a12)'),'Mass_rel','Mom2_rel','Min_rel','Max_rel','CPU(min)'
write(*,'(5f12.5)'),res

! write output to file

write(21,'(5f12.8)')res

end
