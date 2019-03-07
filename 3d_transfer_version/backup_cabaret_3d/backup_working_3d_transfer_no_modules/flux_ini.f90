subroutine flux_ini(Nx,Ny,Nz,Nx1,Ny1,Nz1,u,ux,uy,uz,lambda,phi,z,z_m,dlam,dphi,dlev,cphi,pi)

! interpolate function on the flux variables grid

implicit none

logical,parameter::mcorr=.true.
integer i,j,k,i0,i1,Nx,Ny,Nz,Nx1,Ny1,Nz1
real pi,mr
real lambda(1:Nx),phi(1:Ny1),z_m(1:Nz),z(1:Nz1),dphi(1:Ny1),dlev(1:Nz1),dlam(1:Nx),cphi(1:Ny1),&
& u(1:Nx1,1:Ny1,1:Nz1),ux(1:Nx,1:Ny1,1:Nz1),uy(1:Nx,1:Ny,1:Nz1),uz(1:Nx,1:Ny1,1:Nz),mass(1:2)


do k=1,Nz1
  do i=1,Nx
    do j=1,Ny1
       i1=modulo(i-2,Nx1)+1
       i0=modulo(i-1,Nx1)+1
       ux(i,j,k)=(u(i0,j,k)+u(i1,j,k))/2.
    end do
end do
end do

do k=1,Nz1
  do i=1,Nx1
    do j=2,Ny1
          uy(i,j,k)=(u(i,j,k)+u(i,j-1,k))/2.
    end do
  end do
  
!   to define tracer at the pole point 
   j=Ny
   uy(1:Nx1,j,k)=sum(uy(1:Nx1,j-1,k)*(lambda(2:Nx)-lambda(1:Nx1)))/(2*pi) 
   j=1
   uy(1:Nx1,j,k)=sum(uy(1:Nx1,j+1,k)*(lambda(2:Nx)-lambda(1:Nx1)))/(2*pi)
end do

do i=1,Nx1
  do j=1,Ny1
    do k=2,Nz1
       uz(i,j,k)=((z(k)-z_m(k))*u(i,j,k-1)+(z_m(k)-z(k-1))*u(i,j,k))/(z(k)-z(k-1))
    end do
    uz(i,j,1)=u(i,j,1)
    uz(i,j,Nz)=u(i,j,Nz1)
  end do
end do

if (mcorr) then
   call mass_calc(Nx,Ny,Nz,Nx1,Ny1,Nz1,u,dlam,dphi,dlev,cphi,mass(1))      
end if

do k=1,Nz1
  do j=1,Ny1
     do i=1,Nx1
        u(i,j,k)=1/6.*(ux(i,j,k)+ux(i+1,j,k)+uy(i,j,k)+uy(i,j+1,k)+uz(i,j,k+1)+uz(i,j,k))
    end do   
  end do
end do

if (mcorr) then
   call mass_calc(Nx,Ny,Nz,Nx1,Ny1,Nz1,u,dlam,dphi,dlev,cphi,mass(2))
mr=mass(1)/mass(2)
u=mr*u
ux=mr*ux
uy=mr*uy
uz=mr*uz

print*,'mass',mass

end if



end
