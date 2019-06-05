module wind_ini_mod
implicit none

contains

subroutine wind_ini(Nx,Ny,Nz,Nx1,Ny1,Nz1,v,c,lambda,lambda_m,phi,z,z_m,u0,w0,alpha,Omega,eps,a,dt,pi,n)
use restore_pwind_mod
implicit none

integer i,j,k,i1,i2,Nx,Ny,Nz,Nx1,Ny1,Nz1,n
real a,pi,u0,w0,alpha,Omega,eps,dt,s1,s2,lambda_0(1:2),b(1:2)
real lambda(1:Nx),lambda_m(1:Nx+1),phi(1:Ny1),z_m(1:Nz),z(1:Nz1),v(1:Nx1,1:Ny1,1:Nz1,1:3),c(1:Nx,1:Ny,1:Nz,1:3)

! define velocities at t=n*dt

!do k=1,Nz1
!  do i=1,Nx1
!     do j=1,Ny1
!      v(i,j,k,1)=u0*(cos(phi(j))*cos(alpha)+sin(phi(j))*cos(lambda(i))*sin(alpha))
!      v(i,j,k,2)=-u0*sin(lambda(i))*sin(alpha)
!      v(i,j,k,3)=w0*cos(Omega*n*dt)
!     end do
!  end do
!end do

do k=1,Nz1
  do i=1,Nx
    do j=1,Ny1
       i2=modulo(i-2,Nx1)+1
       i1=modulo(i-1,Nx1)+1
       c(i,j,k,1)=(v(i1,j,k,1)+v(i2,j,k,1))/2.
    end do
  end do
end do

do k=1,Nz1
  do i=1,Nx1
    do j=2,Ny1
       c(i,j,k,2)=(v(i,j,k,2)+v(i,j-1,k,2))/2.
    end do
  end do
end do

call restore_pwind(Nx,Ny,Nz,Nx1,Ny1,Nz1,c,lambda_m,eps,pi)

  do i=1,Nx1
    do j=1,Ny1
       do k=2,Nz1
          c(i,j,k,3)=((z(k)-z_m(k))*v(i,j,k-1,3)+(z_m(k)-z(k-1))*v(i,j,k,3))/(z(k)-z(k-1))
!           c(i,j,k,3)=(v(i,j,k-1,3)+v(i,j,k,3))/2.
       end do
       c(i,j,1,3)=0.
       c(i,j,Nz,3)=0.
   end do
end do

end subroutine

end module

