module wind_ini_mod
implicit none

contains

subroutine wind_ini(Nx,Ny,Nz,Nx1,Ny1,Nz1,v,c,lambda_m,z,z_m,eps,pi)
use restore_pwind_mod
implicit none

integer i,j,k,i1,i2,Nx,Ny,Nz,Nx1,Ny1,Nz1
real pi,eps
real lambda_m(1:Nx+1),v(1:Nx1,1:Ny1,1:Nz1,1:3),c(1:Nx,1:Ny,1:Nz,1:3),z(1:Nz1),z_m(1:Nz)

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
       c(i,j,1,3)=v(i,j,1,3)
       c(i,j,Nz,3)=v(i,j,Nz1,3)
   end do
end do

end subroutine

end module

