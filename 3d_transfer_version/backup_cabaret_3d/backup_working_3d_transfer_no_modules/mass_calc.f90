subroutine mass_calc(Nx,Ny,Nz,Nx1,Ny1,Nz1,u,dlam,dphi,dlev,cphi,mass)

! calculation of the tracer mass

implicit none

integer i,j,k,Nx,Ny,Nz,Nx1,Ny1,Nz1
real lambda(1:Nx),phi(1:Ny1),dphi(1:Ny1),dlev(1:Nz1),dlam(1:Nx),cphi(1:Ny1),&
& u(1:Nx1,1:Ny1,1:Nz1),mass

mass=0.
do k=1,Nz1
   do j=1,Ny1
      do i=1,Nx1
        mass=mass+u(i,j,k)*cphi(j)*dphi(j)*dlam(i)*(-dlev(k))
      end do
    end do
end do

!print*,mass

end 
