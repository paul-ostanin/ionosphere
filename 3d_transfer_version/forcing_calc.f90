module forcing_calc_mod
implicit none

contains

subroutine forcing_calc(Nx,Ny,Nz,Nx1,Ny1,Nz1,c,u,fu,dlam,dphi,dlev,cphi,cphi_m,a,dt, arr_input, arr_prev, tau)

! calculation of the rhs-tendency of mass continuity equation

implicit none

integer i,i1,j,k,Nx,Ny,Nz,Nx1,Ny1,Nz1
real lambda(1:Nx),phi(1:Ny1),dphi(1:Ny1),dlev(1:Nz1),dlam(1:Nx),cphi(1:Ny1),cphi_m(1:Ny),&
& u(1:Nx1,1:Ny1,1:Nz1),fu(1:Nx1,1:Ny1,1:Nz1),c(1:Nx,1:Ny,1:Nz,1:3,1:2)
real a,dt

    real*8 arr_input(:, :, :), arr_prev(:, :, :)
    real*8 tau

do k=1,Nz1
   do j=1,Ny1
      do i=1,Nx1
        i1=modulo(i,Nx1)+1
        fu(i,j,k)=(arr_input(k, j, Nx1-i+1) - arr_prev(k, j, Nx1-i+1))/tau!0.
      end do
    end do
end do

end subroutine

end module
