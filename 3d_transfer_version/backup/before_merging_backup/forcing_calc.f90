module forcing_calc_mod

contains

subroutine forcing_calc(Nx,Ny,Nz,Nx1,Ny1,Nz1,fu,arr_inc, dt)

! calculation of the rhs-tendency of mass continuity equation

implicit none

integer i,j,k,Nx,Ny,Nz,Nx1,Ny1,Nz1
real fu(1:Nx1,1:Ny1,1:Nz1),dt

    real*8 arr_inc(:, :, :)

do k=1,Nz1
   do j=1,Ny1
      do i=1,Nx1
        fu(i,j,k)=arr_inc(Nz-k, j, i)/dt
      end do
    end do
end do

end subroutine

end module
