module forcing_calc_mod

use functions
use constants

contains

subroutine forcing_calc(Nx,Ny,Nz,Nx1,Ny1,Nz1,fu,arr_inc, tau, z, phi)

! calculation of the rhs-tendency of mass continuity equation

implicit none

integer i,j,k,Nx,Ny,Nz,Nx1,Ny1,Nz1
real fu(1:Nx1,1:Ny1,1:Nz1),dt
real z(1:Nz1), phi(1:Ny1)

real*8 arr_inc(:, :, :), tau, z_curr, phi_curr, n_mod, dn_dz, dn_dphi, vphi, dvphi_dphi, vz, dvz_dz

do k=1,Nz1
   do j=1,Ny1
      do i=1,Nx1
        fu(i,j,k)=arr_inc(Nz-k, j, i)/tau


        ! z_curr = z(k)
        ! phi_curr = phi(j)

        ! n_mod = f_m(z_curr, phi_curr)
        ! dn_dz = df_dz(z_curr, phi_curr)
        ! dn_dphi = df_dphi(z_curr, phi_curr)

        ! vphi = vel_phi(z_curr, phi_curr)
        ! dvphi_dphi = dvelphi_dphi(z_curr, phi_curr)
        ! vz = vel_z(z_curr, phi_curr)
        ! dvz_dz = dvelz_dz(z_curr, phi_curr)

        ! fu(i,j,k) = fu(i,j,k) + (dvphi_dphi*n_mod*cos(phi_curr) + vphi*dn_dphi*cos(phi_curr) - vphi*n_mod*sin(phi_curr))/a/cos(phi_curr) + &
        !                          (dvz_dz*n_mod + vz*dn_dz)
      end do
    end do
end do

end subroutine

end module
