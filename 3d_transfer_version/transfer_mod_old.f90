module transfer_mod

use constants
use wind_ini_mod
use forcing_calc_mod
use mass_calc_mod
use cabadv_mod
use flux_ini_mod
use functions

real, allocatable :: lambda(:), lambda_m(:), phi(:), phi_m(:), cphi_m(:), cphi(:), z_m(:), z(:), t(:), ue(:, :, :), c(:,:,:,:,:), v(:,:,:,:), u(:,:,:,:), fu(:,:,:), ux(:,:,:,:), uy(:,:,:,:), uz(:,:,:,:), uef(:,:,:), kurant(:), dphi(:), dlev(:), dlam(:)
real*8, allocatable :: uout(:,:,:)
real*4, allocatable :: uinp(:,:,:)
integer i, j, k, n, m, p, r, nrec

contains

subroutine init_transfer(z1,arr)

    real z1(:)
    real*8 arr(:, :, :)
    integer nrec
    real*8 z_curr, phi_curr

    allocate(lambda(1:Nx), lambda_m(1:Nx+1), phi(1:Ny1), phi_m(1:Ny), cphi_m(1:Ny), cphi(1:Ny1), z_m(1:Nz))
    allocate(c(1:Nx, 1:Ny, 1:Nz, 1:3, 1:2), v(1:Nx1, 1:Ny1, 1:Nz1, 3))
    allocate(u(1:Nx1, 1:Ny1, 1:Nz1, 1:2), fu(1:Nx1, 1:Ny1, 1:Nz1), ux(1:Nx, 1:Ny1, 1:Nz1, 1:2), uy(1:Nx, 1:Ny, 1:Nz1, 1:2))
    allocate(uz(1:Nx, 1:Ny1, 1:Nz, 1:2), kurant(1:3), dphi(1:Ny1), dlev(1:Nz1), dlam(1:Nx))
    allocate(uout(1:Nx1, 1:Ny1, 1:Nz1), uinp(1:Nx1,1:Ny1,1:Nz1))
    allocate(z(1:Nz1))

    z = z1

    nrec=1
    do i = 1, Nx+1
        lambda_m(i) = (i-1)*dx
    end do

    do i = 1, Nx   
        lambda(i) = 0.5*(lambda_m( i ) + lambda_m(i+1))
        dlam(i)   =      lambda_m(i+1) - lambda_m( i )
    end do

    do j = 1, Ny
        phi_m(j) = -pi/2 + (j-1)*dy
    end do

    do j = 1, Ny1
        phi(j)  = 0.5*(phi_m( j ) + phi_m(j+1))
        dphi(j) =      phi_m(j+1) - phi_m( j )
    end do
    z_m(1) = (3*z(1) - z(2)) / 2.

    do k = 2, Nz
        z_m(k) = 2*z(k-1) - z_m(k-1)
    end do

    do k = 1, Nz1
        dlev(k) = z_m(k+1) - z_m(k)
    end do

    cphi   = cos( phi )
    cphi_m = cos(phi_m)

    !open(10, file = 'ion_wind/upb.std', status = 'old', access = 'direct', form = 'unformatted', recl = Nx1*Ny1*Nz1)
    !open(11, file = 'ion_wind/vpb.std', status = 'old', access = 'direct', form = 'unformatted', recl = Nx1*Ny1*Nz1)
    !open(12, file = 'ion_wind/wpb.std', status = 'old', access = 'direct', form = 'unformatted', recl = Nx1*Ny1*Nz1)

    v(:, :, :, 1) = 0
    do k = 1, Nz1
        do j = 1, Ny1
            do i = 1, Nx1
                z_curr = z(k)
                phi_curr = phi(j)
                v(i, j, k, 2) = vel_phi(z_curr, phi_curr)
                v(i, j, k, 3) =   vel_z(z_curr, phi_curr)
            enddo
        enddo
    enddo


    
    do i = 1, Nx1
        do j = 1, Ny1
            do k = 1, Nz1
                u(i, j, k,1) = arr(Nz-k, j, i)
            enddo
        enddo
    enddo
    u(1:Nx1, 1:Ny1, 1:Nz1, 2) = u(1:Nx1, 1:Ny1, 1:Nz1, 1)

    call flux_ini(Nx, Ny, Nz, Nx1, Ny1, Nz1, u(:, :, :, 1), ux(:, :, :, 1), uy(:, :, :, 1), uz(:, :, :, 1), lambda, phi, z, z_m, dlam, dphi, dlev, cphi, pi)

    call wind_ini(Nx, Ny, Nz, Nx1, Ny1, Nz1, v, c(:, :, :, :, 1), lambda_m, z, z_m, eps, pi)
    print*,z
    
    do i = 1, Nx1
        do j = 1, Ny1
            do k = 1, Nz1
                arr(Nz-k, j, i) = u(i, j, k,1)
            enddo
        enddo
    enddo

end subroutine init_transfer




subroutine step_of_transfer( arr_inc,  tau)

    real*8 arr_inc(:, :, :)
    real*8 tau



    dt=tau/Nadv

     call wind_ini(Nx, Ny, Nz, Nx1, Ny1, Nz1, v, c(:, :, :, :, 2), lambda_m, z, z_m, eps, pi)
     ! define Courant numbers

     kurant(1) = dt*maxval(maxval(maxval(abs(c(1:Nx1, 1:Ny1, 1:Nz1, 1, 1)), 1), 2)/(a*cphi(1:Ny1)*dx), 1)
     kurant(2) = dt*maxval(maxval(maxval(abs(c(1:Nx1, 1:Ny1, 1:Nz1, 2, 1)), 1), 2)/(a*dy), 1)
     kurant(3) = dt*maxval(maxval(maxval(abs(c(1:Nx1, 1:Ny1, 1:Nz1, 3, 1)), 1), 1)/(-dlev), 1)

     write(*, '(7a7)'), 'dt(s)', 'dlam', 'dphi', 'dz(km)', 'cx', 'cy', 'cz'
     write(*, '(7f7.3)'), dt, dx*180/pi, dy*180/pi, maxval(-dlev)/1.e3, kurant
     call forcing_calc(Nx, Ny, Nz, Nx1, Ny1, Nz1, fu, arr_inc, tau, z, phi)

    do n = 1, Nadv

        call cab_adv(Nx, Ny, Nz, Nx1, Ny1, Nz1, c, u, fu, ux, uy, uz, dlam, dphi, dlev, phi_m, cphi, cphi_m, lambda, lambda_m, z_m, a, pi, dt)

        u( 1:Nx1, 1:Ny1, 1:Nz1, 1) = u( 1:Nx1, 1:Ny1, 1:Nz1, 2)
        ux(1:Nx1, 1:Ny1, 1:Nz1, 1) = ux(1:Nx1, 1:Ny1, 1:Nz1, 2)
        uy(1:Nx1, 1:Ny , 1:Nz1, 1) = uy(1:Nx1, 1:Ny , 1:Nz1, 2)
        uz(1:Nx1, 1:Ny1, 1:Nz , 1) = uz(1:Nx1, 1:Ny1, 1:Nz , 2)


    end do
    uout(1:Nx1, 1:Ny1, 1:Nz1) = u(1:Nx1, 1:Ny1, 1:Nz1, 1)



end subroutine step_of_transfer

end module transfer_mod
