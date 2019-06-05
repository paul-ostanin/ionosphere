module transfer_mod

use constants
use wind_ini_mod
use forcing_calc_mod
use mass_calc_mod
use cabadv_mod
use flux_ini_mod

real, allocatable :: lambda(:), lambda_m(:), phi(:), phi_m(:), cphi_m(:), cphi(:), z_m(:), z(:), t(:), ue(:,:,:), c(:,:,:,:,:), v(:,:,:,:), u(:,:,:,:), fu(:,:,:), ux(:,:,:,:), uy(:,:,:,:), uz(:,:,:,:), uef(:,:,:), kurant(:), dphi(:), dlev(:), dlam(:)
real*4, allocatable :: res(:), uout(:,:,:), uinp(:)
integer i, j, k, n, m, p, nrec, r
real max_val, min_val, mass(3), mom2(2), mass_rel, mom2_rel, h, rad, rad_tot, s1, s2, cgrid, mass_rel1

character(2) spc
character(1) phic
character(3) lamc


contains

subroutine init_transfer()

    allocate(lambda(1:Nx), lambda_m(1:Nx+1), phi(1:Ny1), phi_m(1:Ny), cphi_m(1:Ny), cphi(1:Ny1), z_m(1:Nz), z(1:Nz1))
    allocate(t(1:Nt), ue(1:Nx1, 1:Ny1, 1:Nz1), c(1:Nx, 1:Ny, 1:Nz, 1:3, 1:2), v(1:Nx1, 1:Ny1, 1:Nz1, 3))
    allocate(u(1:Nx1, 1:Ny1, 1:Nz1, 1:2), fu(1:Nx1, 1:Ny1, 1:Nz1), ux(1:Nx, 1:Ny1, 1:Nz1, 1:2), uy(1:Nx, 1:Ny, 1:Nz1, 1:2))
    allocate(uz(1:Nx, 1:Ny1, 1:Nz, 1:2), uef(1:Nx1, 1:Ny1, 1:Nz1), kurant(1:3), dphi(1:Ny1), dlev(1:Nz1), dlam(1:Nx))
    allocate(res(1:Nnorms), uout(1:Nx1, 1:Ny1, 1:Nz1), uinp(1:Nz1))



end subroutine init_transfer

subroutine step_of_transfer()


    write(spc, '(i2)') ic
    write(phic, '(i1)') nint(phi_c*180/pi)
    write(lamc, '(i3)') nint(lambda_c*180/pi)

    print*, maska

    open (20, file = 'res/res_'//maska//'.dat', status = 'replace', access = 'direct', form = 'unformatted', recl = Nx1*Ny1*Nz1)
    open (21, file = 'res/norms_'//maska//'.txt', status = 'replace', access = 'sequential', form = 'formatted')

    ! initialize grid

    ! read vertical levels
    open (22, file = 'ion_wind/zet.dat', status = 'old', access = 'sequential', form = 'formatted')
    do i = 1, Nz1
        read(22, '(i12, f11.2)') k, z(i)
    end do
    print*, z


    do i = 1, Nx+1
        lambda_m(i)  =  (i-1)*dx
    end do

    do i = 1, Nx   
        lambda(i)  =  0.5*(lambda_m(i)+lambda_m(i+1))
        dlam(i)  =  lambda_m(i+1)  -  lambda_m(i)
    end do

    cgrid = 2.

    do j = 1, Ny
        phi_m(j)  =  -pi/2+(j-1)*dy
    end do

    do j = 1, Ny1
        phi(j)  =  0.5*(phi_m(j)+phi_m(j+1))
        dphi(j) = phi_m(j+1)-phi_m(j)
    end do
    z_m(1) = (3*z(1)-z(2))/2.

    do k = 2, Nz
        z_m(k)  =  2*z(k-1)-z_m(k-1)
    end do

    do k = 1, Nz1
        dlev(k) = z_m(k+1)-z_m(k)
    end do

    do n = 1, Nt
        t(n) = (n-1)*dt
    end do

    cphi   = cos( phi )
    cphi_m = cos(phi_m)

    ! initialize tracer

    open(9, file = 'ion_wind/res_Nz_80_Nphi_90_tau_5sec-1.txt', status = 'old', access = 'sequential', form = 'formatted')

    do j = 1, Ny1
        read(9, '(80f10.7)')uinp
        do k = 1, Nz1
            ue(:, j, k) = uinp(Nz1-k+1)
        end do    
    end do
    do i = 1, Nx1
        ue(i, :, :) = ue(i, :, :)*sin(lambda(i)-lambda_c)**2
    end do

    u(1:Nx1, 1:Ny1, 1:Nz1, 1) = ue(1:Nx1, 1:Ny1, 1:Nz1)
    u(1:Nx1, 1:Ny1, 1:Nz1, 2) = u(1:Nx1, 1:Ny1, 1:Nz1, 1)

    ! read wind

    open(10, file = 'ion_wind/upb.std', status = 'old', access = 'direct', form = 'unformatted', recl = Nx1*Ny1*Nz1)
    open(11, file = 'ion_wind/vpb.std', status = 'old', access = 'direct', form = 'unformatted', recl = Nx1*Ny1*Nz1)
    open(12, file = 'ion_wind/wpb.std', status = 'old', access = 'direct', form = 'unformatted', recl = Nx1*Ny1*Nz1)

    read(10, rec = 1)uout
    v(:, :, :, 1) = uout
    read(11, rec = 1)uout
    v(:, :, :, 2) = uout
    read(12, rec = 1)uout
    v(:, :, :, 3) = uout

    !initialization of tracer variables and velocity components
    n = 0

    ! define velocities
     call wind_ini(Nx, Ny, Nz, Nx1, Ny1, Nz1, v, c(:, :, :, :, 1), lambda, lambda_m, phi, z, z_m, u0, w0, alpha, Omega, eps, a, dt, pi, n)            

    !   to define tracer at cell boundaries (flux variables)

    call flux_ini(Nx, Ny, Nz, Nx1, Ny1, Nz1, u(:, :, :, 1), ux(:, :, :, 1), uy(:, :, :, 1), uz(:, :, :, 1), lambda, phi, z, z_m, dlam, dphi, dlev, cphi, pi)

    uef(1:Nx1, 1:Ny1, 1:Nz1) = u(1:Nx1, 1:Ny1, 1:Nz1, 1)

    ! define Courant numbers

    kurant(1) = dt*maxval(maxval(maxval(abs(c(1:Nx1, 1:Ny1, 1:Nz1, 1, 1)), 1), 2)/(a*cphi(1:Ny1)*dx), 1)
    kurant(2) = dt*maxval(maxval(maxval(abs(c(1:Nx1, 1:Ny1, 1:Nz1, 2, 1)), 1), 2)/(a*dy), 1)
    kurant(3) = dt*maxval(maxval(maxval(abs(c(1:Nx1, 1:Ny1, 1:Nz1, 3, 1)), 1), 1)/(-dlev), 1)

    write(*, '(7a7)'), 'dt(s)', 'dlam', 'dphi', 'dz(km)', 'cx', 'cy', 'cz'
    write(*, '(7f7.3)'), dt, dx*180/pi, dy*180/pi, maxval(-dlev)/1.e3, kurant

    write(*, '(3a10)'), 'lambda0', 'phi0', 'z0(km)'
    write(*, '(3f10.1)'), lambda_c*180/pi, phi_c*180/pi, z_c/1.e3

    ! print*, dt, dx, dy, dz

    ! write output to file

    nrec = 1
    uout(1:Nx1, 1:Ny1, 1:Nz1) = u(1:Nx1, 1:Ny1, 1:Nz1, 1)
    write(20, rec = nrec)uout

    ! main loop on time
    do n = 1, Nt

        call wind_ini(Nx, Ny, Nz, Nx1, Ny1, Nz1, v, c(:, :, :, :, 2), lambda, lambda_m, phi, z, z_m, u0, w0, alpha, Omega, eps, a, dt, pi, n)
        call forcing_calc(Nx, Ny, Nz, Nx1, Ny1, Nz1, c, u, fu, dlam, dphi, dlev, cphi, cphi_m, a, dt)
        call cab_adv(Nx, Ny, Nz, Nx1, Ny1, Nz1, c, u, fu, ux, uy, uz, dlam, dphi, dlev, phi_m, cphi, cphi_m, lambda, lambda_m, z_m, a, pi, dt)

        u(1:Nx1, 1:Ny1, 1:Nz1, 1) = u(1:Nx1, 1:Ny1, 1:Nz1, 2)
        ux(1:Nx1, 1:Ny1, 1:Nz1, 1) = ux(1:Nx1, 1:Ny1, 1:Nz1, 2)
        uy(1:Nx1, 1:Ny, 1:Nz1, 1) = uy(1:Nx1, 1:Ny, 1:Nz1, 2)
        uz(1:Nx1, 1:Ny1, 1:Nz, 1) = uz(1:Nx1, 1:Ny1, 1:Nz, 2)
        c(1:Nx1, 1:Ny, 1:Nz, 1:3, 1) = c(1:Nx1, 1:Ny, 1:Nz, 1:3, 2)

        ! write output to file

        if (mod(n, Nt/Nout) == 0) then
            nrec = nrec+1
            uout(1:Nx1, 1:Ny1, 1:Nz1) = u(1:Nx1, 1:Ny1, 1:Nz1, 1)
            write(20, rec = nrec)uout
            print*, 'Nout = ', n*Nout/Nt
        end if

    end do

    ! estimation norms of error field

    max_val = maxval(u(1:Nx1, 1:Ny1, 1:Nz1, 1))/maxval(ue(1:Nx1, 1:Ny1, 1:Nz1))
    min_val = minval(u(1:Nx1, 1:Ny1, 1:Nz1, 1))/maxval(ue(1:Nx1, 1:Ny1, 1:Nz1))
    call mass_calc(Nx, Ny, Nz, Nx1, Ny1, Nz1, ue, dlam, dphi, dlev, cphi, mass(1))
    call mass_calc(Nx, Ny, Nz, Nx1, Ny1, Nz1, u(:, :, :, 1), dlam, dphi, dlev, cphi, mass(2))
    call mass_calc(Nx, Ny, Nz, Nx1, Ny1, Nz1, ue**2, dlam, dphi, dlev, cphi, mom2(1))
    call mass_calc(Nx, Ny, Nz, Nx1, Ny1, Nz1, u(:, :, :, 1)**2, dlam, dphi, dlev, cphi, mom2(2))

    mass_rel = mass(2)/mass(1)
    mom2_rel = mom2(2)/mom2(1)

    res(1) = mass_rel
    res(2) = mom2_rel
    res(3) = min_val
    res(4) = max_val


    write(*, '(5a12)'), 'Mass_rel', 'Mom2_rel', 'Min_rel', 'Max_rel', 'CPU(min)'
    write(*, '(5f12.5)'), res

    ! write output to file

    write(21, '(5f12.8)')res

end subroutine step_of_transfer

end module transfer_mod