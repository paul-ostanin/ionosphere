program mat

use tridi_matrix
use vector

implicit none

type (tridiagonal_matrix) S_y, S_z
type (vect) F_ub, rhs_z, rhs_y, z, h, hmid, nO, nO2, nN2, gradnO, k, p, pcurr, m, n_old, n_new, delta, D, D_node, gradD, u, gradu, Tn, Ti, Te, Tr, Tp, gradTp, gradTn, gradTr, gradgradTp, n_day, tau0, model_sol(1801), n_new_z(1801), p_mod(1801), n_old_z(1801), n_new_y(701), n_old_y(701), n_new_z_1(1801), n_old_z_1(1801), n_new_y_1(701), n_old_y_1(701), error(1801)
integer i, j, t, Te0, Tn0, Ti0, day, nonlinear_scheme_type, profile_output, diurnal_on, convergence_test, Nphi, Nz, pk_switch, mixed_z_switch, mixed_y_switch, transf_yz_switch, transf_y_switch, upper_bound_type, monotonizator
real (8) u_analytical, D_analytical, df_dz, df_dphi, ddf_dzdz, ddf_dphidphi, ddf_dzdphi, f, paramet, tau, tau_1, Hmax, h0, F_z, delta_norm, eps, tgdelta, sindelta, cosdelta, dphi, phi, integral, n_max, h_max, coschi, pi, omega, sigma_O2, sigma_N2, sigma_O, sI, cI, R, u_phi, u_phi_1, u_phi_2, u_phi_3, u_phi_4, u_z, u_z_mh, u_z_ph, u_z_m1, u_z_p1, x, A, B, dA_dphi, dB_dphi, u_phi_ph, u_phi_mh, Ndays, Niter, delta_err, sIph, cIph, sImh, cImh, l1_norm, l1_norm_err

!opening files for writing the output

!    open(unit=100, name='an.txt')
!    open(unit=110, name='num.txt')

!    open(unit=101, name='df_dphi_gnp.txt')
!    open(unit=102, name='ddf_dzdz_gnp.txt')
!    open(unit=103, name='ddf_dphidphi_gnp.txt')
!    open(unit=104, name='ddf_dzdphi_gnp.txt')
!    open(unit=21, name='rhs_components.txt')
!    open(unit=22, name='approx_components.txt')

paramet = 1

pi = 3.141592653589793238462643


 
!number of nodes in phi
Nphi = 180
!number of nodes in z
Nz = 81
!maximum altitude in km
Hmax = 500
!latitude
dphi = pi / Nphi
!angle velocity of the Earth 
omega = 2*pi/24/60/60
!magnetic inclination sin I
sI = 1
!Earth radius
R = 637100000
!number of calculation days
Ndays = 20
Niter = 800
!upper boundary electron flow
F_z = 0

!Time step (in seconds) 5 min
tau = 1

!photochemistry switcher
pk_switch = 1
!mixed derivative u_phi switcher
mixed_z_switch = 1
!mixed derivative u_z switcher
mixed_y_switch = 1
!transfer d/dz(u n) and d/dphi(B(phi) n) switcher (multiplyer of u)
transf_yz_switch = 1
!transfer d/dphi(B(phi) n) switcher
transf_y_switch = 1
!switcher for the upper boundary approximation of the first splitting step
upper_bound_type = 4
!switcher for the monotonizator after both steps - all less or equal than zero elements are made 1
monotonizator = 1
diurnal_on = 0
convergence_test = 0
profile_output = 0



!Initialization block
    !Vector of altitudes. Step h_i = z(i) - z(i - 1). Counting from 100 km to 500 km. z.d(i) is in metres.
    call z.init(Nz)

    !Space step (in cm) 5 km
    h0 = (Hmax - 100) * 1E+5 / (z.n - 1)
    do i = 1, z.n
        z.d(i) = 100E+3 + h0 * (i-1)/100
    end do

    !Vector of middles of [z(i), z(i+1)].  m.d(i) = z_{i+1/2} (in metres)
    call m.init(z.n)
    do i = 1, z.n - 1
        m.d(i) = 0.5 * (z.d(i) + z.d(i + 1))
    end do
    m.d(z.n) = z.d(z.n) + (z.d(z.n) - m.d(z.n - 1))

    call h.init(z.n - 1)
    do i = 1, z.n - 1
        h.d(i) = 100 * (z.d(i + 1) - z.d(i))
    end do

    call hmid.init(z.n - 1)
    do i = 1, z.n - 2
         hmid.d(i) = 100 * (m.d(i + 1) - m.d(i))
    end do
    hmid.d(z.n - 1) = hmid.d(z.n - 2)

    call Ti.init(z.n+1)
    call Tn.init(z.n+1)
    call Te.init(z.n+1)
    call Tr.init(z.n+1)
    call Tp.init(z.n+1)
    Ti0 = 950 
    Tn0 = 800
    Te0 = 2200
    do i = 1, z.n
        Ti.d(i) = Ti0 - (Ti0 - 200) * exp(-9.8 / (287 * Ti0) * (z.d(i) - 100000))
        Tn.d(i) = Tn0 - (Tn0 - 200) * exp(-9.8 / (287 * Tn0) * (z.d(i) - 100000))
        Te.d(i) = Te0 - (Te0 - 200) * exp(-9.8 / (287 * Te0) * (z.d(i) - 100000))
        Tr.d(i) = 0.5 * (Tn.d(i) + Ti.d(i))
        Tp.d(i) = 0.5 * (Te.d(i) + Ti.d(i))
    end do
        Ti.d(z.n+1) = Ti0 - (Ti0 - 200) * exp(-9.8 / (287 * Ti0) * (z.d(z.n) + h0/100 - 100000))
        Tn.d(z.n+1) = Tn0 - (Tn0 - 200) * exp(-9.8 / (287 * Tn0) * (z.d(z.n) + h0/100 - 100000))
        Te.d(z.n+1) = Te0 - (Te0 - 200) * exp(-9.8 / (287 * Te0) * (z.d(z.n) + h0/100 - 100000))
        Tr.d(z.n+1) = 0.5 * (Tn.d(z.n+1) + Ti.d(z.n+1))
        Tp.d(z.n+1) = 0.5 * (Te.d(z.n+1) + Ti.d(z.n+1))

    !dT/dz = s(T_\infty - T_0)exp(-s(z-z0)) = (T_\infty - T)*s - for Te, Tn, Ti. This derivative is in [K/cm].
    call gradTp.init(z.n+1)
    call gradgradTp.init(z.n+1)
    call gradTr.init(z.n+1)
    call gradTn.init(z.n+1)
    do i = 1, z.n+1
        gradTp.d(i)     =  0.5 * ((Te0 - Te.d(i)) * 9.8 / (287 * Te0) + (Ti0 - Ti.d(i)) * 9.8 / (287 * Ti0) ) * 1E-2
        gradgradTp.d(i) = -0.5 * ((Te0 - Te.d(i)) * (9.8 / (287 * Te0))**2 + (Ti0 - Ti.d(i)) * (9.8 / (287 * Ti0))**2 ) * 1E-4
        gradTr.d(i)     =  0.5 * ((Tn0 - Tn.d(i)) * 9.8 / (287 * Tn0) + (Ti0 - Ti.d(i)) * 9.8 / (287 * Ti0) ) * 1E-2       
        gradTn.d(i)     =         (Tn0 - Tn.d(i)) * 9.8 / (287 * Tn0) * 1E-2
    end do

    call nO.init(z.n+1)
    call nO2.init(z.n+1)
    call nN2.init(z.n+1)
    do i = 1, z.n
        nO.d(i)  = 2.8E+10 * exp(-9.8 * 16E-3 / (8.31 * Tn.d(i)) * (z.d(i) - 140000))
        nO2.d(i) = 5.6E+9  * exp(-9.8 * 32E-3 / (8.31 * Tn.d(i)) * (z.d(i) - 140000))
        nN2.d(i) = 5.2E+10 * exp(-9.8 * 28E-3 / (8.31 * Tn.d(i)) * (z.d(i) - 140000))
    end do
        nO.d(z.n+1)  = 2.8E+10 * exp(-9.8 * 16E-3 / (8.31 * Tn.d(z.n+1)) * (z.d(z.n) + h0/100 - 140000))
        nO2.d(z.n+1) = 5.6E+9  * exp(-9.8 * 32E-3 / (8.31 * Tn.d(z.n+1)) * (z.d(z.n) + h0/100 - 140000))
        nN2.d(z.n+1) = 5.2E+10 * exp(-9.8 * 28E-3 / (8.31 * Tn.d(z.n+1)) * (z.d(z.n) + h0/100 - 140000))

    ! nO(z) = nO(z) * gm/RTn ((z-140000)*1/Tn*dTn/dz - 1/100), 1/100 = 1/100 m/cm, derivative is in cm^-4
    call gradnO.init(z.n+1)
    do i = 1, z.n
        gradnO.d(i) = nO.d(i) * (9.8*16E-3/(8.31*Tn.d(i))) * ((z.d(i)-140000)/Tn.d(i)*gradTn.d(i) - 1E-2)
    end do
        gradnO.d(z.n+1) = nO.d(z.n+1) * (9.8*16E-3/(8.31*Tn.d(z.n+1))) * ((z.d(z.n) + h0/100 - 140000)/Tn.d(z.n+1)*gradTn.d(z.n+1) - 1E-2)

    call p.init(z.n+1)
    call k.init(z.n+1)
    do i = 1, z.n+1
            p.d(i) = ( 4E-7 * nO.d(i) ) * pk_switch
            k.d(i) = ( 1.2E-12 * nN2.d(i) + 2.1E-11 * nO2.d(i) ) * pk_switch
    end do

    !Diffusion coefficients vector. D.d(i) = D_{i+1/2}
    call D.init(z.n)
    call D_node.init(z.n+1)
    do i = 1, z.n
        D.d(i) = D_analytical(z.d(i)+h0/200)
        D_node.d(i) = D_analytical(z.d(i))
    end do
        D_node.d(z.n+1) = D_analytical(z.d(z.n) + h0/100)

    call gradD.init(z.n+1)
    do i = 1, z.n+1
        gradD.d(i) = D_node.d(i) * ( gradTp.d(i)/Tp.d(i) - (gradnO.d(i))/nO.d(i) - gradTr.d(i)/(2*Tr.d(i)) )
    end do

    !Effective velocity vector
    !u.d(i) = u_i = D_i * (dTp/dz + mg/2k)/Tp, mg/2k ~ 5.6*10^{-5} [K/cm]
    call u.init(z.n+1)
    do i = 1, z.n+1    
        u.d(i) = ( 3E+17 / (nO.d(i) * sqrt(Tr.d(i))) * (56E-6 + gradTp.d(i)) ) * transf_yz_switch
    end do

    call gradu.init(z.n+1)
    do i = 1, z.n+1
        gradu.d(i) = 3E+17 * ( gradgradTp.d(i)/(nO.d(i) * sqrt(Tr.d(i))) - &
                               (56E-6 + gradTp.d(i)) * gradnO.d(i)/ (nO.d(i)**2 * sqrt(Tr.d(i))) - &
                               (56E-6 + gradTp.d(i)) * gradTr.d(i)/2/(nO.d(i)**2 * sqrt(Tr.d(i))) )
    end do


    !System matrix (1-st step)
    call S_z.init(z.n)

    !System matrix (2-nd step)
    call S_y.init(Nphi)

    call n_old.init(z.n)
    call n_new.init(z.n)
    call delta.init(z.n)
    call n_old.gen() !n_old = (1, 1, ..., 1)
    call n_day.init(z.n)
    call rhs_z.init(z.n) !generating RHS for S_z * n_new = rhs_z
    call rhs_y.init(Nphi)
    call F_ub.init(Nphi)
    call pcurr.init(z.n)

    do j = 1, Nphi
        call n_new_z(j).init(z.n)
        call n_old_z(j).init(z.n)
        call model_sol(j).init(z.n)
        call n_new_z_1(j).init(z.n)
        call n_old_z_1(j).init(z.n)
        call error(j).init(z.n)
        call n_old_z(j).gen()
        call n_old_z(j).gen()
        call p_mod(j).init(z.n)
    end do

    do i = 1, z.n
        call n_new_y(i).init(Nphi)
        call n_old_y(i).init(Nphi)
        call n_new_y_1(i).init(Nphi)
        call n_old_y_1(i).init(Nphi)
    end do

    do j = 1, Nphi
        phi = (j-0.5)*dphi-pi/2
        sI   = sin(atan(2*tan(-pi/2+(j-0.5)*dphi)))
        cI   = cos(atan(2*tan(-pi/2+(j-0.5)*dphi)))
        !F_ub.d(j) = D_node.d(z.n)*(sI**2)*df_dz(z.d(z.n), phi) + u.d(z.n)*(sI**2)*f(z.d(z.n), phi) - 1/R*sI*cI * df_dphi(z.d(z.n), phi)
        F_ub.d(j) =  D.d(z.n-1)*(sI**2)*df_dz(z.d(z.n)-h0/200, phi) + u_analytical(z.d(z.n)-h0/200)*(sI**2)*f(z.d(z.n)-h0/200, phi) - 1/R*sI*cI*D.d(z.n-1)*df_dphi(z.d(z.n)-h0/200, phi)*mixed_z_switch
        do i = 1, z.n
            
            p_mod(j).d(i) = k.d(i)*f(z.d(i), phi) - ( &
                            sI**2 * (gradD.d(i)*df_dz(z.d(i), phi) + ddf_dzdz(z.d(i), phi)*D_node.d(i)) + &
                            sI**2 * (gradu.d(i)*f(z.d(i), phi) + u.d(i)*df_dz(z.d(i), phi)) - &
                            mixed_z_switch*1/R*sI*cI * (gradD.d(i)*df_dphi(z.d(i), phi) + D_node.d(i)*ddf_dzdphi(z.d(i), phi)) + &
                            D_node.d(i)/(R**2)/cos(phi) * (dA_dphi(phi)*df_dphi(z.d(i), phi) + A(phi)*ddf_dphidphi(z.d(i), phi)) - &
                            mixed_y_switch*D_node.d(i)/R/cos(phi)/2 * (dB_dphi(phi)*df_dz(z.d(i), phi) + B(phi)*ddf_dzdphi(z.d(i), phi)) - &
                            u.d(i)/R/cos(phi)/2 * (dB_dphi(phi)*f(z.d(i), phi) + B(phi)*df_dphi(z.d(i), phi)) &
                            )
            !if(j == 20) then
            !    write(21, *) sI**2 * (gradD.d(i)*df_dz(z.d(i), phi) + ddf_dzdz(z.d(i), phi)*D_node.d(i))
            !end if
        end do
    end do


!stationary model
do t = 0, Ndays*86400/tau
    if(mod(t, 100) .eq. 0) then
        print *, t
    end if
    do j = 1, Nphi
    ! first splitting step

        !sinus and cosinus of magnetic inclination angle I
        sI   = sin(atan(2*tan(-pi/2+(j-0.5)*dphi)))
        cI   = cos(atan(2*tan(-pi/2+(j-0.5)*dphi)))
        sIph = sin(atan(2*tan(-pi/2+(j-0.5)*dphi+dphi/2)))
        cIph = cos(atan(2*tan(-pi/2+(j-0.5)*dphi+dphi/2)))
        sImh = sin(atan(2*tan(-pi/2+(j-0.5)*dphi-dphi/2)))
        cImh = cos(atan(2*tan(-pi/2+(j-0.5)*dphi-dphi/2)))

        !lower boundary condition: n_1 = P_1/k_1
        S_z.d(1, 2) = 1
        rhs_z.d(1) = p_mod(j).d(1)/k.d(1)


        if (j .eq. 1 .or. j .eq. Nphi) then
        !upper boundary type 0: no mixed derivative in the upper boundary condition
            if(sI .ge. 0) then
                S_z.d(z.n, 1) = (-D.d(z.n-1)*tau/(h0**2) + 0.5*u.d(z.n-1)*tau/h0) * sI**2 - 0.5*tau/(R*h0*dphi)*D.d(z.n-1)*sI*cI
                S_z.d(z.n, 2) = ( D.d(z.n-1)*tau/(h0**2) + 0.5*u.d( z.n )*tau/h0) * sI**2 !+ 0.5*tau/(R*h0*dphi)*D.d(z.n-1)*sI*cI
                rhs_z.d(z.n)  = tau/h0*F_ub.d(j) + 0.5*tau/(R*h0*dphi)*D.d(z.n-1)*sI*cI * (0*n_old_z(j+1).d(z.n) - n_old_z(j-1).d(z.n-1))
            else
                S_z.d(z.n, 1) = (-D.d(z.n-1)*tau/(h0**2) + 0.5*u.d(z.n-1)*tau/h0) * sI**2 + 0.5*tau/(R*h0*dphi)*D.d(z.n-1)*sI*cI
                S_z.d(z.n, 2) = ( D.d(z.n-1)*tau/(h0**2) + 0.5*u.d( z.n )*tau/h0) * sI**2 !- 0.5*tau/(R*h0*dphi)*D.d(z.n-1)*sI*cI
                rhs_z.d(z.n)  = tau/h0*F_ub.d(j) + 0.5*tau/(R*h0*dphi)*D.d(z.n-1)*sI*cI * (n_old_z(j+1).d(z.n-1) - 0*n_old_z(j-1).d(z.n))
            end if    
        else 
            if(sI .ge. 0) then
                S_z.d(z.n, 1) = (-D.d(z.n-1)*tau/(h0**2) + 0.5*u.d(z.n-1)*tau/h0) * sI**2 - 0.5*tau/(R*h0*dphi)*D.d(z.n-1)*sI*cI
                S_z.d(z.n, 2) = ( D.d(z.n-1)*tau/(h0**2) + 0.5*u.d( z.n )*tau/h0) * sI**2 + 0.5*tau/(R*h0*dphi)*D.d(z.n-1)*sI*cI
                rhs_z.d(z.n)  = tau/h0*F_ub.d(j) + 0.5*tau/(R*h0*dphi)*D.d(z.n-1)*sI*cI * (n_old_z(j+1).d(z.n) - n_old_z(j-1).d(z.n-1))
            else
                S_z.d(z.n, 1) = (-D.d(z.n-1)*tau/(h0**2) + 0.5*u.d(z.n-1)*tau/h0) * sI**2 + 0.5*tau/(R*h0*dphi)*D.d(z.n-1)*sI*cI
                S_z.d(z.n, 2) = ( D.d(z.n-1)*tau/(h0**2) + 0.5*u.d( z.n )*tau/h0) * sI**2 - 0.5*tau/(R*h0*dphi)*D.d(z.n-1)*sI*cI
                rhs_z.d(z.n)  = tau/h0*F_ub.d(j) + 0.5*tau/(R*h0*dphi)*D.d(z.n-1)*sI*cI * (n_old_z(j+1).d(z.n-1) - n_old_z(j-1).d(z.n))
            end if
        end if

        do i = 2, z.n - 1
            if (j .eq. 1 .or. j .eq. Nphi) then
                if(sI .ge. 0) then
                    S_z.d(i, 1) =  (-D.d(i-1)*tau/(h0**2) + u.d(i-1)*tau/(2*h0)) * sI**2 - &
                        0.5*D_node.d(i-1)*tau*sImh*cImh/(R*h0*dphi)
                    S_z.d(i, 2) = 1 + k.d(i)*tau + (D.d(i-1) + D.d(i)) * tau/(h0**2) * sI**2 + &
                        0.5*D_node.d(i)*tau*sImh*cImh/(R*h0*dphi)
                    S_z.d(i, 3) =  (-D.d( i )*tau/(h0**2) - u.d(i+1)*tau/(2*h0)) * sI**2

                    rhs_z.d(i) = n_old_z(j).d(i) + tau*p_mod(j).d(i) - &
                        0.5*tau/(R*h0*dphi)*(-D_node.d( i ) * sImh*cImh * n_old_z(j-1).d( i ) + &
                                              D_node.d(i-1) * sImh*cImh * n_old_z(j-1).d(i-1) )
                else
                    S_z.d(i, 1) =  (-D.d(i-1)*tau/(h0**2) + u.d(i-1)*tau/(2*h0)) * sI**2 + &
                        0.5*D_node.d(i-1)*tau*sIph*cIph/(R*h0*dphi)
                    S_z.d(i, 2) = 1 + k.d(i)*tau + (D.d(i-1) + D.d(i))*tau/(h0**2) * sI**2 - &
                        0.5*D_node.d(i)*tau*sIph*cIph/(R*h0*dphi)
                    S_z.d(i, 3) =  (-D.d( i )*tau/(h0**2) - u.d(i+1)*tau/(2*h0)) * sI**2

                    rhs_z.d(i) = n_old_z(j).d(i) + tau*p_mod(j).d(i) - &
                        0.5*tau/(R*h0*dphi)*( D_node.d( i ) * sIph*cIph * n_old_z(j+1).d( i ) - &
                                              D_node.d(i-1) * sIph*cIph * n_old_z(j+1).d(i-1) )
                end if

            else
                if(sI .ge. 0) then
                    S_z.d(i, 1) =  (-D.d(i-1)*tau/(h0**2) + u.d(i-1)*tau/(2*h0)) * sI**2 - &
                        0.5*D_node.d(i-1)*tau*sImh*cImh/(R*h0*dphi)
                    S_z.d(i, 2) = 1 + k.d(i)*tau + (D.d(i-1) + D.d(i)) * tau/(h0**2) * sI**2 + &
                        0.5*D_node.d(i)*tau*(sIph*cIph+sImh*cImh)/(R*h0*dphi)
                    S_z.d(i, 3) =  (-D.d( i )*tau/(h0**2) - u.d(i+1)*tau/(2*h0)) * sI**2 - &
                        0.5*D_node.d(i+1)*tau*sIph*cIph/(R*h0*dphi)

                    rhs_z.d(i) = n_old_z(j).d(i) + tau*p_mod(j).d(i) - &
                        0.5*tau/(R*h0*dphi)*(-D_node.d( i ) * sImh*cImh * n_old_z(j-1).d( i ) - &
                                              D_node.d( i ) * sIph*cIph * n_old_z(j+1).d( i ) + &
                                              D_node.d(i-1) * sImh*cImh * n_old_z(j-1).d(i-1) + &
                                              D_node.d(i+1) * sIph*cIph * n_old_z(j+1).d(i+1) )
                else
                    S_z.d(i, 1) =  (-D.d(i-1)*tau/(h0**2) + u.d(i-1)*tau/(2*h0)) * sI**2 + &
                        0.5*D_node.d(i-1)*tau*sIph*cIph/(R*h0*dphi)
                    S_z.d(i, 2) = 1 + k.d(i)*tau + (D.d(i-1) + D.d(i))*tau/(h0**2) * sI**2 - &
                        0.5*D_node.d(i)*tau*(sImh*cImh+sIph*cIph)/(R*h0*dphi)
                    S_z.d(i, 3) =  (-D.d( i )*tau/(h0**2) - u.d(i+1)*tau/(2*h0)) * sI**2 + &
                        0.5*D_node.d(i+1)*tau*sImh*cImh/(R*h0*dphi)

                    rhs_z.d(i) = n_old_z(j).d(i) + tau*p_mod(j).d(i) - &
                        0.5*tau/(R*h0*dphi)*( D_node.d( i ) * sImh*cImh * n_old_z(j-1).d( i ) + &
                                              D_node.d( i ) * sIph*cIph * n_old_z(j+1).d( i ) - &
                                              D_node.d(i+1) * sImh*cImh * n_old_z(j-1).d(i+1) - &
                                              D_node.d(i-1) * sIph*cIph * n_old_z(j+1).d(i-1) )
                end if
            end if
        end do

        n_new_z(j) = tridiagonal_matrix_algorithm(S_z, rhs_z)

    end do

    !sending 1st step result to the 2nd step as the initial condition
    do i = 1, z.n
        do j = 1, Nphi
            if(n_new_z(j).d(i) .ge. 0 .or. monotonizator .eq. 0) then
                n_old_y(i).d(j) = n_new_z(j).d(i)
            else
                n_old_y(i).d(j) = 1
            end if
        end do
    end do

    !boundary conditions at z = 100 and z = 500
    do j = 1, Nphi
        n_new_y(1).d(j) = n_new_z(j).d(1)
        !n_new_y(z.n).d(j) = n_new_z(j).d(z.n)
    end do

    !simulating block for the 2nd step
        ! do i = 2, z.n-1
        ! do j = 1, Nphi
        !    n_new_y(i).d(j) = n_old_y(i).d(j)
        ! end do
        ! end do

    do i = 2, z.n
    !second splitting step
        if (i .ne. z.n) then
            do j = 1, Nphi    
                phi = (j-0.5)*dphi-pi/2
                if (j .ne. 1 .and. j .ne. Nphi) then
                    if(B(phi) .le. 0) then
                        S_y.d(j, 1) =  (-D_node.d(i)/R * A(phi-dphi/2)/(dphi**2) + (-u.d(i)/2)*B(phi-dphi)/(2*dphi) + &
                            0.5*mixed_y_switch*D.d( i )/2*B(phi-dphi)/h0/dphi)*tau/(R)/cos(phi)
                        S_y.d(j, 2) = 1+(D_node.d(i)/R *(A(phi-dphi/2) + A(phi+dphi/2))/(dphi**2) - &
                            0.5*mixed_y_switch*(D.d(i)+D.d(i-1))/2*B(phi)/h0/dphi) * tau/(R)/cos(phi)
                        S_y.d(j, 3) =  (-D_node.d(i)/R * A(phi+dphi/2)/(dphi**2) - (-u.d(i)/2)*B(phi+dphi)/(2*dphi) + &
                            0.5*mixed_y_switch*D.d(i-1)/2*B(phi+dphi)/h0/dphi)*tau/(R)/cos(phi)
                        rhs_y.d(j) = n_old_y(i).d(j)-&
                            0.5*mixed_y_switch*tau/(2*R*cos(phi)*h0*dphi)*( &
                                B(   phi  ) * D.d(i-1) * n_old_y(i-1).d( j ) + &
                                B(   phi  ) * D.d( i ) * n_old_y(i+1).d( j ) - &
                                B(phi+dphi) * D.d(i-1) * n_old_y(i-1).d(j+1) - &
                                B(phi-dphi) * D.d( i ) * n_old_y(i+1).d(j-1) )
                    else
                        S_y.d(j, 1) =  (-D_node.d(i)/R * A(phi-dphi/2)/(dphi**2) + (-u.d(i)/2)*B(phi-dphi)/(2*dphi)- &
                            0.5*mixed_y_switch*D.d(i-1)/2*B(phi-dphi)/h0/dphi)*tau/(R)/cos(phi)
                        S_y.d(j, 2) = 1+(D_node.d(i)/R *(A(phi-dphi/2) + A(phi+dphi/2))/(dphi**2) + &
                            0.5*mixed_y_switch*(D.d(i)+D.d(i-1))/2*B(phi)/h0/dphi) * tau/(R)/cos(phi)
                        S_y.d(j, 3) =  (-D_node.d(i)/R * A(phi+dphi/2)/(dphi**2) - (-u.d(i)/2)*B(phi+dphi)/(2*dphi)- &
                            0.5*mixed_y_switch*D.d( i )/2*B(phi+dphi)/h0/dphi)*tau/(R)/cos(phi)
                        rhs_y.d(j) = n_old_y(i).d(j)-&
                            0.5*mixed_y_switch*tau/(2*R*cos(phi)*h0*dphi)*(  - &
                                B(   phi  ) * D.d(i-1) * n_old_y(i-1).d( j ) - &
                                B(   phi  ) * D.d( i ) * n_old_y(i+1).d( j ) + &
                                B(phi-dphi) * D.d(i-1) * n_old_y(i-1).d(j-1) + &
                                B(phi+dphi) * D.d( i ) * n_old_y(i+1).d(j+1) )
                    end if    
                else
                    if(B(phi) .ge. 0) then
                        S_y.d(j, 1) =  (-D_node.d(i)/R * A(phi-dphi/2)/(dphi**2) - u.d(i)/2*B(phi-dphi)/(2*dphi)- &
                            0.5*mixed_y_switch*D.d(i-1)/2*B(phi-dphi)/h0/dphi)*tau/(R)/cos(phi)
                        S_y.d(j, 2) = 1+(D_node.d(i)/R * A(phi-dphi/2)/(dphi**2) + u.d(i)/2*B(phi+dphi)/(2*dphi) + &
                            0.5*mixed_y_switch*D.d(i-1)/2*B(phi)/h0/dphi) * tau/(R)/cos(phi)
                        S_y.d(j, 3) = 0
                        rhs_y.d(j) = n_old_y(i).d(j)-&
                            0.5*mixed_y_switch*tau/(2*R*cos(phi)*h0*dphi)*(  - &
                                B(   phi  ) * D.d(i-1) * n_old_y(i-1).d( j ) + &
                                B(phi-dphi) * D.d(i-1) * n_old_y(i-1).d(j-1) )
                    else
                        S_y.d(j, 1) = 0
                        S_y.d(j, 2) = 1+(D_node.d(i)/R * A(phi+dphi/2)/(dphi**2) - u.d(i)/2*B(phi-dphi)/(2*dphi) - &
                            0.5*mixed_y_switch*D.d(i-1)/2*B(phi)/h0/dphi) * tau/(R)/cos(phi)
                        S_y.d(j, 3) =  (-D_node.d(i)/R * A(phi+dphi/2)/(dphi**2) + u.d(i)/2*B(phi+dphi)/(2*dphi) + &
                            0.5*mixed_y_switch*D.d(i-1)/2*B(phi+dphi)/h0/dphi)*tau/(R)/cos(phi)
                        rhs_y.d(j) = n_old_y(i).d(j)-&
                            0.5*mixed_y_switch*tau/(2*R*cos(phi)*h0*dphi)*( &
                                B(   phi  ) * D.d(i-1) * n_old_y(i-1).d( j ) - &
                                B(phi+dphi) * D.d(i-1) * n_old_y(i-1).d(j+1) )
                    end if
                end if
            end do
        else
            do j = 2, Nphi-1    
                phi = (j-0.5)*dphi-pi/2
                if(B(phi) .le. 0) then
                    S_y.d(j, 1) = 0.5*tau/(R*h0*dphi)*D.d(z.n-1)*sI*cI
                    S_y.d(j, 2) = (D.d(z.n-1)*tau/(h0**2) + 0.5*u.d( z.n )*tau/h0) * sI**2 - 0.5*tau/(R*h0*dphi)*D.d(z.n-1)*sI*cI
                    S_y.d(j, 3) = 0
                    rhs_y.d(j)  = tau/h0*F_ub.d(j) + 0.5*tau/(R*h0*dphi)*D.d(z.n-1)*sI*cI * n_old_y(z.n-1).d(j+1) + &
                                    n_old_y(z.n-1).d(j) * ((D.d(z.n-1)*tau/(h0**2) - 0.5*u.d(z.n-1)*tau/h0) * sI**2 - 0.5*tau/(R*h0*dphi)*D.d(z.n-1)*sI*cI)
                else
                    S_y.d(j, 1) = 0
                    S_y.d(j, 2) = (D.d(z.n-1)*tau/(h0**2) + 0.5*u.d( z.n )*tau/h0) * sI**2 + 0.5*tau/(R*h0*dphi)*D.d(z.n-1)*sI*cI
                    S_y.d(j, 3) = -0.5*tau/(R*h0*dphi)*D.d(z.n-1)*sI*cI
                    rhs_y.d(j)  = tau/h0*F_ub.d(j) - 0.5*tau/(R*h0*dphi)*D.d(z.n-1)*sI*cI * n_old_y(z.n-1).d(j-1) + &
                                    n_old_y(z.n-1).d(j) * ((D.d(z.n-1)*tau/(h0**2) - 0.5*u.d(z.n-1)*tau/h0) * sI**2 + 0.5*tau/(R*h0*dphi)*D.d(z.n-1)*sI*cI)
                end if
            end do
            j = 1
                phi = (j-0.5)*dphi-pi/2
                S_y.d(j, 1) = 0
                S_y.d(j, 2) = (D.d(z.n-1)*tau/(h0**2) + 0.5*u.d( z.n )*tau/h0) * sI**2
                S_y.d(j, 3) = 0
                rhs_y.d(j)  = tau/h0*F_ub.d(j) + 0.5*tau/(R*h0*dphi)*D.d(z.n-1)*sI*cI * n_old_y(z.n-1).d(j+1) + &
                                n_old_y(z.n-1).d(j) * ((D.d(z.n-1)*tau/(h0**2) - 0.5*u.d(z.n-1)*tau/h0) * sI**2 - 0.5*tau/(R*h0*dphi)*D.d(z.n-1)*sI*cI)
            j = Nphi
                phi = (j-0.5)*dphi-pi/2
                S_y.d(j, 1) = 0
                S_y.d(j, 2) = (D.d(z.n-1)*tau/(h0**2) + 0.5*u.d( z.n )*tau/h0) * sI**2
                S_y.d(j, 3) = 0
                rhs_y.d(j)  = tau/h0*F_ub.d(j) - 0.5*tau/(R*h0*dphi)*D.d(z.n-1)*sI*cI * n_old_y(z.n-1).d(j-1) + &
                                n_old_y(z.n-1).d(j) * ((D.d(z.n-1)*tau/(h0**2) - 0.5*u.d(z.n-1)*tau/h0) * sI**2 + 0.5*tau/(R*h0*dphi)*D.d(z.n-1)*sI*cI)
        end if

        n_new_y(i) = tridiagonal_matrix_algorithm(S_y, rhs_y)
    end do


    !sending the result back to the 1-st step
    do i = 1, z.n
        do j = 1, Nphi
            if (n_new_y(i).d(j) .ge. 0 .or. monotonizator .eq. 0) then
                n_old_z(j).d(i) = n_new_y(i).d(j)
            else
                n_old_z(j).d(i) = 1
            end if
        end do
    end do

    !output block
    if(mod(t, 100) .eq. 0) then

        open(unit=11, name='res_split_model.txt')
        open(unit=12, name='res_split_gnp_model.txt')
        do j = 1, Nphi
            do i = 1, z.n
                write(12,*) (j-5E-1)*180/(Nphi)-90, 100+(Hmax - 100)/(z.n-1)*(i-1), n_old_z(j).d(i)
            end do
            write (12, *)
        end do
        do j = 1, Nphi
            call n_old_z(j).print(11)
        end do
        close(unit=11)
        close(unit=12)

        if(t .eq. 0) then
            open(unit=98, name='model_solution_gnp.txt')    
            open(unit=99, name='model_solution_matrix.txt')
            do j = 1, Nphi
                do i = 1, z.n
                    write(98,*) (j-5E-1)*180/(Nphi)-90, 100+(Hmax - 100)/(z.n-1)*(i-1), f(z.d(i), (j-5E-1)*pi/(Nphi)-pi/2)
                    model_sol(j).d(i) = f(z.d(i), (j-5E-1)*pi/(Nphi)-pi/2)
                end do
                write (98, *)
            end do
            do j = 1, Nphi
                call model_sol(j).print(99)
            end do
            close(unit=98)
            close(unit=99)
        end if

        open(unit=198,name='diff_res_sol_gnp.txt')
        open(unit=199,name='diff_res_sol.txt')
        do j = 1, Nphi
            do i = 1, z.n
                write(198,*) (j-5E-1)*180/(Nphi)-90, 100+(Hmax - 100)/(z.n-1)*(i-1), n_old_z(j).d(i)-model_sol(j).d(i)
            end do
            write (198, *)
        end do
        do j = 1, Nphi
            call n_old_z(j).print(199)
        end do
        close(unit=198)
        close(unit=199)

    end if
end do

!profiles output block
if(profile_output .eq. 1) then
    open(unit=88, name='step2_88_fz_-1e9.txt')
    open(unit=80, name='step2_80_fz_-1e9.txt')
    open(unit=70, name='step2_70_fz_-1e9.txt')
    open(unit=60, name='step2_60_fz_-1e9.txt')
    open(unit=50, name='step2_50_fz_-1e9.txt')
    open(unit=40, name='step2_40_fz_-1e9.txt')
    open(unit=30, name='step2_30_fz_-1e9.txt')
    open(unit=20, name='step2_20_fz_-1e9.txt')
    open(unit=10, name='step2_10_fz_-1e9.txt')
    open(unit=5,  name='step2_05_fz_-1e9.txt')
    open(unit=2,  name='step2_02_fz_-1e9.txt')
    open(unit=0,  name='step2_00_fz_-1e9.txt')
    do i = 1, z.n
        write(88,*) 100+400/(z.n-1)*(i-1), n_old_z(4 ).d(i)
        write(80,*) 100+400/(z.n-1)*(i-1), n_old_z(20).d(i)
        write(70,*) 100+400/(z.n-1)*(i-1), n_old_z(40).d(i)
        write(60,*) 100+400/(z.n-1)*(i-1), n_old_z(60).d(i)
        write(50,*) 100+400/(z.n-1)*(i-1), n_old_z(80).d(i)
        write(40,*) 100+400/(z.n-1)*(i-1), n_old_z(100).d(i)
        write(30,*) 100+400/(z.n-1)*(i-1), n_old_z(120).d(i)
        write(20,*) 100+400/(z.n-1)*(i-1), n_old_z(140).d(i)
        write(10,*) 100+400/(z.n-1)*(i-1), n_old_z(160).d(i)
        write(5, *) 100+400/(z.n-1)*(i-1), n_old_z(170).d(i)
        write(2, *) 100+400/(z.n-1)*(i-1), n_old_z(176).d(i)
        write(0, *) 100+400/(z.n-1)*(i-1), n_old_z(180).d(i)
    end do
end if


end

real*8 function A(phi)
    implicit none
    real*8 phi 
    A = cos(phi) / (1+4*(tan(phi)**2))
    return
end

real*8 function B(phi)
    implicit none
    real*8 phi
    B = 4*sin(phi) / (1+4*(tan(phi)**2))
    return
end

real*8 function dA_dphi(phi)
    implicit none
    real*8 phi 
    dA_dphi = -8*tan(phi) / cos(phi) / (1+4*(tan(phi)**2))**2 - sin(phi) / (1+4*(tan(phi)**2))
    return
end

real*8 function dB_dphi(phi)
    implicit none
    real*8 phi
    dB_dphi =  -32*tan(phi)**2 / cos(phi) / (1+4*(tan(phi)**2))**2 + 4*cos(phi) / (1+4*(tan(phi)**2))
    return
end

real*8 function f(z, phi)
    implicit none
    real*8 phi, z
    f =  3E+6 * ((z/1000)-100)/133 * exp((100-(z/1000))/133) * cos(phi/2)**2
    !f = (cos(z)*cos(phi))**2
    return
end

real*8 function df_dz(z, phi)
    implicit none
    real*8 phi, z
    df_dz =  3E+6 * exp((100-(z/1000))/133) * cos(phi/2)**2 * (1 - ((z/1000)-100)/133) / 133  * 1E-5
    !df_dz = -sin(2*z)*cos(phi)**2
    return
end

real*8 function df_dphi(z, phi)
    implicit none
    real*8 phi, z
    df_dphi =  3E+6 * ((z/1000)-100)/133 * exp((100-(z/1000))/133) * cos(phi/2)**2 * (-tan(phi/2))
    !df_dphi = -cos(z)**2 * sin(2*phi)
    return
end

real*8 function ddf_dzdz(z, phi)
    implicit none
    real*8 phi, z
    ddf_dzdz =  1E-10 * 3E+6 * cos(phi/2)**2 * &
            (-exp((100-(z/1000))/133) - &
            &exp((100-(z/1000))/133) * (1 - ((z/1000)-100)/133) )/133/133
    !ddf_dzdz = -2*cos(2*z)*cos(phi)**2
    return
end

real*8 function ddf_dphidphi(z, phi)
    implicit none
    real*8 phi, z
    ddf_dphidphi = 3E+6 * ((z/1000)-100)/133 * exp((100-(z/1000))/133) * (-cos(phi)/2)
    !ddf_dphidphi = -2*cos(z)**2 * sin(2*phi)
    return
end

real*8 function ddf_dzdphi(z, phi)
    implicit none
    real*8 phi, z
    ddf_dzdphi = 3E+6 * exp((100-(z/1000))/133) * cos(phi/2)**2 * (1 - ((z/1000)-100)/133) * (-tan(phi/2)) / 133 * 1E-5
    !ddf_dzdphi = sin(2*z)*sin(2*phi)
    return
end

real*8 function u_analytical(z)
    implicit none
    real*8 z, Ti0, Tn0, Te0, nO, Ti, Tn, Te, Tr, gradTp
    Ti0 = 950 
    Tn0 = 800
    Te0 = 2200
    Ti = Ti0 - (Ti0 - 200) * exp(-9.8 / (287 * Ti0) * (z - 100000))
    Tn = Tn0 - (Tn0 - 200) * exp(-9.8 / (287 * Tn0) * (z - 100000))
    Te = Te0 - (Te0 - 200) * exp(-9.8 / (287 * Te0) * (z - 100000))
    Tr = 0.5 * (Tn + Ti)
    nO = 2.8E+10 * exp(-9.8 * 16E-3 / (8.31 * Tn) * (z - 140000))
    gradTp = 0.5 * ((Te0 - Te) * 9.8 / (287 * Te0) + (Ti0 - Ti) * 9.8 / (287 * Ti0) ) * 1E-2
    u_analytical =  3E+17 / (nO * sqrt(Tr)) * (56E-6 + gradTp) 
    return
end

real*8 function D_analytical(z)
    implicit none
    real*8 z, Ti0, Tn0, Te0, nO, Ti, Tn, Te, Tr, Tp, gradTp
    Ti0 = 950 
    Tn0 = 800
    Te0 = 2200
    Ti = Ti0 - (Ti0 - 200) * exp(-9.8 / (287 * Ti0) * (z - 100000))
    Tn = Tn0 - (Tn0 - 200) * exp(-9.8 / (287 * Tn0) * (z - 100000))
    Te = Te0 - (Te0 - 200) * exp(-9.8 / (287 * Te0) * (z - 100000))
    Tr = 0.5 * (Tn + Ti)
    Tp = 0.5 * (Te + Ti)
    nO = 2.8E+10 * exp(-9.8 * 16E-3 / (8.31 * Tn) * (z - 140000))
    D_analytical =  3E+17 * Tp / (nO * sqrt(Tr))
    return
end