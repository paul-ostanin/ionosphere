program mat

use vector
use transfer_mod, only : init_transfer, step_of_transfer, uout
use functions
!implicit none

type (vect) F_ub, z, h, hmid, nO, nO2, nN2, k, p, pcurr, m, D, D_node, u, Tn, Ti, Te, Tr, Tp, n_day, tau0, p_mod(181), model_sol(1801)
integer printing_mode, s_i, s_j, i, j, t, Te0, Tn0, Ti0, day, profile_output, diurnal_on, convergence_test, pk_switch, mixed_z_switch, mixed_y_switch, transf_yz_switch, transf_y_switch, second_step_scheme_type, upper_bound_type, monotonizator, model_solution
real (8) lambda, diurnal_start_time, D_analytical, u_analytical
real (8) n_max, h_max, tau, Hmax, h0, F_z, eps, tgdelta, sindelta, cosdelta, dphi, phi, coschi, pi, omega, sigma_O2, sigma_N2, sigma_O, sI, cI, R, x, A, B, Ndays, Niter, sIph, cIph, sImh, cImh, l1_norm, l1_norm_err
type (vect) gradu, gradD, gradTp, gradTr, gradTi, gradTn, gradTe, gradgradTp, gradnO
real (8) n_mod, dn_dz, dn_dphi, ddn_dzdz, ddn_dzdphi, ddn_dphidphi, vphi, dvphi_dphi, vz, dvz_dz


integer, parameter :: max_size    = 10000000
integer, parameter :: max_nonzero = 100000000
integer, parameter :: maxwr = max_nonzero + 8 * max_size
integer, parameter :: maxwi = max_nonzero + 2 * max_size + 1
integer, parameter :: Nz = 80
integer, parameter :: Nphi = 90!180
integer, parameter :: Nlambda = 145

integer, allocatable:: ia(:), ja(:), iw(:)
double precision, allocatable:: arr(:), f(:, :), n(:, :), rw(:)
real, allocatable:: z_arr(:)

external matvec, prevec0
integer ITER, INFO, NUNIT, ierr, ipalu, ipjlu, ipju, ipiw
real*8 RESID, ddot

external ddot, dcopy

integer imatvec(1), iprevec(1), ipbcg, matrix_size, nonzero
real*8 resinit

integer counter_a, counter_ia, counter_rhs
real*8 C_norm, L2_norm

real*8 stencil(Nz, Nphi, -1:1, -1:1), operator(-1:1, -1:1), diffusion_transfer_z(-1:1, -1:1), diffusion_transfer_y(-1:1, -1:1), mixed_z(-1:1, -1:1), mixed_y(-1:1, -1:1)
real*8 rhs(Nz, Nphi, Nlambda), ans(Nz, Nphi, Nlambda), ans_prev(Nz, Nphi, Nlambda), ans_inc(Nz, Nphi, Nlambda), rhs_analytical(Nz, Nphi)

allocate(ia(max_size + 1), ja(max_nonzero), iw(maxwi), arr(max_nonzero), f((Nz+1)*(Nphi+1), Nlambda), n((Nz+1)*(Nphi+1), Nlambda), rw(maxwr))



pi = 3.141592653589793238462643

!nonlinear_scheme_type variable switches the u_phi-approximation. 
!nonlinear_scheme_type = 8
 
!maximum altitude in m
Hmax = 460000
!latitude
dphi = pi / Nphi
!angle velocity of the Earth 
omega = 2*pi/24/60/60
!magnetic inclination sin I
sI = 1
!Earth radius
R = 6.37e6
!number of calculation days
Ndays = 500
Niter = 800
!upper boundary electron flow
F_z = 0

!Time step (in seconds) 5 min
tau =  10

!Before diurnal_start_time - stationary daytime, after - diurnal evolution, depending on lambda
diurnal_start_time = 0

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
!switcher of schemes for the second step mixed derivative: 0 is nonlinear, 1 is 1st order, 2 is 2nd order
second_step_scheme_type = 2
!switcher for the upper boundary approximation of the first splitting step
upper_bound_type = 1
!switcher for the monotonizator after both steps - all less or equal than zero elements are made 1
monotonizator = 0
model_solution = 1
diurnal_on = 0
convergence_test = 1
profile_output = 0
printing_mode = 1

if(convergence_test .eq. 1) then
    open(unit=50, name='convergence_20.txt')
end if

!Initialization block
    !Vector of altitudes. Step h_i = z(i) - z(i - 1). Counting from 100 km to 500 km. z.d(i) is in metres.
    call z.init(Nz)

    !Space step (in m) 5 km
    h0 = (Hmax - 100000) / (z.n - 1)
    do i = 1, z.n
        z.d(i) = 100E+3 + h0 * (i-1)
    end do

    !Vector of middles of [z(i), z(i+1)].  m.d(i) = z_{i+1/2} (in metres)
    call m.init(z.n)
    do i = 1, z.n - 1
        m.d(i) = 0.5 * (z.d(i) + z.d(i + 1))
    end do
    m.d(z.n) = z.d(z.n) + (z.d(z.n) - m.d(z.n - 1))

    call h.init(z.n - 1)
    do i = 1, z.n - 1
        h.d(i) = (z.d(i + 1) - z.d(i))
    end do

    call hmid.init(z.n - 1)
    do i = 1, z.n - 2
         hmid.d(i) = (m.d(i + 1) - m.d(i))
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
        Ti.d(z.n+1) = Ti0 - (Ti0 - 200) * exp(-9.8 / (287 * Ti0) * (z.d(z.n) + h0 - 100000))
        Tn.d(z.n+1) = Tn0 - (Tn0 - 200) * exp(-9.8 / (287 * Tn0) * (z.d(z.n) + h0 - 100000))
        Te.d(z.n+1) = Te0 - (Te0 - 200) * exp(-9.8 / (287 * Te0) * (z.d(z.n) + h0 - 100000))
        Tr.d(z.n+1) = 0.5 * (Tn.d(z.n+1) + Ti.d(z.n+1))
        Tp.d(z.n+1) = 0.5 * (Te.d(z.n+1) + Ti.d(z.n+1))

    !dT/dz = s(T_\infty - T_0)exp(-s(z-z0)) = (T_\infty - T)*s - for Te, Tn, Ti. This derivative is in [K/m].
    call gradTp.init(z.n+1)
    call gradgradTp.init(z.n+1)
    call gradTr.init(z.n+1)
    call gradTn.init(z.n+1)
    do i = 1, z.n+1
        gradTp.d(i)     =  0.5 * ((Te0 - Te.d(i)) * 9.8 / (287 * Te0) + (Ti0 - Ti.d(i)) * 9.8 / (287 * Ti0) ) 
        gradgradTp.d(i) = -0.5 * ((Te0 - Te.d(i)) * (9.8 / (287 * Te0))**2 + (Ti0 - Ti.d(i)) * (9.8 / (287 * Ti0))**2 ) 
        gradTr.d(i)     =  0.5 * ((Tn0 - Tn.d(i)) * 9.8 / (287 * Tn0) + (Ti0 - Ti.d(i)) * 9.8 / (287 * Ti0) )       
        gradTn.d(i)     =         (Tn0 - Tn.d(i)) * 9.8 / (287 * Tn0)
    end do



    call nO.init(z.n+1)
    call nO2.init(z.n+1)
    call nN2.init(z.n+1)

    do i = 1, z.n
        !*1E+6 - from cm^{-3} to m^{-3}
        nO.d(i)  = 2.8E+10 * exp(-9.8 * 16E-3 / (8.31 * Tn.d(i)) * (z.d(i) - 140000)) * 1E+6
        nO2.d(i) = 5.6E+9  * exp(-9.8 * 32E-3 / (8.31 * Tn.d(i)) * (z.d(i) - 140000)) * 1E+6
        nN2.d(i) = 5.2E+10 * exp(-9.8 * 28E-3 / (8.31 * Tn.d(i)) * (z.d(i) - 140000)) * 1E+6
    end do

        nO.d(z.n+1)  = 2.8E+10 * exp(-9.8 * 16E-3 / (8.31 * Tn.d(z.n+1)) * (z.d(z.n) + h0 - 140000)) * 1E+6
        nO2.d(z.n+1) = 5.6E+9  * exp(-9.8 * 32E-3 / (8.31 * Tn.d(z.n+1)) * (z.d(z.n) + h0 - 140000)) * 1E+6
        nN2.d(z.n+1) = 5.2E+10 * exp(-9.8 * 28E-3 / (8.31 * Tn.d(z.n+1)) * (z.d(z.n) + h0 - 140000)) * 1E+6

    ! nO(z) = nO(z) * gm/RTn ((z-140000)*1/Tn*dTn/dz - 1),  derivative is in m^-4
    call gradnO.init(z.n+1)
    do i = 1, z.n
        gradnO.d(i) = nO.d(i) * (9.8*16E-3/(8.31*Tn.d(i))) * ((z.d(i)-140000)/Tn.d(i)*gradTn.d(i) - 1)
    end do
        gradnO.d(z.n+1) = nO.d(z.n+1) * (9.8*16E-3/(8.31*Tn.d(z.n+1))) * ((z.d(z.n) + h0 - 140000)/Tn.d(z.n+1)*gradTn.d(z.n+1) - 1)



    !P = P_0 exp(tau_0(z)*(1-sec chi))
    !tau_0 = sum_{N2, O2, O} sigma_i R_0*T_n/(M_i g) n_i(z)
    !sigma is in cm^2
            ! sigma_O2 = 2E-17 
            ! sigma_N2 = 15E-18
            ! sigma_O  = 1E-17
            ! call tau0.init(z.n+1)
            ! do i = 1, z.n+1
            !     tau0.d(i) = sigma_N2 * (8.31 * 100 * Tn.d(i))/(28E-3 * 9.8) * nN2.d(i) + &
            !                 sigma_O2 * (8.31 * 100 * Tn.d(i))/(32E-3 * 9.8) * nO2.d(i) + &
            !                 sigma_O  * (8.31 * 100 * Tn.d(i))/(16E-3 * 9.8) * nO.d(i)
            ! end do

    call p.init(z.n+1)
    call k.init(z.n+1)
    do i = 1, z.n+1
            p.d(i) = ( 4E-7 * nO.d(i) ) * pk_switch
            k.d(i) = ( 1.2E-12 * nN2.d(i) + 2.1E-11 * nO2.d(i) ) * 1E-6 * pk_switch
    end do

    !Diffusion coefficients vector. D.d(i) = D_{i+1/2}
    call D.init(z.n)
    do i = 1, z.n
        D.d(i) = D_analytical(z.d(i)+h0/2)
    end do

    call D_node.init(z.n+1)
    do i = 1, z.n
        D_node.d(i) = D_analytical(z.d(i))
    end do
        D_node.d(z.n+1) = D_analytical(z.d(z.n) + h0)

    call gradD.init(z.n+1)
    do i = 1, z.n+1
        gradD.d(i) = D_node.d(i) * ( gradTp.d(i)/Tp.d(i) - (gradnO.d(i))/nO.d(i) - gradTr.d(i)/(2*Tr.d(i)) )
    end do

    !Effective velocity vector
    !u.d(i) = u_i = D_i * (dTp/dz + mg/2k)/Tp, mg/2k ~ 5.6*10^{-3} [K/m]
    call u.init(z.n+1)
    do i = 1, z.n+1    
        u.d(i) = ( 3E+19 / (nO.d(i) * sqrt(Tr.d(i))) * (56E-4 + gradTp.d(i)) ) * transf_yz_switch
    end do

    call gradu.init(z.n+1)
    do i = 1, z.n+1
        gradu.d(i) = 3E+19 * ( gradgradTp.d(i)/(nO.d(i) * sqrt(Tr.d(i))) - &
                               (56E-4 + gradTp.d(i)) * gradnO.d(i)/ (nO.d(i)**2 * sqrt(Tr.d(i))) - &
                               (56E-4 + gradTp.d(i)) * gradTr.d(i)/2/(nO.d(i) * sqrt(Tr.d(i))**3) )
    end do

        ! write(20, *) 'Ti'
        ! call Ti.print_column(20)
        ! write(20, *) 

        ! write(20, *) 'Tn'
        ! call Tn.print_column(20)
        ! write(20, *) 

        ! write(20, *) 'Te'
        ! call Tr.print_column(20)
        ! write(20, *) 

        ! write(20, *) 'Tr'
        ! call Tr.print_column(20)
        ! write(20, *) 

        ! write(20, *) 'Tp'
        ! call Tp.print_column(20)
        ! write(20, *) 

        ! write(20, *) 'gradTp'
        ! call gradTp.print_column(20)
        ! write(20, *) 

        ! write(20, *) 'gradgradTp'
        ! call gradgradTp.print_column(20)
        ! write(20, *) 

        ! write(20, *) 'gradTr'
        ! call gradTr.print_column(20)
        ! write(20, *) 

        ! write(20, *) 'gradTn'
        ! call gradTn.print_column(20)
        ! write(20, *) 

        ! write(20, *) 'nO'
        ! call nO.print_column(20)
        ! write(20, *) 

        ! write(20, *) 'nO2'
        ! call nO2.print_column(20)
        ! write(20, *) 

        ! write(20, *) 'nN2'
        ! call nN2.print_column(20)
        ! write(20, *) 

        ! write(20, *) 'gradnO'
        ! call gradnO.print_column(20)
        ! write(20, *) 

        ! write(20, *) 'P'
        ! call P.print_column(20)
        ! write(20, *) 

        ! write(20, *) 'k'
        ! call k.print_column(20)
        ! write(20, *) 

        ! write(20, *) 'D'
        ! call D.print_column(20)
        ! write(20, *) 

        ! write(20, *) 'D_node'
        ! call D_node.print_column(20)
        ! write(20, *) 

        ! write(20, *) 'u'
        ! call u.print_column(20)
        ! write(20, *) 

        ! write(20, *) 'gradu'
        ! call gradu.print_column(20)
        ! write(20, *) 



    !if(0) then !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

allocate(z_arr(1:z.n))
do i = 1, z.n
    z_arr(z.n-i+1) = z.d(i) !0.5*(z.d(i)+z.d(i+1))
end do

!call init_transfer(z_arr)

call F_ub.init(Nphi)
do j = 1, Nphi
    call p_mod(j).init(z.n)
    call model_sol(j).init(z.n)
end do


do j = 1, Nphi
    phi = (j-0.5)*dphi-pi/2
    do i = 1, z.n
        model_sol(j).d(i) = f_m(z.d(i), phi)
    end do
end do   


!initialization of n before iterations
do i = 1, Nz
    do j = 1, Nphi
        do lambda = 1, Nlambda
            n((j - 1) * Nz + i, lambda) = model_sol(j).d(i) 
            ans(i, j, lambda) = model_sol(j).d(i)
        enddo
    enddo
enddo

call init_transfer(z_arr,ans)

!if(0 .eq. 1) then

do t = 0, Ndays*3600*24/tau
    if(mod(t, 100) .eq. 0) then
        print *, t
        !write (*, *) "Files written."
    end if



    do j = 1, Nphi
        phi = (j-0.5)*dphi-pi/2
        sI   = sin(atan(2*tan(-pi/2+(j-0.5)*dphi)))
        cI   = cos(atan(2*tan(-pi/2+(j-0.5)*dphi)))

        
        if(model_solution .eq. 1) then
            F_ub.d(j) = D.d(z.n-1)*(sI**2)*df_dz(z.d(z.n)-h0/2, phi) + u_analytical(z.d(z.n)-h0/2)*(sI**2)*f_m(z.d(z.n)-h0/2, phi) - 1/R*sI*cI*D.d(z.n-1)*df_dphi(z.d(z.n)-h0/2, phi)*mixed_z_switch
            F_ub.d(j) = F_ub.d(j)! * (2+cos(2*pi*t*tau/8460))
            do i = 1, z.n
                n_mod = f_m(z.d(i), phi)! * (2+cos(2*pi*t*tau/8460))
                dn_dz = df_dz(z.d(i), phi)! * (2+cos(2*pi*t*tau/8460))
                dn_dphi = df_dphi(z.d(i), phi)! * (2+cos(2*pi*t*tau/8460))
                ddn_dzdz = ddf_dzdz(z.d(i), phi)! * (2+cos(2*pi*t*tau/8460))
                ddn_dzdphi = ddf_dzdphi(z.d(i), phi)! * (2+cos(2*pi*t*tau/8460))
                ddn_dphidphi = ddf_dphidphi(z.d(i), phi)! * (2+cos(2*pi*t*tau/8460))

                ! multiplyer = 1

                !  vphi = vel_phi(z.d(i), phi)*multiplyer
                !  dvphi_dphi = dvelphi_dphi(z.d(i), phi)*multiplyer
                !  vz = vel_z(z.d(i), phi)*multiplyer
                !  dvz_dz = dvelz_dz(z.d(i), phi)*multiplyer

                p_mod(j).d(i) = k.d(i)*n_mod - ( &
                                sI**2 * (gradD.d(i)*dn_dz + ddn_dzdz*D_node.d(i)) + &
                                sI**2 * (gradu.d(i)*n_mod + u.d(i)*dn_dz) - &
                                (gradD.d(i)*dn_dphi + D_node.d(i)*ddn_dzdphi)/R*sI*cI + &
                                D_node.d(i)/(R**2)/cos(phi) * (dA_dphi(phi)*dn_dphi + A(phi)*ddn_dphidphi) - &
                                D_node.d(i)/R/cos(phi)/2    * (dB_dphi(phi)*dn_dz   + B(phi)*ddn_dzdphi) - &
                                u.d(i)/R/cos(phi)/2 * (dB_dphi(phi)*n_mod + B(phi)*dn_dphi)) !&
                                ! +f_m(z.d(i), phi) * (-sin(2*pi*t*tau/8460)*2*pi/8460) &

                !p_mod(j).d(i) = 0

                model_sol(j).d(i) = f_m(z.d(i), phi)
            end do
            p_mod(j).d(1) = 0
        else
            F_ub.d(j) = 0
            do i = 1, z.n     
                p_mod(j).d(i) = p.d(i)
            end do
        end if        
    end do   




    !if(0 .eq. 1) then
    !matrix initialization
    stencil = 0
    do i = 1, Nz
        do j = 1, Nphi

            phi = (j-0.5)*dphi-pi/2
            sI  =  sin(atan(2*tan(    phi   )))
            cI  =  cos(atan(2*tan(    phi   )))
            sIph = sin(atan(2*tan(phi+dphi/2)))
            cIph = cos(atan(2*tan(phi+dphi/2)))
            sImh = sin(atan(2*tan(phi-dphi/2)))
            cImh = cos(atan(2*tan(phi-dphi/2)))

            if (i == 1) then

                !lower boundary: n = P/k at z = 100 km
                operator(+1, :) = [0d0,  0d0  , 0d0]
                operator( 0, :) = [0d0, k.d(i)*tau, 0d0]
                operator(-1, :) = [0d0,  0d0  , 0d0]
                   
            else if (j == 1 .and. i /= Nz) then

                !south pole
                diffusion_transfer_z(+1, :) = [0d0, (-D.d( i )/(h0**2) - 1*u.d(i+1)/(2*h0))*tau*sI**2       , 0d0]
                diffusion_transfer_z( 0, :) = [0d0, 1 + k.d(i)*tau + (D.d(i-1) + D.d(i))*tau/(h0**2)*sI**2, 0d0]
                diffusion_transfer_z(-1, :) = [0d0, (-D.d(i-1)/(h0**2) + 1*u.d(i-1)/(2*h0))*tau*sI**2       , 0d0]


                diffusion_transfer_y(+1, :) = [0d0, 0d0, 0d0]
                diffusion_transfer_y( 0, :) = [0d0, & 
                                                 (D_node.d(i)/R*A(phi+dphi/2)/(dphi**2) - 1*u.d(i)/2*B(phi-dphi)/(2*dphi))*tau/R/cos(phi), &
                                                (-D_node.d(i)/R*A(phi+dphi/2)/(dphi**2) + 1*u.d(i)/2*B(phi+dphi)/(2*dphi))*tau/R/cos(phi)  ]
                diffusion_transfer_y(-1, :) = [0d0, 0d0, 0d0]

                mixed_z(+1, :) = [0d0,  0d0                                        ,  0d0                                        ]
                mixed_z( 0, :) = [0d0, -0.5*D_node.d( i )*tau*sIph*cIph/(R*h0*dphi),  0.5*D_node.d( i )*tau*sIph*cIph/(R*h0*dphi)]
                mixed_z(-1, :) = [0d0,  0.5*D_node.d(i-1)*tau*sIph*cIph/(R*h0*dphi), -0.5*D_node.d(i-1)*tau*sIph*cIph/(R*h0*dphi)]

                mixed_y(+1, :) = [0d0,  0d0                                           ,  0d0                                                ]
                mixed_y( 0, :) = [0d0, -0.5*D.d(i-1)/2*B(phi)*tau/(h0*dphi*R*cos(phi)),  0.5*D.d(i-1)/2*B(phi+dphi)*tau/(h0*dphi*R*cos(phi))]
                mixed_y(-1, :) = [0d0,  0.5*D.d(i-1)/2*B(phi)*tau/(h0*dphi*R*cos(phi)), -0.5*D.d(i-1)/2*B(phi+dphi)*tau/(h0*dphi*R*cos(phi))]

                do s_i = -1, 1
                    do s_j = -1, 1
                        operator(s_i, s_j) = diffusion_transfer_z(s_i, s_j) + mixed_z_switch*mixed_z(s_i, s_j) + diffusion_transfer_y(s_i, s_j) + mixed_y_switch*mixed_y(s_i, s_j)
                    enddo
                enddo

            else if (j == Nphi .and. i /= Nz) then

                !north pole
                diffusion_transfer_z(+1, :) = [0d0, (-D.d( i )/(h0**2) - 1*u.d(i+1)/(2*h0))*tau*sI**2       , 0d0]
                diffusion_transfer_z( 0, :) = [0d0, 1 + k.d(i)*tau + (D.d(i-1) + D.d(i))*tau/(h0**2)*sI**2, 0d0]
                diffusion_transfer_z(-1, :) = [0d0, (-D.d(i-1)/(h0**2) + 1*u.d(i-1)/(2*h0))*tau*sI**2       , 0d0]

                diffusion_transfer_y(+1, :) = [0d0, 0d0, 0d0]
                diffusion_transfer_y( 0, :) = [ (-D_node.d(i)/R*A(phi-dphi/2)/(dphi**2) - 1*u.d(i)/2*B(phi-dphi)/(2*dphi))*tau/R/cos(phi), &
                                                 (D_node.d(i)/R*A(phi-dphi/2)/(dphi**2) + 1*u.d(i)/2*B(phi+dphi)/(2*dphi))*tau/R/cos(phi), &
                                               0d0]
                diffusion_transfer_y(-1, :) = [0d0, 0d0, 0d0]

                mixed_z(+1, :) = [ 0d0                                        ,  0d0                                        , 0d0]
                mixed_z( 0, :) = [-0.5*D_node.d( i )*tau*sImh*cImh/(R*h0*dphi),  0.5*D_node.d( i )*tau*sImh*cImh/(R*h0*dphi), 0d0]
                mixed_z(-1, :) = [ 0.5*D_node.d(i-1)*tau*sImh*cImh/(R*h0*dphi), -0.5*D_node.d(i-1)*tau*sImh*cImh/(R*h0*dphi), 0d0]

                mixed_y(+1, :) = [ 0d0                                                ,  0d0                                           ,  0d0]
                mixed_y( 0, :) = [-0.5*D.d(i-1)/2*B(phi-dphi)*tau/(h0*dphi*R*cos(phi)),  0.5*D.d(i-1)/2*B(phi)*tau/(h0*dphi*R*cos(phi)),  0d0]
                mixed_y(-1, :) = [ 0.5*D.d(i-1)/2*B(phi-dphi)*tau/(h0*dphi*R*cos(phi)), -0.5*D.d(i-1)/2*B(phi)*tau/(h0*dphi*R*cos(phi)),  0d0]

                do s_i = -1, 1
                    do s_j = -1, 1
                        operator(s_i, s_j) = diffusion_transfer_z(s_i, s_j) + mixed_z_switch*mixed_z(s_i, s_j) + diffusion_transfer_y(s_i, s_j) + mixed_y_switch*mixed_y(s_i, s_j)
                    enddo
                enddo

            else if (i == Nz .and. upper_bound_type == 1) then

                diffusion_transfer_z(+1, :) = [0d0, 0d0, 0d0]
                diffusion_transfer_z( 0, :) = [0d0, ( D.d(i-1)/(h0**2) + 1*u.d( i )/(2*h0))*tau*sI**2, 0d0] 
                diffusion_transfer_z(-1, :) = [0d0, (-D.d(i-1)/(h0**2) + 1*u.d(i-1)/(2*h0))*tau*sI**2, 0d0] 

                if(j == 1 .or. j == Nphi) then
                    if(sI*cI .ge. 0) then
                        mixed_z(+1, :) = [ 0d0                                        ,  0d0                                        , 0d0]
                        mixed_z( 0, :) = [ 0d0                                        ,  0d0                                        , 0d0]
                        mixed_z(-1, :) = [ 0.5*D.d(i-1)*tau*sImh*cImh/(R*h0*dphi), -0.5*D.d(i-1)*tau*sImh*cImh/(R*h0*dphi), 0d0]
                    else
                        mixed_z(+1, :) = [ 0d0,  0d0                                        ,  0d0                                        ]
                        mixed_z( 0, :) = [ 0d0,  0d0                                        ,  0d0                                        ]
                        mixed_z(-1, :) = [ 0d0,  0.5*D.d(i-1)*tau*sIph*cIph/(R*h0*dphi), -0.5*D.d(i-1)*tau*sIph*cIph/(R*h0*dphi)]
                    end if
                else
                    if(sI*cI .ge. 0) then
                        mixed_z(+1, :) = [ 0d0                                   ,  0d0                                   ,  0d0                                   ]
                        mixed_z( 0, :) = [ 0d0                                   ,  0.5*D.d(i-1)*tau*sIph*cIph/(R*h0*dphi), -0.5*D.d(i-1)*tau*sIph*cIph/(R*h0*dphi)]
                        mixed_z(-1, :) = [ 0.5*D.d(i-1)*tau*sImh*cImh/(R*h0*dphi), -0.5*D.d(i-1)*tau*sImh*cImh/(R*h0*dphi),  0d0                                   ]
                    else
                        mixed_z(+1, :) = [ 0d0                                   ,  0d0                                   ,  0d0                                   ]
                        mixed_z( 0, :) = [ 0.5*D.d(i-1)*tau*sImh*cImh/(R*h0*dphi), -0.5*D.d(i-1)*tau*sImh*cImh/(R*h0*dphi),  0d0                                   ]
                        mixed_z(-1, :) = [ 0d0                                   ,  0.5*D.d(i-1)*tau*sIph*cIph/(R*h0*dphi), -0.5*D.d(i-1)*tau*sIph*cIph/(R*h0*dphi)]
                    end if
                end if

                do s_i = -1, 1
                    do s_j = -1, 1
                        operator(s_i, s_j) = diffusion_transfer_z(s_i, s_j) + diffusion_transfer_y(s_i, s_j) + mixed_z_switch*mixed_z(s_i, s_j) + mixed_y_switch*mixed_y(s_i, s_j)
                    enddo
                enddo

            else

                !the main part of the operator

                diffusion_transfer_z(+1, :) = [0d0, (-D.d( i )/(h0**2) - 1*u.d(i+1)/(2*h0))*tau*sI**2        , 0d0]
                diffusion_transfer_z( 0, :) = [0d0, 1 + k.d(i)*tau + (D.d(i-1) + D.d(i))*tau/(h0**2)*sI**2 , 0d0]
                diffusion_transfer_z(-1, :) = [0d0, (-D.d(i-1)/(h0**2) + 1*u.d(i-1)/(2*h0))*tau*sI**2        , 0d0]

                diffusion_transfer_y(+1, :) = [0d0, 0d0, 0d0]
                diffusion_transfer_y( 0, :) = [(-D_node.d(i)/R*A(phi-dphi/2)/(dphi**2) + (-1*u.d(i)/2)*B(phi-dphi)/(2*dphi))*tau/R/cos(phi), &
                                                (D_node.d(i)/R *(A(phi-dphi/2) + A(phi+dphi/2))/(dphi**2))*tau/(R)/cos(phi), &
                                               (-D_node.d(i)/R*A(phi+dphi/2)/(dphi**2) - (-1*u.d(i)/2)*B(phi+dphi)/(2*dphi))*tau/R/cos(phi)]
                diffusion_transfer_y(-1, :) = [0d0, 0d0, 0d0]

                if(sI .ge. 0) then
                    mixed_z(+1, :) = [ 0d0,                                         -0.5*D_node.d(i+1)*tau*sIph*cIph/(R*h0*dphi)          ,  0.5*D_node.d(i+1)*tau*sIph*cIph/(R*h0*dphi)]
                    mixed_z( 0, :) = [-0.5*D_node.d( i )*tau*sImh*cImh/(R*h0*dphi),  0.5*D_node.d(i)*tau*(sIph*cIph+sImh*cImh)/(R*h0*dphi), -0.5*D_node.d( i )*tau*sIph*cIph/(R*h0*dphi)]
                    mixed_z(-1, :) = [ 0.5*D_node.d(i-1)*tau*sImh*cImh/(R*h0*dphi), -0.5*D_node.d(i-1)*tau*sImh*cImh/(R*h0*dphi)          ,  0d0                                        ]

                    mixed_y(+1, :) = [ 0d0                                                , -0.5*D.d( i )/2*B(phi)*tau/(h0*dphi*R*cos(phi))         ,  0.5*D.d( i )/2*B(phi+dphi)*tau/(h0*dphi*R*cos(phi))]
                    mixed_y( 0, :) = [-0.5*D.d(i-1)/2*B(phi-dphi)*tau/(h0*dphi*R*cos(phi)),  0.5*(D.d(i)+D.d(i-1))/2*B(phi)*tau/(h0*dphi*R*cos(phi)), -0.5*D.d( i )/2*B(phi+dphi)*tau/(h0*dphi*R*cos(phi))]
                    mixed_y(-1, :) = [ 0.5*D.d(i-1)/2*B(phi-dphi)*tau/(h0*dphi*R*cos(phi)), -0.5*D.d(i-1)/2*B(phi)*tau/(h0*dphi*R*cos(phi))         ,  0d0                                                ]
                else
                    mixed_z(+1, :) = [-0.5*D_node.d(i+1)*tau*sImh*cImh/(R*h0*dphi),  0.5*D_node.d(i+1)*tau*sImh*cImh/(R*h0*dphi)          ,  0d0                                          ]
                    mixed_z( 0, :) = [ 0.5*D_node.d( i )*tau*sImh*cImh/(R*h0*dphi), -0.5*D_node.d(i)*tau*(sImh*cImh+sIph*cIph)/(R*h0*dphi),  0.5*D_node.d( i )*tau*sIph*cIph/(R*h0*dphi)]
                    mixed_z(-1, :) = [ 0d0,                                          0.5*D_node.d(i-1)*tau*sIph*cIph/(R*h0*dphi)          , -0.5*D_node.d(i-1)*tau*sIph*cIph/(R*h0*dphi)]

                    mixed_y(+1, :) = [-0.5*D.d( i )/2*B(phi-dphi)*tau/(h0*dphi*R*cos(phi)),  0.5*D.d( i )/2*B(phi)*tau/(h0*dphi*R*cos(phi))         ,  0d0                                                ]
                    mixed_y( 0, :) = [ 0.5*D.d( i )/2*B(phi-dphi)*tau/(h0*dphi*R*cos(phi)), -0.5*(D.d(i)+D.d(i-1))/2*B(phi)*tau/(h0*dphi*R*cos(phi)),  0.5*D.d(i-1)/2*B(phi+dphi)*tau/(h0*dphi*R*cos(phi))]
                    mixed_y(-1, :) = [ 0d0,                                                  0.5*D.d(i-1)/2*B(phi)*tau/(h0*dphi*R*cos(phi)),          -0.5*D.d(i-1)/2*B(phi+dphi)*tau/(h0*dphi*R*cos(phi))]
                end if

                do s_i = -1, 1
                    do s_j = -1, 1
                        operator(s_i, s_j) = diffusion_transfer_z(s_i, s_j) + diffusion_transfer_y(s_i, s_j) + mixed_z_switch*mixed_z(s_i, s_j) + mixed_y_switch*mixed_y(s_i, s_j)
                    enddo
                enddo

            end if

            do s_i = -1, 1
                do s_j = -1, 1
                    if (i + s_i >= 1 .and. i + s_i <= Nz .and. j + s_j >= 1 .and. j + s_j <= Nphi) then
                        stencil(i, j, s_i, s_j) = operator(s_i, s_j)
                    endif
                enddo
            enddo

        end do 
    end do

    !forming the RHS
    if(t*tau < diurnal_start_time .or. diurnal_on .eq. 0) then
        do lambda = 1, Nlambda-1
            do i = 1, Nz
                do j = 1, Nphi
                    if (i == 1) then
                        rhs(i, j, lambda) = p_mod(j).d(i)*tau
                    else if (i == Nz) then
                        rhs(i, j, lambda) = F_ub.d(j) * tau/h0
                    else
                        rhs(i, j, lambda) = tau * p_mod(j).d(i) + n((j - 1) * Nz + i, lambda)
                    end if
                enddo
            enddo
        enddo
    else
        ! do lambda = 1, Nlambda-1
        !     do i = 1, Nz
        !         do j = 1, Nphi
        !             day = (t*tau)/86400 + 1/2 + (lambda-1)/(1.*(Nlambda-1)) !starting from the middle of the 1-st day
        !             tgdelta = tan(pi/180*23.5) * sin(2*pi/365 * (day - 80))
        !             sindelta = tgdelta/sqrt(1+tgdelta**2)
        !             cosdelta = sqrt(1-sindelta**2) !cos of the zenith angle is > 0
        !             coschi = sin(-pi/2+(j-0.5)*dphi)*sindelta - &
        !                      cos(-pi/2+(j-0.5)*dphi)*cosdelta*cos(omega*( tau*t + 86400/2 + 86400*(lambda-1)/(1.*(Nlambda-1)) )) !start: middle of the 1 da
        !             if(coschi > 1E-6) then
        !                 if (i == 1) then
        !                     rhs(i, j, lambda) = p_mod(j).d(i) * exp(tau0.d(i) * (1-1/coschi)) 
        !                 else if (i == Nz) then
        !                     rhs(i, j, lambda) = F_ub.d(j) * tau/h0
        !                 else
        !                     rhs(i, j, lambda) = tau * p_mod(j).d(i) * exp(tau0.d(i) * (1-1/coschi)) + ans(i, j, lambda)
        !                 end if
        !             else
        !                 if (i == 1) then
        !                     rhs(i, j, lambda) = 1 
        !                 else if (i == Nz) then
        !                     rhs(i, j, lambda) = F_ub.d(j) * tau/h0
        !                 else
        !                     rhs(i, j, lambda) = ans(i, j, lambda)
        !                 end if
        !             end if
        !         enddo
        !     enddo
        ! enddo
    end if

    !CSR matrix forming + RHS for the 2d system
    do lambda = 1, Nlambda-1    
        counter_a = 1
        counter_ia = 1
        counter_rhs = 1
        ia(1) = 1
        do j = 1, Nphi
            do i = 1, Nz

                do s_j = -1, 1
                    do s_i = -1, 1
                        if (stencil(i, j, s_i, s_j) /= 0d0) then
                            arr(counter_a) = stencil(i, j, s_i, s_j)
                            ja(counter_a) = (j + s_j - 1) * Nz + (i + s_i - 1) + 1
                            counter_a = counter_a + 1
                        endif
                    end do 
                end do

                ia(counter_rhs + 1) = counter_a

                f(counter_rhs, lambda) = rhs(i, j, lambda)
                counter_rhs = counter_rhs + 1

            end do
        end do
        matrix_size = Nz * Nphi
        nonzero = counter_a - 1
    enddo

    !Preconditioning
    if (t == 0) then
        ipalu = 1
        ipbcg = ipalu + nonzero
        ipju = 1
        ipjlu = ipju + matrix_size + 1
        ipiw = ipjlu + nonzero

        call ilu0(matrix_size, arr, ja, ia, rw(ipalu), iw(ipjlu), iw(ipju), iw(ipiw), ierr)
        !write(*, *) ierr
    end if

    !Iterative solving (bi-conjugate gradients)
    do lambda = 1, Nlambda-1
        resinit = 0
        !n = 0d0
        do i = 1, Nz
            do j = 1, Nphi
                eps = rhs(i, j, lambda)
                do s_i = -1, 1
                    do s_j = -1, 1
                        if (stencil(i, j, s_i, s_j) /= 0) then
                            eps = eps - stencil(i, j, s_i, s_j) * n((j+s_j - 1) * Nz + i+s_i, lambda)
                        endif
                    enddo
                enddo
                if (abs(eps) > resinit) resinit = abs(eps);
            enddo
        enddo
        if(printing_mode .eq. 1) then
            write(*, *) "resinit", resinit
        end if

        ITER = 1000
        RESID = 1d-10 * resinit
        INFO = 0
        if(printing_mode .eq. 1) then
            NUNIT = 6 ! 6 to output
        else
            NUNIT = 0
        end if

        iprevec(1) = matrix_size
        imatvec(1) = matrix_size

        call slpbcgs(prevec0, iprevec, iw, rw,   &
                     matvec, imatvec, ia, ja, arr, &
                     rw(ipbcg), matrix_size, 8, matrix_size, f(:, lambda), n(:, lambda),   &
                     ITER, RESID, INFO, NUNIT)
        if(printing_mode .eq. 1) then
            print *, 'INFO =', INFO
        end if



        do i = 1, Nz
            do j = 1, Nphi
                ans_prev(i, j, lambda) = ans(i, j, lambda) 
                ans(i, j, lambda)      = n((j - 1) * Nz + i, lambda)
            enddo
        enddo
        if(model_solution .eq. 1) then
            do j = 1, Nphi
                ans(1, j, lambda) = 0
                n((j - 1) * Nz + 1, lambda) = 0
                ans_prev(1, j, lambda) = 0
            end do
        end if


        resinit = 0

        do i = 1, Nz
            do j = 1, Nphi
                eps = rhs(i, j, lambda)
                do s_i = -1, 1
                    do s_j = -1, 1
                        if (stencil(i, j, s_i, s_j) /= 0) then
                            eps = eps - stencil(i, j, s_i, s_j) * ans(i + s_i, j + s_j, lambda)
                        endif
                    enddo
                enddo
                if (abs(eps) > resinit) resinit = abs(eps);
            enddo
        enddo


        if(printing_mode .eq. 1) then
            print *, 'Residual = ', resinit
        end if
    end do

    do i = 1, Nz
        do j = 1, Nphi
            ans_prev(i, j, Nlambda) = ans(i, j, Nlambda)
            ans(i, j, Nlambda)      = ans(i, j, 1)
        enddo
    enddo

    if(model_solution .eq. 1) then
        do j = 1, Nphi
            ans(1, j, Nlambda) = 0
            n((j - 1) * Nz + 1, Nlambda) = 0
            ans_prev(1, j, Nlambda) = 0
        end do
    end if

    !print*,'in diffusion: Nz=',Nz,'Nphi=',Nphi,'Nlambda=',Nlambda
    ans_inc=ans-ans_prev
    call step_of_transfer(ans_inc(1:Nz, 1:Nphi, 1:Nlambda-1), tau) !ans at z = 500 is not sent to the transfer code
    do i = 1, Nz
        do j = 1, Nphi
            do lambda = 1, Nlambda-1
                ans(i, j, lambda) = uout(lambda, j, Nz-i+1)
            enddo
        enddo
    enddo

    do j = 1, Nphi
        do lambda = 1, Nlambda
            ans(1, j, lambda) = 0
        end do
    end do

    do i = 1, Nz
        do j = 1, Nphi
            ans(i, j, Nlambda) = ans(i, j, 1)
        enddo
    enddo

    !endif

    !Printing output
    if(mod(t, 2) == 0 .and. printing_mode == 1) then

        open(unit=11, name='res_implicit.txt')
        open(unit=12, name='res_implicit_gnp.txt')
        open(unit=30, name='diff.txt')
        open(unit=35, name='vel_z.txt')

        !open(unit=13, name='res_diff_gnp.txt')
        open(unit=36, name='p_mod.txt')

        do j = 1, Nphi
            do i = 1, z.n
                write(12,*) (j-5E-1)*180/(Nphi)-90, 100+(Hmax - 100)/(z.n-1)*(i-1), ans(i, j, 1)
            end do
            write (12, *)
        end do

        do j = 1, Nphi
            do i = 1, z.n
                write(30,*) (j-5E-1)*180/(Nphi)-90, 100+(Hmax - 100)/(z.n-1)*(i-1), model_sol(j).d(i) - ans(i, j, 1)
            end do
            write (30, *)
        end do

        do j = 1, Nphi
            do i = 1, z.n-1
                write(35,*) (j-5E-1)*180/(Nphi)-90, 100+(Hmax - 100)/(z.n-1)*(i-1), vel_z(z.d(i), (j-0.5)*dphi-pi/2)
            end do
            write (35, *)
        end do

        do j = 1, Nphi
            phi = (j-0.5)*dphi-pi/2
            do i = 1, z.n-1
                write(36,*) (j-5E-1)*180/(Nphi)-90, 100+(Hmax - 100)/(z.n-1)*(i-1), p_mod(j).d(i)
            end do
            write (36, *)
        end do

        do j = 1, Nphi
            write (11, '(1000(e10.3))') (ans(i, j, 1), i = 1, z.n)
        end do

        close(11)
        close(12)
        close(30)
        close(35)
        close(36)
        close(37)
    end if

    if(convergence_test .eq. 1 .and. mod(t, 10) == 0) then
        do i = 1, z.n-1
            write(50,*) t, 100+(Hmax - 100)/(z.n-1)*(i-1), ans(i, 65, 1)
        enddo
        write (50, *)
    end if

    if(convergence_test .eq. 1 .and. mod(t, 1) == 0) then
        do j = 1, Nphi
            write(60,*) t, j, ans(1, j, 1)
        enddo
        write (60, *)
    end if

end do

!end if!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



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
    nO = 2.8E+10 * exp(-9.8 * 16E-3 / (8.31 * Tn) * (z - 140000)) * 1E+6
    gradTp = 0.5 * ((Te0 - Te) * 9.8 / (287 * Te0) + (Ti0 - Ti) * 9.8 / (287 * Ti0) )
    u_analytical =  3E+17*1E+2 / (nO * sqrt(Tr)) * (56E-4 + gradTp) 
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
    nO = 2.8E+10 * exp(-9.8 * 16E-3 / (8.31 * Tn) * (z - 140000)) * 1E+6
    D_analytical =  3E+17 * Tp / (nO * sqrt(Tr)) * 1E+2
    return
end




