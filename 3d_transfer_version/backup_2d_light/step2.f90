program mat

use vector

implicit none

type (vect) F_ub, z, h, hmid, nO, nO2, nN2, k, p, pcurr, m, D, D_node, u, Tn, Ti, Te, Tr, Tp, tau0, model_sol(181), p_mod(181)
integer convergence_profile_output, latitude_node, continuous_output, model_solution, printing_mode, s_i, s_j, i, j, t, Te0, Tn0, Ti0, day, profile_output, diurnal_on, convergence_test, pk_switch, mixed_z_switch, mixed_y_switch, transf_yz_switch, transf_y_switch, upper_bound_type
real (8) coef, D_analytical, u_analytical, df_dz, df_dphi, ddf_dzdz, ddf_dphidphi, ddf_dzdphi, f_m, dA_dphi, dB_dphi
real (8) n_max, h_max, tau, Hmax, h0, F_z, eps, tgdelta, sindelta, cosdelta, dphi, phi, coschi, pi, omega, sigma_O2, sigma_N2, sigma_O, sI, cI, R, A, B, Ndays, Niter, sIph, cIph, sImh, cImh
type (vect) gradu, gradD, gradTp, gradTr, gradTi, gradTn, gradTe, gradgradTp, gradnO

integer, parameter ::   max_size  =  10000000
integer, parameter :: max_nonzero = 100000000
integer, parameter :: maxwr = max_nonzero + 8 * max_size
integer, parameter :: maxwi = max_nonzero + 2 * max_size + 1
integer, parameter ::  Nz  =  81
integer, parameter :: Nphi = 180


integer, allocatable:: ia(:), ja(:), iw(:)
double precision, allocatable:: arr(:), f(:), n(:), rw(:)


external matvec, prevec0
integer ITER, INFO, NUNIT, ierr, ipalu, ipjlu, ipju, ipiw
real*8 RESID, ddot

external ddot, dcopy

integer imatvec(1), iprevec(1), ipbcg, matrix_size, nonzero
real*8 resinit


integer counter_a, counter_ia, counter_rhs
real*8 C_norm, L2_norm


real*8 stencil(Nz, Nphi, -1:1, -1:1), operator(-1:1, -1:1), diffusion_transfer_z(-1:1, -1:1), diffusion_transfer_y(-1:1, -1:1), mixed_z(-1:1, -1:1), mixed_y(-1:1, -1:1)
real*8 rhs_z_phi(Nz, Nphi), ans(Nz, Nphi), n_mean(Nz, Nphi), rhs_z_phi_analytical(Nz, Nphi)

allocate(ia(max_size + 1))
allocate(ja(max_nonzero))
allocate(iw(maxwi))
allocate(arr(max_nonzero))
allocate(f(max_size))
allocate(n(max_size))
allocate(rw(maxwr))


pi = 3.141592653589793238462643

!nonlinear_scheme_type variable switches the u_phi-approximation. 
!nonlinear_scheme_type = 8
 
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
Ndays = 500
Niter = 800
!upper boundary electron flow
F_z = 0

!Time step (in seconds) 5 min
tau =  60

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
upper_bound_type = 1
diurnal_on = 1
convergence_test = 0
profile_output = 0
printing_mode = 1
model_solution = 0
continuous_output = 0
convergence_profile_output = 0

    if(continuous_output .eq. 1) then
        !open(unit=988, name='evolution_matrix_88.txt')
        open(unit=930, name='evolution_matrix_30.txt')
        open(unit=960, name='evolution_matrix_60.txt')
        open(unit=900, name='evolution_matrix_0.txt')
        latitude_node = 2 !this is -60
    end if

    if(convergence_profile_output .eq. 1) then
        open(unit=1030, name='convergence_C_norm_gnp.txt')
        open(unit=1040, name='convergence_C_norm_60_gnp.txt')
    end if




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

    !P = P_0 exp(tau_0(z)*(1-sec chi))
    !tau_0 = sum_{N2, O2, O} sigma_i R_0*T_n/(M_i g) n_i(z)
    !sigma is in cm^2
    sigma_O2 = 2E-17 
    sigma_N2 = 15E-18
    sigma_O  = 1E-17
    call tau0.init(z.n+1)
    do i = 1, z.n+1
        tau0.d(i) = sigma_N2 * (8.31 * 100 * Tn.d(i))/(28E-3 * 9.8) * nN2.d(i) + &
                    sigma_O2 * (8.31 * 100 * Tn.d(i))/(32E-3 * 9.8) * nO2.d(i) + &
                    sigma_O  * (8.31 * 100 * Tn.d(i))/(16E-3 * 9.8) * nO.d(i)
    end do


    call p.init(z.n+1)
    call k.init(z.n+1)
    do i = 1, z.n+1
            p.d(i) = ( 4E-7 * nO.d(i) ) * pk_switch
            k.d(i) = ( 1.2E-12 * nN2.d(i) + 2.1E-11 * nO2.d(i) ) * pk_switch
    end do


    !Diffusion coefficients vector. D.d(i) = D_{i+1/2}
    call D.init(z.n)
    do i = 1, z.n
        D.d(i) = D_analytical(z.d(i)+h0/200)
    end do

    call D_node.init(z.n+1)
    do i = 1, z.n
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


call F_ub.init(Nphi)
do j = 1, Nphi
    call p_mod(j).init(z.n)
    call model_sol(j).init(z.n)
end do

do j = 1, Nphi
    phi = (j-0.5)*dphi-pi/2
    sI   = sin(atan(2*tan(-pi/2+(j-0.5)*dphi)))
    cI   = cos(atan(2*tan(-pi/2+(j-0.5)*dphi)))

    
    if(model_solution .eq. 1) then
        F_ub.d(j) = D.d(z.n-1)*(sI**2)*df_dz(z.d(z.n)-h0/200, phi) + u_analytical(z.d(z.n)-h0/200)*(sI**2)*f_m(z.d(z.n)-h0/200, phi) - 1/R*sI*cI*D.d(z.n-1)*df_dphi(z.d(z.n)-h0/200, phi)*mixed_z_switch
        F_ub.d(j) = F_ub.d(j)
        do i = 1, z.n     
            p_mod(j).d(i) = k.d(i)*f_m(z.d(i), phi) - ( &
                            sI**2 * (gradD.d(i)*df_dz(z.d(i), phi) + ddf_dzdz(z.d(i), phi)*D_node.d(i)) + &
                            sI**2 * (gradu.d(i)*f_m(z.d(i), phi) + u.d(i)*df_dz(z.d(i), phi)) - &
                            mixed_z_switch * 1/R*sI*cI * (gradD.d(i)*df_dphi(z.d(i), phi) + D_node.d(i)*ddf_dzdphi(z.d(i), phi)) + &
                            D_node.d(i)/(R**2)/cos(phi) * (dA_dphi(phi)*df_dphi(z.d(i), phi) + A(phi)*ddf_dphidphi(z.d(i), phi)) - &
                            mixed_y_switch * D_node.d(i)/R/cos(phi)/2 * (dB_dphi(phi)*df_dz(z.d(i), phi) + B(phi)*ddf_dzdphi(z.d(i), phi)) - &
                            u.d(i)/R/cos(phi)/2 * (dB_dphi(phi)*f_m(z.d(i), phi) + B(phi)*df_dphi(z.d(i), phi)) &
                            +0)
            p_mod(j).d(i) = p_mod(j).d(i)*1
        end do
    else
        F_ub.d(j) = 0
        do i = 1, z.n     
            p_mod(j).d(i) = p.d(i)
        end do
    end if        
end do    




!initialization of n before iterations
if(model_solution .eq. 1) then
    do i = 1, Nz
        do j = 1, Nphi
            n((j - 1) * Nz + i) = 1
            ans(i, j) = 1
        enddo
    enddo
else
    do i = 1, Nz
        do j = 1, Nphi
            n((j - 1) * Nz + i) = 1 
            ans(i, j) = 1
        enddo
    enddo
end if


!    n = 1d0

if(0 .eq. 0) then

do t = 0, Ndays*86400/tau
    if(mod(t, 100) .eq. 0) then
        print *, t
    end if


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
                operator( 0, :) = [0d0, k.d(i), 0d0]
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
    if(t*tau < 86400 .or. diurnal_on .eq. 0) then
        do i = 1, Nz
            do j = 1, Nphi
                if (i == 1) then
                    rhs_z_phi(i, j) = p_mod(j).d(i)
                else if (i == Nz) then
                    rhs_z_phi(i, j) = F_ub.d(j) * tau/h0
                else
                    rhs_z_phi(i, j) = tau * p_mod(j).d(i) + n((j - 1) * Nz + i)
                end if
            enddo
        enddo
    else
        do i = 1, Nz
            do j = 1, Nphi
                day = (t*tau)/86400 + 1/2 !starting from the middle of the 1-st day
                tgdelta = tan(pi/180*23.5) * sin(2*pi/365 * (day - 80))
                sindelta = tgdelta/sqrt(1+tgdelta**2)
                cosdelta = sqrt(1-sindelta**2) !cos of the zenith angle is > 0
                coschi = sin(-pi/2+(j-0.5)*dphi)*sindelta - &
                         cos(-pi/2+(j-0.5)*dphi)*cosdelta*cos(omega*(tau*t+86400/2)) !start: middle of the 1 da
                if(coschi > 1E-6) then
                    if (i == 1) then
                        rhs_z_phi(i, j) = p_mod(j).d(i) * exp(tau0.d(i) * (1-1/coschi)) 
                    else if (i == Nz) then
                        rhs_z_phi(i, j) = F_ub.d(j) * tau/h0
                    else
                        rhs_z_phi(i, j) = tau * p_mod(j).d(i) * exp(tau0.d(i) * (1-1/coschi)) + ans(i, j)
                    end if
                else
                    if (i == 1) then
                        rhs_z_phi(i, j) = 1 
                    else if (i == Nz) then
                        rhs_z_phi(i, j) = F_ub.d(j) * tau/h0
                    else
                        rhs_z_phi(i, j) = ans(i, j)
                    end if
                end if
            enddo
        enddo
    end if

    if(t == 0 .and. model_solution == 1) then 
        do j = 1, Nphi
            do i = 1, Nz
                if (i == 1) then
                    rhs_z_phi_analytical(i, j) = p_mod(j).d(i)
                else if (i == Nz) then
                    rhs_z_phi_analytical(i, j) = F_ub.d(j) * tau/h0
                else
                    rhs_z_phi_analytical(i, j) = tau * p_mod(j).d(i) + f_m(z.d(i), (j-0.5)*dphi-pi/2)
                end if
            enddo
        enddo

        open(unit=4, name='residue.txt')
        do j = 1, Nphi
            do i = 1, Nz
                eps = 0
                do s_i = -1, 1
                    do s_j = -1, 1
                        if (stencil(i, j, s_i, s_j) /= 0) then
                            eps = eps + stencil(i, j, s_i, s_j) * f_m(z.d(i+s_i), (j+s_j-0.5)*dphi-pi/2)
                        endif
                    enddo
                enddo
                write(4,*) (j-5E-1)*180/(Nphi)-90, 100+(Hmax - 100)/(z.n-1)*(i-1), eps - rhs_z_phi_analytical(i, j)
            enddo
            write (4, *)
        enddo
        close(unit=4)
    end if

    !CSR matrix forming + RHS for the 2d system
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

            f(counter_rhs) = rhs_z_phi(i, j)
            counter_rhs = counter_rhs + 1

        end do
    end do
    matrix_size = Nz * Nphi
    nonzero = counter_a - 1


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
        resinit = 0
        do i = 1, Nz
            do j = 1, Nphi
                eps = rhs_z_phi(i, j)
                do s_i = -1, 1
                    do s_j = -1, 1
                        if (stencil(i, j, s_i, s_j) /= 0) then
                            eps = eps - stencil(i, j, s_i, s_j) * n((j+s_j - 1) * Nz + i+s_i)
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
                     rw(ipbcg), matrix_size, 8, matrix_size, f, n,   &
                     ITER, RESID, INFO, NUNIT)
        if(printing_mode .eq. 1) then
            print *, 'INFO =', INFO
        end if

        do i = 1, Nz
            do j = 1, Nphi
                ans(i, j) = n((j - 1) * Nz + i)
            enddo
        enddo

        resinit = 0
        do i = 1, Nz
            do j = 1, Nphi
                eps = rhs_z_phi(i, j)
                do s_i = -1, 1
                    do s_j = -1, 1
                        if (stencil(i, j, s_i, s_j) /= 0) then
                            eps = eps - stencil(i, j, s_i, s_j) * ans(i + s_i, j + s_j)
                        endif
                    enddo
                enddo
                if (abs(eps) > resinit) resinit = abs(eps);
            enddo
        enddo


        if(printing_mode .eq. 1) then
            print *, 'Residual = ', resinit
        end if

    !Printing output
    if(mod(t, 100) == 0 .and. printing_mode == 1) then


        if(model_solution == 1) then
            open(unit=98, name='model_solution_gnp.txt')    
            open(unit=99, name='model_solution_matrix.txt')
            open(unit=198,name='diff_res_sol_gnp.txt')
            do j = 1, Nphi
                do i = 1, z.n
                    write(98,*) (j-5E-1)*180/(Nphi)-90, 100+(Hmax - 100)/(z.n-1)*(i-1), f_m(z.d(i), (j-5E-1)*pi/(Nphi)-pi/2)
                    model_sol(j).d(i) = f_m(z.d(i), (j-5E-1)*pi/(Nphi)-pi/2)
                end do
                write (98, *)
            end do
            do j = 1, Nphi
                call model_sol(j).print(99)
            end do

            do j = 1, Nphi
                do i = 1, z.n
                    write(198,*) (j-5E-1)*180/(Nphi)-90, 100+(Hmax - 100)/(z.n-1)*(i-1), ans(i, j)-model_sol(j).d(i)
                end do
                write (198, *)
            end do
            close(98)
            close(99)
            close(198)


            open(unit=11, name='res_implicit_model.txt')
            open(unit=12, name='res_implicit_model_gnp.txt')

            do j = 1, Nphi
                do i = 1, z.n
                    write(12,*) (j-5E-1)*180/(Nphi)-90, 100+(Hmax - 100)/(z.n-1)*(i-1), ans(i, j)
                end do
                write (12, *)
            end do

            do j = 1, Nphi
                write (11, '(1000(e10.3))') (ans(i, j), i = 1, z.n)
            end do

            close(11)
            close(12)
        else
            open(unit=11, name='res_implicit.txt')
            open(unit=12, name='res_implicit_gnp.txt')

            do j = 1, Nphi
                do i = 1, z.n
                    write(12,*) (j-5E-1)*180/(Nphi)-90, 100+(Hmax - 100)/(z.n-1)*(i-1), ans(i, j)
                end do
                write (12, *)
            end do

            do j = 1, Nphi
                write (11, '(1000(e10.3))') (ans(i, j), i = 1, z.n)
            end do

            close(11)
            close(12)
        end if

        !profiles output block
        if(profile_output .eq. 1) then
            open(unit=88, name='step2_impl_88.txt')
            open(unit=80, name='step2_impl_80.txt')
            open(unit=70, name='step2_impl_70.txt')
            open(unit=60, name='step2_impl_60.txt')
            open(unit=50, name='step2_impl_50.txt')
            open(unit=40, name='step2_impl_40.txt')
            open(unit=30, name='step2_impl_30.txt')
            open(unit=20, name='step2_impl_20.txt')
            open(unit=10, name='step2_impl_10.txt')
            open(unit=5,  name='step2_impl_05.txt')
            open(unit=2,  name='step2_impl_02.txt')
            open(unit=0,  name='step2_impl_00.txt')
            do i = 1, z.n
                write(88,*) 100+400/(z.n-1)*(i-1), ans(i,   4)
                write(80,*) 100+400/(z.n-1)*(i-1), ans(i,  20)
                write(70,*) 100+400/(z.n-1)*(i-1), ans(i,  40)
                write(60,*) 100+400/(z.n-1)*(i-1), ans(i,  60)
                write(50,*) 100+400/(z.n-1)*(i-1), ans(i,  80) 
                write(40,*) 100+400/(z.n-1)*(i-1), ans(i, 100)
                write(30,*) 100+400/(z.n-1)*(i-1), ans(i, 120)
                write(20,*) 100+400/(z.n-1)*(i-1), ans(i, 140)
                write(10,*) 100+400/(z.n-1)*(i-1), ans(i, 160)
                write(5, *) 100+400/(z.n-1)*(i-1), ans(i, 170)
                write(2, *) 100+400/(z.n-1)*(i-1), ans(i, 176)
                write(0, *) 100+400/(z.n-1)*(i-1), ans(i, 180)
            end do
            close(88)
            close(80)
            close(70)
            close(60)
            close(50)
            close(40)
            close(30)
            close(20)
            close(10)
            close(5 )
            close(2 )
            close(0 )

            !maximum altitude and maximum concentration
            open(unit=900, name='n_max.txt')
            open(unit=901, name='h_max.txt')
            do j = 1, Nphi
                n_max = 1
                h_max = 1
                do i = 1, z.n
                    if (ans(i, j) .ge. n_max) then
                        n_max = ans(i, j)
                        h_max = z.d(i)
                    end if
                end do
                write(900, *) (j-5E-1)*180/(Nphi)-90, n_max
                write(901, *) (j-5E-1)*180/(Nphi)-90, h_max/1000
            end do
            close(900)
            close(901)
        end if

        if(convergence_profile_output .eq. 1) then ! .and. t*tau .ge. 86400*2
            C_norm = 0
            do j = 1, Nphi
                if(C_norm .le. ans(1, j)) then
                    C_norm = ans(1, j)
                end if
                do i = 1, z.n
                    if(C_norm .le. ans(i, j)) then
                        C_norm = ans(i, j)
                    end if
                end do
            end do
            write (1030, *) t, C_norm
            C_norm = ans(1, 30)
            do i = 1, z.n
                if(C_norm .le. ans(i, 30)) then
                    C_norm = ans(i, 30)
                end if
            end do
            write (1040, *) t*tau/86400, C_norm
        end if


    end if

    if(continuous_output .eq. 1 .and. diurnal_on .eq. 0) then
        write (900, '(1000(e10.3))') (ans(i, latitude_node), i = 1, z.n)
    end if

    if(continuous_output .eq. 1 .and. mod(t, 10) .eq. 0 .and. diurnal_on .eq. 1) then
            do i = 1, z.n
                write(900,*) t*tau,  100+(Hmax - 100)/(z.n-1)*(i-1), ans(i, 90)
                write(930,*) t*tau,  100+(Hmax - 100)/(z.n-1)*(i-1), ans(i, 60)
                write(960,*) t*tau,  100+(Hmax - 100)/(z.n-1)*(i-1), ans(i, 30)
            end do
            write (900, *)
            write (930, *)
            write (960, *)            
    end if
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

real*8 function f_m(z, phi)
    implicit none
    real*8 phi, z
    f_m =  3E+6 * ((z/1000)-100)/133 * exp((100-(z/1000))/133) * cos(phi/2)**2
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

