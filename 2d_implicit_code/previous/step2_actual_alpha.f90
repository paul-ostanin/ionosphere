program mat

use tridi_matrix
use vector

implicit none

type (tridiagonal_matrix) S_y, S_z
type (vect) rhs_z, rhs_y, z, h, hmid, nO, nO2, nN2, k, p, pcurr, m, n_old, n_new, delta, D, D_node, u, Tn, Ti, Te, Tr, Tp, gradTp, n_day, tau0, n_new_z(1441), n_old_z(1441), n_new_y(401), n_old_y(401), n_new_z_1(1441), n_old_z_1(1441), n_new_y_1(401), n_old_y_1(401), error(1441)
integer s_i, s_j, i, j, t, Te0, Tn0, Ti0, day, nonlinear_scheme_type, profile_output, diurnal_on, convergence_test, pk_switch, mixed_z_switch, mixed_y_switch, transf_yz_switch, transf_y_switch, second_step_scheme_type, upper_bound_type, monotonizator
real (8) tau, tau_1, h0, F_z, delta_norm, eps, tgdelta, sindelta, cosdelta, dphi, phi, coschi, pi, omega, sigma_O2, sigma_N2, sigma_O, sI, cI, R, u_phi, u_phi_1, u_phi_2, u_phi_3, u_phi_4, u_z, u_z_mh, u_z_ph, u_z_m1, u_z_p1, x, A, B, u_phi_ph, u_phi_mh, Ndays, Niter, delta_err, sIph, cIph, sImh, cImh, l1_norm, l1_norm_err


integer, parameter :: max_size = 1000000
integer, parameter :: max_nonzero = 10000000
integer, parameter :: maxwr = max_nonzero + 8 * max_size
integer, parameter :: maxwi = max_nonzero + 2 * max_size + 1
integer, parameter :: Nz = 201
integer, parameter :: Nphi = 180


integer ia(max_size + 1), ja(max_nonzero), iw(maxwi)
double precision arr(max_nonzero), f(max_size), n(max_size), prev_n(max_size), rw(maxwr), v(max_size)

external matvec, prevec0
integer ITER, INFO, NUNIT, ierr, ipalu, ipjlu, ipju, ipiw
real*8 RESID, ddot

external ddot, dcopy

integer imatvec(1), iprevec(1), ipbcg, matrix_size, nonzero
real*8 resinit


integer counter_a, counter_ia, counter_rhs
real*8 alpha1, alpha2, alpha3, alpha4, alpha5, beta1
real*8 currt
real*8 C_norm, L2_norm

real*8 prev, analytical_solution, cubature_formula

real*8 stencil(Nz, Nphi, -1:1, -1:1), operator(-1:1, -1:1), diffusion_transfer_z(-1:1, -1:1), diffusion_transfer_y(-1:1, -1:1), mixed_z(-1:1, -1:1), mixed_y(-1:1, -1:1)
real*8 rhs_z_phi(Nz, Nphi), ans(Nz, Nphi)
real*8 o

!opening files for writing the output
!    open(unit=1, name='res.txt')
    open(unit=12, name='res_gnp.txt')

pi = 3.141592653589793238462643

!nonlinear_scheme_type variable switches the u_phi-approximation. 
!nonlinear_scheme_type = 8
 

!latitude
dphi = pi / Nphi
!angle velocity of the Earth 
omega = 2*pi/24/60/60
!magnetic inclination sin I
sI = 1
!Earth radius
R = 637100000
!number of calculation days
Ndays = 0.25
Niter = 800
!upper boundary electron flow
F_z = 0

!Time step (in seconds) 5 min
tau = 100

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
upper_bound_type = 4
!switcher for the monotonizator after both steps - all less or equal than zero elements are made 1
monotonizator = 1
diurnal_on = 0
convergence_test = 0
profile_output = 1

!Initialization block
    !Vector of altitudes. Step h_i = z(i) - z(i - 1). Counting from 100 km to 500 km. z.d(i) is in metres.
    call z.init(Nz)

    !Space step (in cm) 5 km
    h0 = 400E+5 / (z.n - 1)
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

    call Ti.init(z.n)
    call Tn.init(z.n)
    call Te.init(z.n)
    call Tr.init(z.n)
    call Tp.init(z.n)

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

    !dT/dz = s(T_\infty - T_0)exp(-s(z-z0)) = (T_\infty - T)*s - for Te, Tn, Ti. This derivative is in [K/cm].
    call gradTp.init(z.n)
    do i = 1, z.n
        gradTp.d(i) = 0.5 * ((Te0 - Te.d(i)) * 9.8 / (287 * Te0) + (Ti0 - Ti.d(i)) * 9.8 / (287 * Ti0) ) * 1E-2
    end do

    call nO.init(z.n)
    call nO2.init(z.n)
    call nN2.init(z.n)

    do i = 1, z.n
        nO.d(i)  = 2.8E+10 * exp(-9.8 * 16E-3 / (8.31 * Tn.d(i)) * (z.d(i) - 140000))
        nO2.d(i) = 5.6E+9  * exp(-9.8 * 32E-3 / (8.31 * Tn.d(i)) * (z.d(i) - 140000))
        nN2.d(i) = 5.2E+10 * exp(-9.8 * 28E-3 / (8.31 * Tn.d(i)) * (z.d(i) - 140000))
    end do


    call p.init(z.n)
    call k.init(z.n)
    do i = 1, z.n
            p.d(i) = ( 4E-7 * nO.d(i) ) * pk_switch
            k.d(i) = ( 1.2E-12 * nN2.d(i) + 2.1E-11 * nO2.d(i) ) * pk_switch
    end do


    !Diffusion coefficients vector. D.d(i) = D_{i+1/2}
    call D.init(z.n)
    do i = 1, z.n - 1
        D.d(i) = 3E+17 * Tp.interp(z, m.d(i)) / (nO.interp(z, m.d(i)) * sqrt(Tr.interp(z, m.d(i))))
    end do
        D.d(z.n) = D.d(z.n - 1) + (D.d(z.n - 1) - D.d(z.n - 2)) !extrapolation

    call D_node.init(z.n)
    do i = 1, z.n
        D_node.d(i) = 3E+17 * Tp.d(i) / (nO.d(i) * sqrt(Tr.d(i)))
    end do

    !Effective velocity vector
    !u.d(i) = u_i = D_i * (dTp/dz + mg/2k)/Tp, mg/2k ~ 5.6*10^{-5} [K/cm]
    call u.init(z.n)
    do i = 1, z.n    
        u.d(i) = ( 3E+17 / (nO.d(i) * sqrt(Tr.d(i))) * (56E-6 + gradTp.d(i)) ) * transf_yz_switch
    end do


!initialization of n before iterations
    n = 1d0
    prev_n = 0d0

do t = 0, 1000
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
                diffusion_transfer_z(+1, :) = [0d0, (-D.d( i )*tau/(hmid.d(i)*h.d( i )) - u.d(i+1)*tau/(h.d(i)+h.d(i-1)))*sI**2, 0d0]
                diffusion_transfer_z( 0, :) = [0d0, 1 + k.d(i)*tau + (D.d(i-1)/h.d(i-1) + D.d(i)/h.d(i))*tau/hmid.d(i)*sI**2   , 0d0]
                diffusion_transfer_z(-1, :) = [0d0, (-D.d(i-1)*tau/(hmid.d(i)*h.d(i-1)) + u.d(i-1)*tau/(h.d(i)+h.d(i-1)))*sI**2, 0d0]


                diffusion_transfer_y(+1, :) = [0d0, 0d0, 0d0]
                diffusion_transfer_y( 0, :) = [0d0, & 
                                               1+(D_node.d(i)/R*A(phi+dphi/2)/(dphi**2) - u.d(i)/2*B(phi-dphi)/(2*dphi))*tau/(R)/cos(phi), &
                                                (-D_node.d(i)/R*A(phi+dphi/2)/(dphi**2) + u.d(i)/2*B(phi+dphi)/(2*dphi))*tau/(R)/cos(phi)  ]
                diffusion_transfer_y(-1, :) = [0d0, 0d0, 0d0]

                mixed_z(+1, :) = [0d0,  0d0                                        ,  0d0                                          ]
                mixed_z( 0, :) = [0d0, -0.5*D_node.d( i )*tau*sIph*cIph/(R*h0*dphi),  0.5*tau/(2*R*h0*dphi)*D_node.d( i )*sIph*cIph]
                mixed_z(-1, :) = [0d0,  0.5*D_node.d(i-1)*tau*sIph*cIph/(R*h0*dphi), -0.5*tau/(2*R*h0*dphi)*D_node.d(i-1)*sIph*cIph]

                mixed_y(+1, :) = [0d0,  0d0                                           ,  0d0                                                ]
                mixed_y( 0, :) = [0d0, -0.5*D.d(i-1)/2*B(phi)*tau/(h0*dphi*R*cos(phi)),  0.5*D.d(i-1)/2*B(phi+dphi)*tau/(h0*dphi*R*cos(phi))]
                mixed_y(-1, :) = [0d0,  0.5*D.d(i-1)/2*B(phi)*tau/(h0*dphi*R*cos(phi)), -0.5*D.d(i-1)/2*B(phi+dphi)*tau/(h0*dphi*R*cos(phi))]


                do s_i = -1, 1
                    do s_j = -1, 1
                        operator(s_i, s_j) = diffusion_transfer_z(s_i, s_j) + diffusion_transfer_y(s_i, s_j) + mixed_z(s_i, s_j) + mixed_y(s_i, s_j)
                    enddo
                enddo

            else if (j == Nphi .and. i /= Nz) then

                !north pole
                diffusion_transfer_z(+1, :) = [0d0, (-D.d( i )*tau/(hmid.d(i)*h.d( i )) - u.d(i+1)*tau/(h.d(i)+h.d(i-1)))*sI**2, 0d0]
                diffusion_transfer_z( 0, :) = [0d0, 1 + k.d(i)*tau + (D.d(i-1)/h.d(i-1) + D.d(i)/h.d(i))*tau/hmid.d(i)*sI**2   , 0d0]
                diffusion_transfer_z(-1, :) = [0d0, (-D.d(i-1)*tau/(hmid.d(i)*h.d(i-1)) + u.d(i-1)*tau/(h.d(i)+h.d(i-1)))*sI**2, 0d0]

                diffusion_transfer_y(+1, :) = [0d0, 0d0, 0d0]
                diffusion_transfer_y( 0, :) = [ (-D_node.d(i)/R*A(phi-dphi/2)/(dphi**2) - u.d(i)/2*B(phi-dphi)/(2*dphi))*tau/(R)/cos(phi), &
                                               1+(D_node.d(i)/R*A(phi-dphi/2)/(dphi**2) + u.d(i)/2*B(phi+dphi)/(2*dphi))*tau/(R)/cos(phi), &
                                               0d0]
                diffusion_transfer_y(-1, :) = [0d0, 0d0, 0d0]

                mixed_z(+1, :) = [ 0d0                                          ,  0d0                                        , 0d0]
                mixed_z( 0, :) = [-0.5*tau/(2*R*h0*dphi)*D_node.d( i )*sImh*cImh,  0.5*D_node.d( i )*tau*sImh*cImh/(R*h0*dphi), 0d0]
                mixed_z(-1, :) = [ 0.5*tau/(2*R*h0*dphi)*D_node.d(i-1)*sImh*cImh, -0.5*D_node.d(i-1)*tau*sImh*cImh/(R*h0*dphi), 0d0]

                mixed_y(+1, :) = [ 0d0                                                ,  0d0                                           ,  0d0]
                mixed_y( 0, :) = [-0.5*D.d(i-1)/2*B(phi-dphi)*tau/(h0*dphi*R*cos(phi)),  0.5*D.d(i-1)/2*B(phi)*tau/(h0*dphi*R*cos(phi)),  0d0]
                mixed_y(-1, :) = [ 0.5*D.d(i-1)/2*B(phi-dphi)*tau/(h0*dphi*R*cos(phi)), -0.5*D.d(i-1)/2*B(phi)*tau/(h0*dphi*R*cos(phi)),  0d0]

                do s_i = -1, 1
                    do s_j = -1, 1
                        operator(s_i, s_j) = diffusion_transfer_z(s_i, s_j) + diffusion_transfer_y(s_i, s_j) + mixed_z(s_i, s_j) + mixed_y(s_i, s_j)
                    enddo
                enddo

!            else if ((j == 1 .or. j == Nphi) .and. i == Nz) then

!                !left and right upper corners
!                diffusion_transfer_z(+1, :) = [0d0, 0d0, 0d0]
!                diffusion_transfer_z( 0, :) = [0d0, 1 + (D.d(z.n-1)*tau/(h.d(z.n-1)**2) + 0.5 * u.d( z.n )*tau/h.d(z.n-1)) * sI**2, 0d0]
!                diffusion_transfer_z(-1, :) = [0d0,   (- D.d(z.n-1)*tau/(h.d(z.n-1)**2) + 0.5 * u.d(z.n-1)*tau/h.d(z.n-1)) * sI**2, 0d0]


            else if (i == Nz) then

                !upper boundary
                diffusion_transfer_z(+1, :) = [0d0, 0d0                                                                       , 0d0]
                diffusion_transfer_z( 0, :) = [0d0, 1 + (D.d(z.n-1)*tau/(h.d(z.n-1)**2) + 0.5*u.d( z.n )*tau/h.d(z.n-1))*sI**2, 0d0]
                diffusion_transfer_z(-1, :) = [0d0, (-D.d(z.n-1)*tau/(h.d(z.n-1)**2) + 0.5*u.d(z.n-1)*tau/h.d(z.n-1))*sI**2   , 0d0]

                if(sI*cI .ge. 0) then
                    mixed_z(+1, :) = [0d0,                                     0d0,                                             0d0]
                    mixed_z( 0, :) = [-0.5*tau*sI*cI/(2*R*h0*dphi)*D_node.d( z.n ),  0.5*D_node.d( z.n )*tau*sI*cI/(R*h0*dphi), 0d0]
                    mixed_z(-1, :) = [ 0.5*tau*sI*cI/(2*R*h0*dphi)*D_node.d(z.n-1), -0.5*D_node.d(z.n-1)*tau*sI*cI/(R*h0*dphi), 0d0]
                else
                    mixed_z(+1, :) = [0d0,                                        0d0,                                          0d0]
                    mixed_z( 0, :) = [0d0, -0.5*D_node.d( z.n )*tau*sI*cI/(R*h0*dphi),  0.5*tau*sI*cI/(2*R*h0*dphi)*D_node.d( z.n )]
                    mixed_z(-1, :) = [0d0,  0.5*D_node.d(z.n-1)*tau*sI*cI/(R*h0*dphi), -0.5*tau*sI*cI/(2*R*h0*dphi)*D_node.d(z.n-1)]
                end if
                do s_i = -1, 1
                    do s_j = -1, 1
                        operator(s_i, s_j) = diffusion_transfer_z(s_i, s_j) + mixed_z(s_i, s_j)
                    enddo
                enddo

            else

                !the main part of the operator

                diffusion_transfer_z(+1, :) = [0d0, (-D.d( i )*tau/(hmid.d(i)*h.d( i )) - u.d(i+1)*tau/(h.d(i)+h.d(i-1)))*sI**2, 0d0]
                diffusion_transfer_z( 0, :) = [0d0, 1 + k.d(i)*tau + (D.d(i-1)/h.d(i-1) + D.d(i)/h.d(i))*tau/hmid.d(i)*sI**2   , 0d0]
                diffusion_transfer_z(-1, :) = [0d0, (-D.d(i-1)*tau/(hmid.d(i)*h.d(i-1)) + u.d(i-1)*tau/(h.d(i)+h.d(i-1)))*sI**2, 0d0]

                diffusion_transfer_y(+1, :) = [0d0, 0d0, 0d0]
                diffusion_transfer_y( 0, :) = [(-D_node.d(i)/R*A(phi-dphi/2)/(dphi**2) + (-u.d(i)/2)*B(phi-dphi)/(2*dphi))*tau/(R)/cos(phi), &
                                                1+(D_node.d(i)/R *(A(phi-dphi/2) + A(phi+dphi/2))/(dphi**2))*tau/(R)/cos(phi), &
                                               (-D_node.d(i)/R*A(phi+dphi/2)/(dphi**2) - (-u.d(i)/2)*B(phi+dphi)/(2*dphi))*tau/(R)/cos(phi)]
                diffusion_transfer_y(-1, :) = [0d0, 0d0, 0d0]

                if(sI .ge. 0) then
                    mixed_z(+1, :) = [ 0d0,                                           -0.5*D_node.d(i+1)*tau*sIph*cIph/(R*h0*dphi)          ,  0.5*tau/(2*R*h0*dphi)*D_node.d(i+1)*sIph*cIph]
                    mixed_z( 0, :) = [-0.5*tau/(2*R*h0*dphi)*D_node.d( i )*sImh*cImh,  0.5*D_node.d(i)*tau*(sIph*cIph+sImh*cImh)/(R*h0*dphi), -0.5*tau/(2*R*h0*dphi)*D_node.d( i )*sIph*cIph]
                    mixed_z(-1, :) = [ 0.5*tau/(2*R*h0*dphi)*D_node.d(i-1)*sImh*cImh, -0.5*D_node.d(i-1)*tau*sImh*cImh/(R*h0*dphi)          ,  0d0                                          ]

                    mixed_y(+1, :) = [ 0d0                                                , -0.5*D.d( i )/2*B(phi)*tau/(h0*dphi*R*cos(phi))         ,  0.5*D.d( i )/2*B(phi+dphi)*tau/(h0*dphi*R*cos(phi))]
                    mixed_y( 0, :) = [-0.5*D.d(i-1)/2*B(phi-dphi)*tau/(h0*dphi*R*cos(phi)),  0.5*(D.d(i)+D.d(i-1))/2*B(phi)*tau/(h0*dphi*R*cos(phi)), -0.5*D.d( i )/2*B(phi+dphi)*tau/(h0*dphi*R*cos(phi))]
                    mixed_y(-1, :) = [ 0.5*D.d(i-1)/2*B(phi-dphi)*tau/(h0*dphi*R*cos(phi)), -0.5*D.d(i-1)/2*B(phi)*tau/(h0*dphi*R*cos(phi))         ,  0d0                                                ]
                else
                    mixed_z(+1, :) = [-0.5*tau/(2*R*h0*dphi)*D_node.d(i+1)*sImh*cImh,  0.5*D_node.d(i+1)*tau*sImh*cImh/(R*h0*dphi)          ,  0d0                                          ]
                    mixed_z( 0, :) = [ 0.5*tau/(2*R*h0*dphi)*D_node.d( i )*sImh*cImh, -0.5*D_node.d(i)*tau*(sImh*cImh+sIph*cIph)/(R*h0*dphi),  0.5*tau/(2*R*h0*dphi)*D_node.d( i )*sIph*cIph]
                    mixed_z(-1, :) = [ 0d0,                                            0.5*D_node.d(i-1)*tau*sIph*cIph/(R*h0*dphi)          , -0.5*tau/(2*R*h0*dphi)*D_node.d(i-1)*sIph*cIph]

                    mixed_y(+1, :) = [-0.5*D.d( i )/2*B(phi-dphi)*tau/(h0*dphi*R*cos(phi)),  0.5*D.d( i )/2*B(phi)*tau/(h0*dphi*R*cos(phi))         ,  0d0                                                ]
                    mixed_y( 0, :) = [ 0.5*D.d( i )/2*B(phi-dphi)*tau/(h0*dphi*R*cos(phi)), -0.5*(D.d(i)+D.d(i-1))/2*B(phi)*tau/(h0*dphi*R*cos(phi)),  0.5*D.d(i-1)/2*B(phi+dphi)*tau/(h0*dphi*R*cos(phi))]
                    mixed_y(-1, :) = [ 0d0,                                                  0.5*D.d(i-1)/2*B(phi)*tau/(h0*dphi*R*cos(phi)),          -0.5*D.d(i-1)/2*B(phi+dphi)*tau/(h0*dphi*R*cos(phi))]
                end if

                do s_i = -1, 1
                    do s_j = -1, 1
                        operator(s_i, s_j) = diffusion_transfer_z(s_i, s_j) + diffusion_transfer_y(s_i, s_j) + mixed_z(s_i, s_j) + mixed_y(s_i, s_j)
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
    do i = 1, Nz
        do j = 1, Nphi
            if (i == 1) then
                rhs_z_phi(i, j) = p.d(i)
            else if ((j == 1 .or. j == Nphi) .and. i == Nz) then
                rhs_z_phi(i, j) = tau/h.d(z.n-1) * F_z + ans(i, j)
            else if (i == Nz) then
                rhs_z_phi(i, j) = tau/h.d(z.n-1) * F_z + ans(i, j)
            else
                rhs_z_phi(i, j) = tau * p.d(i) + ans(i, j)
            end if
        enddo
    enddo

    !CSR matrix forming + RHS for the 2d system
    !currt = t*tau
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
        resinit = dsqrt(ddot(matrix_size, f, 1, f, 1))
        write(*, *) "resinit", resinit
        n = 0d0

        ITER = 1000
        RESID = 1d-30 * resinit
        INFO = 0
        NUNIT = 6 ! 6 to output
        iprevec(1) = matrix_size
        imatvec(1) = matrix_size

        call slpbcgs(prevec0, iprevec, iw, rw,   &
                     matvec, imatvec, ia, ja, arr, &
                     rw(ipbcg), matrix_size, 8, matrix_size, f, n,   &
                     ITER, RESID, INFO, NUNIT)
        print *, 'INFO =', INFO

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

       print *, 'Residual = ', resinit

end do


!Printing output
do j = 1, Nphi
    do i = 1, z.n
        write(12,*) (j-5E-1)*180/(Nphi)-90, 100+400/(z.n-1)*(i-1), ans(i, j)
    end do
    write (12, *)
end do






!destructor
!    close(unit=1)
    close(unit=12)
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
