program mat

use tridi_matrix
use vector

implicit none

type (tridiagonal_matrix) A, S
type (vect) b, z, h, hmid, nO, nO2, nN2, k, p, pcurr, m, njold, njnew, delta, D, D_node, u, Tn, Ti, Te, Tr, Tp, gradTp, nday, tau0, nnew(721), nold(721)
integer i, j, q, Te0, Tn0, Ti0, day, diurnal_on, Nphi, mixed_switch
real (8) tau, h0, Fub, delta_norm, eps, tgdelta, sindelta, cosdelta, dphi, phi, coschi, pi, omega, sigma_O2, sigma_N2, sigma_O, sI, cI, R, u_phi, u_phi_N, u_phi_Nm1, u_phi_mh, u_phi_ph, u_phi_p1, u_phi_m1, u_phi_Np1

!opening files for writing the output
open(unit=10, name='res_gnp.txt')
!open(unit=105, name='new_flux_scheme_Nphi=180_tau=1_day05.txt')
!open(unit=11, name='new_flux_scheme_Nphi=180_tau=1_day1.txt')



pi = 3.141592653589793238462643

!number of nodes in phi
Nphi = 181
!latitude
dphi = pi / (Nphi-1)
!angle velocity of the Earth
omega = 2*pi/24/60/60
!magnetic inclination sin I
sI = 1
!Earth radius
R = 637100000

mixed_switch = 0



!Time step (in seconds) 5 min
tau = 50

!Vector of altitudes. Step h_i = z(i) - z(i - 1). Counting from 100 km to 500 km. z.d(i) is in metres.
call z.init(81)

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

!P = P_0 exp(tau_0(z)*(1-sec chi))
!tau_0 = sum_{N2, O2, O} sigma_i R_0*T_n/(M_i g) n_i(z)
!sigma is in cm^2
sigma_O2 = 2E-17
sigma_N2 = 15E-18
sigma_O  = 1E-17
call tau0.init(z.n)

do i = 1, z.n
    tau0.d(i)  = sigma_N2 * (8.31 * 100 * Tn.d(i))/(28E-3 * 9.8) * nN2.d(i) + &
             sigma_O2 * (8.31 * 100 * Tn.d(i))/(32E-3 * 9.8) * nO2.d(i) + &
             sigma_O  * (8.31 * 100 * Tn.d(i))/(16E-3 * 9.8) * nO.d(i)
end do

call p.init(z.n)
call k.init(z.n)
do i = 1, z.n
    p.d(i) = 4E-7 * nO.d(i)
    k.d(i) = 1.2E-12 * nN2.d(i) + 2.1E-11 * nO2.d(i)
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
    u.d(i) = 3E+17 / (nO.d(i) * sqrt(Tr.d(i))) * (56E-6 + gradTp.d(i))
end do

!System matrix for the poles boundary (phi = +- 90). n(first) and n(last) will be known from boundary conditions.
call S.init(z.n)
!lower boundary condition: n_1 = P_1/k_1
S.d(1, 2) = 1

!upper boundary condition:
S.d(z.n, 1) =    - D.d(z.n-1)*tau/(h.d(z.n-1)**2) + 0.5 * u.d(z.n-1)*tau/h.d(z.n-1)
S.d(z.n, 2) = +1 + D.d(z.n-1)*tau/(h.d(z.n-1)**2) + 0.5 * u.d( z.n )*tau/h.d(z.n-1)

!middle matrix rows are similar
do i = 2, z.n - 1

! symmetric scheme
    S.d(i, 1) = -D.d(i-1)*tau/(hmid.d(i) * h.d(i-1)) + u.d(i-1)*tau/(h.d(i) + h.d(i-1))
    S.d(i, 2) = 1 + k.d(i)*tau + (D.d(i-1)/h.d(i-1) + D.d(i)/h.d(i)) * tau / hmid.d(i)
    S.d(i, 3) = -D.d( i )*tau/(hmid.d(i) * h.d( i )) - u.d(i+1)*tau/(h.d(i) + h.d(i-1))

end do

call njold.init(z.n)
call njnew.init(z.n)
call delta.init(z.n)
call njold.gen() !njold = (1, 1, ..., 1)
call nday.init(z.n)
call b.init(z.n) !generating RHS for S * njnew = b
call pcurr.init(z.n)

Fub = 0
delta_norm = 1


!calculating the solution at the north and the south poles: phi = +- 90
    do while (delta_norm > 1E-5)

        !Setting RHS in the middle
            do i = 2, z.n - 1
                    b.d(i) = njold.d(i) + tau * p.d(i)
            end do

        !Setting boundary conditions for the RHS
        b.d(z.n) = tau/h.d(z.n-1) * Fub + njold.d(z.n)
        b.d(1) = P.d(1)/ k.d(1)

        !Solving the system
            njnew = tridiagonal_matrix_algorithm(S, b)

        delta = njnew - njold
        delta_norm = delta.norm()
        njold = njnew

    end do

    nday = njold

do j = 1, Nphi
    call nnew(j).init(z.n)
    call nold(j).init(z.n)
    call nold(j).gen()
end do

    nold( 1) = nday
    nold(Nphi) = nday

diurnal_on = 0

b.d(1) = p.d(1)/k.d(1)

do j = 0, 86400/tau*10
!print *, j
if(mod(j, 100) .eq. 0) then
    print *, j
end if



do q = 2, Nphi-1
! angles phi from -90 to 90; conditions in -90 and 90 are set


    !sinus and cosinus of magnetic inclination angle I
    sI = sin(atan(2*tan(-pi/2+(q-1)*dphi)))
    cI = cos(atan(2*tan(-pi/2+(q-1)*dphi)))

    !lower boundary condition: n_1 = P_1/k_1
    S.d(1, 2) = 1
    !upper boundary condition:
    if(sI*cI .le. 0) then

        S.d(z.n, 1) =  (-D.d(z.n-1)*tau/(h.d(z.n-1)**2) + 0.5 * u.d(z.n-1)*tau/h.d(z.n-1)) * sI**2 - &
                0.5*D_node.d(z.n-1)*tau*sI*cI/(R*h0*dphi)*mixed_switch
        S.d(z.n, 2) = +1 + (D.d(z.n-1)*tau/(h.d(z.n-1)**2) + 0.5 * u.d( z.n )*tau/h.d(z.n-1)) * sI**2 + &
                0.5*D_node.d( z.n )*tau*sI*cI/(R*h0*dphi)*mixed_switch

        b.d(z.n) = +tau/h.d(z.n-1) * Fub + nold(q).d(z.n) - &
                0.5*tau*sI*cI/(2*R*h0*dphi)*(-D_node.d(z.n)*nold(q-1).d(z.n) + D_node.d(z.n-1)*nold(q-1).d(z.n-1))*mixed_switch


    else

        S.d(z.n, 1) =  (-D.d(i-1)*tau/(h.d(z.n-1)**2) + 0.5 * u.d(z.n-1)*tau/h.d(z.n-1)) * sI**2 + &
                0.5*D_node.d(z.n-1)*tau*sI*cI/(R*h0*dphi)*mixed_switch
        S.d(z.n, 2) = 1 + (D.d(z.n-1)*tau/(h.d(z.n-1)**2) + 0.5 * u.d( z.n )*tau/h.d(z.n-1)) * sI**2 - &
                0.5*D_node.d( z.n )*tau*sI*cI/(R*h0*dphi)*mixed_switch

        b.d(z.n) = +tau/h.d(z.n-1) * Fub + nold(q).d(z.n) - &
                0.5*tau*sI*cI/(2*R*h0*dphi)*(D_node.d(z.n)*nold(q+1).d(z.n) - D_node.d(z.n-1)*nold(q+1).d(z.n-1))*mixed_switch
    end if

    do i = 2, z.n - 1
    if(sI*cI .le. 0) then

        S.d(i, 1) =  (-D.d(i-1)*tau/(hmid.d(i) * h.d(i-1)) + u.d(i-1)*tau/(h.d(i) + h.d(i-1))) * sI**2 - &
                0.5*D_node.d(i-1)*tau*sI*cI/(R*h0*dphi)*mixed_switch
        S.d(i, 2) = 1 + k.d(i)*tau + (D.d(i-1)/h.d(i-1) + D.d(i)/h.d(i)) * tau / hmid.d(i) * sI**2 + &
                    D_node.d(i)*tau*sI*cI/(R*h0*dphi)*mixed_switch
        S.d(i, 3) =  (-D.d( i )*tau/(hmid.d(i) * h.d( i )) - u.d(i+1)*tau/(h.d(i) + h.d(i-1))) * sI**2 - &
                0.5*D_node.d(i+1)*tau*sI*cI/(R*h0*dphi)*mixed_switch

        b.d(i) = nold(q).d(i) + tau * p.d(i) - &
                0.5*tau*sI*cI/(2*R*h0*dphi)*(-D_node.d(i)*(nold(q-1).d(i)+nold(q+1).d(i)) + &
                                         D_node.d(i-1)*nold(q-1).d(i-1) + &
                                         D_node.d(i+1)*nold(q+1).d(i+1))*mixed_switch
    else

        S.d(i, 1) =  (-D.d(i-1)*tau/(hmid.d(i) * h.d(i-1)) + u.d(i-1)*tau/(h.d(i) + h.d(i-1))) * sI**2 + &
                0.5*D_node.d(i-1)*tau*sI*cI/(R*h0*dphi)*mixed_switch
        S.d(i, 2) = 1 + k.d(i)*tau + (D.d(i-1)/h.d(i-1) + D.d(i)/h.d(i)) * tau / hmid.d(i) * sI**2 - &
                    D_node.d(i)*tau*sI*cI/(R*h0*dphi)*mixed_switch
        S.d(i, 3) =  (-D.d( i )*tau/(hmid.d(i) * h.d( i )) - u.d(i+1)*tau/(h.d(i) + h.d(i-1))) * sI**2 + &
                0.5*D_node.d(i+1)*tau*sI*cI/(R*h0*dphi)*mixed_switch

        b.d(i) = nold(q).d(i) + tau * p.d(i) - &
                0.5*tau*sI*cI/(2*R*h0*dphi)*(D_node.d(i)*(nold(q-1).d(i)+nold(q+1).d(i)) - &
                                         D_node.d(i+1)*nold(q-1).d(i+1) - &
                                         D_node.d(i-1)*nold(q+1).d(i-1))*mixed_switch
    end if

    end do

    nnew(q) = tridiagonal_matrix_algorithm(S, b)


        if(j*tau .eq. 86400*10) then
        do i = 1, z.n
            write(10,*) (q)*180/(Nphi-1)-90, 100+400/(z.n-1)*(i-1), nnew(q).d(i)
        end do
        write (10, *)
        end if


end do

    do i = 2, Nphi-1
        nold(i) = nnew(i)
    end do

end do


 close(unit=10)
 close(unit=11)
do j = 1, Nphi
 call nnew(j).destroy()
 call nold(j).destroy()
end do
call b.destroy()
call z.destroy()
call h.destroy()
call m.destroy()
call hmid.destroy()
call nO.destroy()
call nO2.destroy()
call nN2.destroy()
call Te.destroy()
call Tn.destroy()
call Ti.destroy()
call Tr.destroy()
call Tp.destroy()
call gradTp.destroy()
call p.destroy()
call k.destroy()
call D.destroy()
call u.destroy()
call njold.destroy()
call delta.destroy()
call njnew.destroy()
call S.destroy()
end