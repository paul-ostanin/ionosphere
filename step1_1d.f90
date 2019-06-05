program mat

use tridi_matrix
use vector

implicit none

type (tridiagonal_matrix) A, S
type (vect) b, z, h, hmid, nO, nO2, nN2, k, k1, p, pcurr, m, njold, njnew, n_zeroflow, dn, dn1, delta, D, D_node, u, Tn, Ti, Te, Tr, Tp, gradTp, nday, tau0, nnew(721), nold(721)
integer Nz, i, j, q, Te0, Tn0, Ti0, day, nonlinear_scheme_type, diurnal_on, Nphi, ub_type
real (8) dn1xdn1, n_k1_comp, sum, dn1xk1, k1xk1, dn1_orth_k1, int_k, int_n, D_analytical, Hmax, tau, h0, Fub, delta_norm, eps, tgdelta, sindelta, cosdelta, dphi, phi, coschi, pi, omega, sigma_O2, sigma_N2, sigma_O, sI, cI, R, u_phi, u_phi_N, u_phi_Nm1, u_phi_mh, u_phi_ph, u_phi_p1, u_phi_m1, u_phi_Np1

!opening files for writing the output
open(unit=10, name='res_gnp.txt')



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
!number of nodes in z
Nz = 81
!maximum altitude in km
Hmax = 500
!upper boundary condition: 1 - 3-rd order; 2 - Neumann; 3 - Dirichlet;
ub_type = 1


!Time step (in seconds) 5 min
tau = 1

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
    do i = 1, z.n+1
        gradTp.d(i)     =  0.5 * ((Te0 - Te.d(i)) * 9.8 / (287 * Te0) + (Ti0 - Ti.d(i)) * 9.8 / (287 * Ti0) ) * 1E-2
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

    call p.init(z.n+1)
    call k.init(z.n+1)
    do i = 1, z.n+1
            p.d(i) = ( 4E-7 * nO.d(i) )
            k.d(i) = ( 1.2E-12 * nN2.d(i) + 2.1E-11 * nO2.d(i) )
    end do

    !Diffusion coefficients vector. D.d(i) = D_{i+1/2}
    call D.init(z.n)
    call D_node.init(z.n+1)
    do i = 1, z.n
        D.d(i) = D_analytical(z.d(i)+h0/200)
        D_node.d(i) = D_analytical(z.d(i))
    end do
        D_node.d(z.n+1) = D_analytical(z.d(z.n) + h0/100)



    !Effective velocity vector
    !u.d(i) = u_i = D_i * (dTp/dz + mg/2k)/Tp, mg/2k ~ 5.6*10^{-5} [K/cm]
    call u.init(z.n+1)
    do i = 1, z.n+1    
        u.d(i) = ( 3E+17 / (nO.d(i) * sqrt(Tr.d(i))) * (56E-6 + gradTp.d(i)) )
    end do


!System matrix for the poles boundary (phi = +- 90). n(first) and n(last) will be known from boundary conditions.
call S.init(z.n)
!lower boundary condition: n_1 = P_1/k_1
S.d(1, 2) = 1


if(ub_type .eq. 1) then
    !upper boundary condition (3-rd order condition):
    S.d(z.n, 1) =    - D.d(z.n-1)*tau/(h0**2) + 0.5 * u.d(z.n-1)*tau/h0
    S.d(z.n, 2) = +1 + D.d(z.n-1)*tau/(h0**2) + 0.5 * u.d( z.n )*tau/h0
else if(ub_type .eq. 2) then
    !upper boundary condition (Neumann condition):
    S.d(z.n, 1) =    - D.d(z.n-1)*tau/(h0**2) + 0.5 * u.d(z.n-1)*tau/h0
    S.d(z.n, 2) = +1 + D.d(z.n-1)*tau/(h0**2) - 0.5 * u.d(z.n+1)*tau/h0
else if(ub_type .eq. 3) then
    !upper boundary condition (Dirichlet condition):
    S.d(z.n, 1) =    - D.d(z.n-1)*tau/(h0**2)  + 0.5 * u.d(z.n-1)*tau/h0
    S.d(z.n, 2) = +1 + (D.d(z.n-1)+D.d(z.n))*tau/(h0**2) 
end if

!middle matrix rows are similar
do i = 2, z.n - 1

! symmetric scheme
    S.d(i, 1) = -D.d(i-1)*tau/(h0**2)  + u.d(i-1)*tau/(2*h0)
    S.d(i, 2) = 1 + k.d(i)*tau + (D.d(i-1) + D.d(i))/(h0**2) * tau
    S.d(i, 3) = -D.d( i )*tau/(h0**2)  - u.d(i+1)*tau/(2*h0)

end do

do i = 1, z.n
!    print *, abs(S.d(i, 2)) - abs(S.d(i, 1)) - abs(S.d(i, 3))
end do 

!call S.print()

call njold.init(z.n)
call njnew.init(z.n)
call n_zeroflow.init(z.n)
call k1.init(z.n)
call dn1.init(z.n)
call dn.init(z.n)
call delta.init(z.n)
call njold.gen() !njold = (1, 1, ..., 1)
call nday.init(z.n)
call b.init(z.n) !generating RHS for S * njnew = b
call pcurr.init(z.n)

Fub = 0
delta_norm = 1


!calculating the solution at the north and the south poles: phi = +- 90
    do while (delta_norm > 1E-5)
        !print *, delta_norm

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
    n_zeroflow = njold

    do i = 1, z.n
        write(10,*) 100+(Hmax - 100)/(z.n-1)*(i-1), nday.d(i)
    end do


!     Fub = 1000000000
!     delta_norm = 1


! !calculating the solution at the north and the south poles: phi = +- 90
!     do while (delta_norm > 1E-5)
!         !print *, delta_norm

!         !Setting RHS in the middle
!             do i = 2, z.n - 1
!                     b.d(i) = njold.d(i) + tau * p.d(i)
!             end do

!         !Setting boundary conditions for the RHS
!         b.d(z.n) = tau/h.d(z.n-1) * Fub + njold.d(z.n)
!         b.d(1) = P.d(1)/ k.d(1)

!         !Solving the system
!             njnew = tridiagonal_matrix_algorithm(S, b)

!         delta = njnew - njold
!         delta_norm = delta.norm()
!         njold = njnew

!     end do

!     dn = njnew - n_zeroflow

!     int_k = 0
!     int_n = 0
!     do i = 1, z.n
!         int_k = int_k + h0* k.d(i)
!         int_n = int_n + h0*dn.d(i)
!     enddo
!     int_k = int_k / (z.n*h0)
!     int_n = int_n / (z.n*h0)
!     print*, int_k, int_n

!     do i = 1, z.n
!          k1.d(i) =  k.d(i) - int_k
!         dn1.d(i) = dn.d(i) - int_n
!     enddo

!     do i = 1, z.n
!         write(10,*) 100+(Hmax - 100)/(z.n-1)*(i-1), dn1.d(i)
!     end do

!     sum = 0
!     do i = 1, z.n
!         sum = sum + h0*dn1.d(i)
!     enddo
!     print*, sum


!     dn1xk1      = 0
!     dn1xdn1     = 0
!     k1xk1       = 0
!     dn1_orth_k1 = 0
!     sum         = 0
!     n_k1_comp   = 0

!     do i = 1, z.n
!         dn1xk1  = dn1xk1  + dn1.d(i)*k1.d(i)
!         dn1xdn1 = dn1xdn1 + dn1.d(i)*dn1.d(i)
!         k1xk1   = k1xk1   + k1.d(i)*k1.d(i)
!     enddo

!     do i = 1, z.n
!         dn1_orth_k1 = dn1_orth_k1 + (dn1.d(i) - k1.d(i)*dn1xk1/k1xk1)**2
!         n_k1_comp   = n_k1_comp + k.d(i)**2 * dn1xk1**2/(k1xk1**2)
!     enddo

!     ! do i = 1, z.n
!     !     dn1xk1 = dn1xk1 + k1.d(i)*dn1.d(i)
!     !     k1xk1  = k1xk1  + k1.d(i)*k1.d(i)
!     ! end do

!     ! do i = 1, z.n
!     !     dn1_orth_k1 = dn1_orth_k1 + (dn1.d(i)-k1.d(i)*dn1xk1/k1xk1)**2
!     !     sum = sum + dn1.d(i)**2
!     !     n_k1_comp = n_k1_comp + (k.d(i)*dn1xk1/k1xk1)**2
!     ! end do

!     print*, n_k1_comp, dn1_orth_k1, dn1xdn1
!     print*, dn1xdn1 - n_k1_comp - dn1_orth_k1

!     print*, n_k1_comp/dn1xdn1









 close(unit=10)

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