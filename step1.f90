program mat

use tridi_matrix
use vector

implicit none

type (tridiagonal_matrix) A, S
type (vect) b, z, h, hmid, nO, nO2, nN2, k, p, pcurr, m, njold, njnew, delta, D, u, Tn, Ti, Te, Tr, Tp, gradTp, nday, tau0, nnew(721), nold(721)
integer i, j, q, Te0, Tn0, Ti0, day, nonlinear_scheme_type, diurnal_on, Nphi
real (8) tau, h0, Fub, delta_norm, eps, tgdelta, sindelta, cosdelta, dphi, phi, coschi, pi, omega, sigma_O2, sigma_N2, sigma_O, sI, cI, R, u_phi, u_phi_N, u_phi_Nm1, u_phi_mh, u_phi_ph, u_phi_p1, u_phi_m1, u_phi_Np1

!opening files for writing the output
open(unit=10, name='res.txt')
open(unit=105, name='new_flux_scheme_Nphi=180_tau=1_day05.txt')
open(unit=11, name='new_flux_scheme_Nphi=180_tau=1_day1.txt')
open(unit=12, name='new_flux_scheme_Nphi=180_tau=1_day2.txt')
open(unit=15, name='new_flux_scheme_Nphi=180_tau=1_day5.txt')
open(unit=110, name='new_flux_scheme_Nphi=180_tau=1_day10.txt')
open(unit=115, name='new_flux_scheme_Nphi=180_tau=1_day15.txt')
open(unit=120, name='new_flux_scheme_Nphi=180_tau=1_day20.txt')


pi = 3.141592653589793238462643

!nonlinear_scheme_type variable switches the u_phi-approximation. 
!If n_s_t = 1, u_phi = 1/n d(ln n)/dphi
!If n_s_t = 2, u_phi = (n(phi+dphi)-n(phi-dphi))/(n(phi+dphi)+n(phi-dphi)) * 1/dphi
!If n_s_t = 3, u_phi 1/n d(ln n)/dphi and the directed difference scheme is used in the equation approximation
!If n_s_t = 4, u_phi = 1/n d(ln n)/dphi, bnd_cond and equation are approximated with directed difference
nonlinear_scheme_type = 7

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



!Time step (in seconds) 5 min
tau = 80

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

!	print *
!	call S.print()
!	print *


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

do j = 0, 86400/tau*20
!print *, j
do q = 2, Nphi-1
! angles phi from -90 to 90; conditions in -90 and 90 are set

diurnal_on = 0 !switcher for the diurnal evolution mode. 0 corresponds to stationary solutions.

if(diurnal_on .eq. 1) then

	day = (j*tau)/86400 + 1/2 !starting from the middle of the 1-st day
	tgdelta = tan(pi/180*23.5) * sin(2*pi/365 * (day - 80))
	sindelta = tgdelta/sqrt(1+tgdelta**2)
	cosdelta = sqrt(1-sindelta**2) !cos of the zenith angle is > 0

        do i = 2, z.n - 1
		coschi = sin(-pi/2+q*dphi)*sindelta - cos(-pi/2+q*dphi)*cosdelta*cos(omega * (tau * j + 86400/2)) !start: middle of the 1 day
		if(coschi > 1E-6) then 
		        b.d(i) = nold(q).d(i) + tau * p.d(i) * exp(tau0.d(i) * (1-1/coschi))
        	else 
			b.d(i) = nold(q).d(i)
        	end if
	end do

	b.d(z.n) = tau/h.d(z.n-1) * Fub + nold(q).d(z.n)

	if(coschi > 1E-6) then 
	        b.d(1) = p.d(1)/k.d(1) * exp(tau0.d(1) * (1-1/coschi))
       	else 
		b.d(1) = 0
       	end if

else if(diurnal_on .eq. 0) then

        do i = 2, z.n - 1
		b.d(i) = nold(q).d(i) + tau * p.d(i)
	end do

	b.d(z.n) = +tau/h.d(z.n-1) * Fub + nold(q).d(z.n)

	b.d(1) = p.d(1)/k.d(1)

end if       


	!sinus and cosinus of magnetic inclination angle I
	sI = sin(atan(2*tan(-pi/2+(q-1)*dphi)))
	cI = cos(atan(2*tan(-pi/2+(q-1)*dphi)))

	!lower boundary condition: n_1 = P_1/k_1
	S.d(1, 2) = 1
	!upper boundary condition:
	if (nold(q+1).d(z.n) .eq. 0 .or. nold(q-1).d(z.n) .eq. 0 .or. nold(q+1).d(z.n-1) .eq. 0 .or. nold(q+1).d(z.n-1) .eq. 0) then
		S.d(z.n, 1) =    (- D.d(z.n-1)*tau/(h.d(z.n-1)**2) + 0.5 * u.d(z.n-1)*tau/h.d(z.n-1)) * sI**2
		S.d(z.n, 2) = +1 + (D.d(z.n-1)*tau/(h.d(z.n-1)**2) + 0.5 * u.d( z.n )*tau/h.d(z.n-1)) * sI**2
	else if (nonlinear_scheme_type .eq. 1 .or. nonlinear_scheme_type .eq. 3) then
		u_phi = -1/R * D.d(i) * sI * cI * log(nold(q+1).d(z.n-1)/nold(q-1).d(z.n-1))/(2*dphi)
		S.d(z.n, 1) =    (- D.d(z.n-1)*tau/(h.d(z.n-1)**2) + 0.5 * u.d(z.n-1)*tau/h.d(z.n-1)) * sI**2 + 0.5*u_phi*tau/h.d(z.n-1)
		u_phi = -1/R * D.d(i) * sI * cI * log(nold(q+1).d(z.n)/nold(q-1).d(z.n))/(2*dphi)
		S.d(z.n, 2) = +1 + (D.d(z.n-1)*tau/(h.d(z.n-1)**2) + 0.5 * u.d( z.n )*tau/h.d(z.n-1)) * sI**2 + 0.5*u_phi*tau/h.d(z.n-1)
	else if (nonlinear_scheme_type .eq. 2) then
		u_phi = -1/R * D.d(i) * sI * cI * (nold(q+1).d(z.n-1)-nold(q-1).d(z.n-1))/(nold(q+1).d(z.n-1)+nold(q-1).d(z.n-1))/dphi
		S.d(z.n, 1) =    (- D.d(z.n-1)*tau/(h.d(z.n-1)**2) + 0.5 * u.d(z.n-1)*tau/h.d(z.n-1)) * sI**2 + 0.5*u_phi*tau/h.d(z.n-1)
		u_phi = -1/R * D.d(i) * sI * cI * (nold(q+1).d(z.n)-nold(q-1).d(z.n))/(nold(q+1).d(i)+nold(q-1).d(i))/dphi
		S.d(z.n, 2) = +1 + (D.d(z.n-1)*tau/(h.d(z.n-1)**2) + 0.5 * u.d( z.n )*tau/h.d(z.n-1)) * sI**2 + 0.5*u_phi*tau/h.d(z.n-1)
	else if (nonlinear_scheme_type .eq. 4 .or. nonlinear_scheme_type .eq. 5) then
		u_phi = -1/R * D.d(i) * sI * cI * (nold(q+1).d(z.n-1)-nold(q-1).d(z.n-1))/(nold(q+1).d(z.n-1)+nold(q-1).d(z.n-1))/dphi
		S.d(z.n, 1) =    (- D.d(z.n-1)*tau/(h.d(z.n-1)**2) + 0.5 * u.d(z.n-1)*tau/h.d(z.n-1)) * sI**2 + 0.5*(abs(u_phi)+u_phi)*tau/h.d(z.n-1)
		u_phi = -1/R * D.d(i) * sI * cI * (nold(q+1).d(z.n)-nold(q-1).d(z.n))/(nold(q+1).d(i)+nold(q-1).d(i))/dphi
		S.d(z.n, 2) = +1 + (D.d(z.n-1)*tau/(h.d(z.n-1)**2) + 0.5 * u.d( z.n )*tau/h.d(z.n-1)) * sI**2 + 0.5*(abs(u_phi)-u_phi)*tau/h.d(z.n-1)

	else if (nonlinear_scheme_type .eq. 6) then
		!u_{phi(N+1/2)}
		u_phi = -1/R * D.d(i) * sI * cI * (nold(q+1).interp(z, m.d(z.n-1))-nold(q-1).interp(z, m.d(z.n-1)))/(nold(q+1).interp(z, m.d(z.n-1))+nold(q-1).interp(z, m.d(z.n-1)))/dphi
		!u_phi(N)
		u_phi_N = -1/R * D.d(i) * sI * cI * (nold(q+1).d(z.n-1)-nold(q-1).d(z.n-1))/(nold(q+1).d(z.n-1)+nold(q-1).d(z.n-1))/dphi
		!u_phi(N+1)
		u_phi_Np1 = -1/R * D.d(i) * sI * cI * (nold(q+1).d(z.n)-nold(q-1).d(z.n))/(nold(q+1).d(z.n)+nold(q-1).d(z.n))/dphi
		S.d(z.n, 1) =    (- D.d(z.n-1)*tau/(h.d(z.n-1)**2) + 0.5 * u.d(z.n-1)*tau/h.d(z.n-1)) * sI**2 + 0.5*(u_phi_N-abs(u_phi))*tau/h.d(z.n-1)
		S.d(z.n, 2) = +1 + (D.d(z.n-1)*tau/(h.d(z.n-1)**2) + 0.5 * u.d( z.n )*tau/h.d(z.n-1)) * sI**2 + 0.5*(u_phi_Np1+abs(u_phi))*tau/h.d(z.n-1)

	else if (nonlinear_scheme_type .eq. 7) then
		!u_{phi(N+1/2)}
		if(sI .ge. 0) then
			u_phi_ph = -2/R * D.d(i) * sI * cI * (nold(q+1).interp(z, m.d(z.n-1))-nold(q).interp(z, m.d(z.n-1)))/(nold(q+1).interp(z, m.d(z.n-1))+nold(q-1).interp(z, m.d(z.n-1)))/dphi
		else
			u_phi_ph = -2/R * D.d(i) * sI * cI * (nold(q).interp(z, m.d(z.n-1))-nold(q-1).interp(z, m.d(z.n-1)))/(nold(q+1).interp(z, m.d(z.n-1))+nold(q-1).interp(z, m.d(z.n-1)))/dphi
		end if
		S.d(z.n, 1) =    (- D.d(z.n-1)*tau/(h.d(z.n-1)**2) + 0.5 * u.d(z.n-1)*tau/h.d(z.n-1)) * sI**2 + 0.5*(u_phi_ph+abs(u_phi_ph))*tau/h.d(z.n-1)
		S.d(z.n, 2) = +1 + (D.d(z.n-1)*tau/(h.d(z.n-1)**2) + 0.5 * u.d( z.n )*tau/h.d(z.n-1)) * sI**2 + 0.5*(u_phi_ph-abs(u_phi_ph))*tau/h.d(z.n-1)

	end if

	do i = 2, z.n - 1
	if (nold(q+1).d(i) .eq. 0 .or. nold(q-1).d(i) .eq. 0) then

	! symmetric scheme without u_{phi} = 1/a D sinIcosI ln(n_{q+1}/n_{q-1}) / 2dphi

		S.d(i, 1) = (-D.d(i-1)*tau/(hmid.d(i) * h.d(i-1)) + u.d(i-1)*tau/(h.d(i) + h.d(i-1))) * sI**2
		S.d(i, 2) = 1 + k.d(i)*tau + (D.d(i-1)/h.d(i-1) + D.d(i)/h.d(i)) * tau / hmid.d(i) * sI**2
		S.d(i, 3) = (-D.d( i )*tau/(hmid.d(i) * h.d( i )) - u.d(i+1)*tau/(h.d(i) + h.d(i-1))) * sI**2

	else if (nonlinear_scheme_type .eq. 1) then

	! symmetric scheme with u_{phi} ~ 1/n d(ln n)/dphi

		u_phi = -1/R * D.d(i) * sI * cI * log(nold(q+1).d(i-1)/nold(q-1).d(i-1))/(2*dphi)
		S.d(i, 1) = (-D.d(i-1)*tau/(hmid.d(i) * h.d(i-1)) + u.d(i-1)*tau/(h.d(i) + h.d(i-1))) * sI**2 + u_phi*tau/(h.d(i)+h.d(i-1))

		S.d(i, 2) = 1 + k.d(i)*tau + (D.d(i-1)/h.d(i-1) + D.d(i)/h.d(i)) * tau / hmid.d(i) * sI**2
 
		u_phi = -1/R * D.d(i) * sI * cI * log(nold(q+1).d(i+1)/nold(q-1).d(i+1))/(2*dphi)
		S.d(i, 3) = (-D.d( i )*tau/(hmid.d(i) * h.d( i )) - u.d(i+1)*tau/(h.d(i) + h.d(i-1))) * sI**2 - u_phi*tau/(h.d(i)+h.d(i-1)) 

	else if (nonlinear_scheme_type .eq. 2) then

	! symmetric scheme with u_{phi} ~ 1/n d(ln n)/dphi

		u_phi = -1/R * D.d(i) * sI * cI * (nold(q+1).d(i-1)-nold(q-1).d(i-1))/(nold(q+1).d(i-1)+nold(q-1).d(i-1))/dphi
		S.d(i, 1) = (-D.d(i-1)*tau/(hmid.d(i) * h.d(i-1)) + u.d(i-1)*tau/(h.d(i) + h.d(i-1))) * sI**2 + u_phi*tau/(h.d(i)+h.d(i-1))

		S.d(i, 2) = 1 + k.d(i)*tau + (D.d(i-1)/h.d(i-1) + D.d(i)/h.d(i)) * tau / hmid.d(i) * sI**2 

		u_phi = -1/R * D.d(i) * sI * cI * (nold(q+1).d(i+1)-nold(q-1).d(i+1))/(nold(q+1).d(i+1)+nold(q-1).d(i+1))/dphi
		S.d(i, 3) = (-D.d( i )*tau/(hmid.d(i) * h.d( i )) - u.d(i+1)*tau/(h.d(i) + h.d(i-1))) * sI**2 - u_phi*tau/(h.d(i)+h.d(i-1))

	else if (nonlinear_scheme_type .eq. 3 .or. nonlinear_scheme_type .eq. 4) then
	!directed difference scheme

		!u_phi = -1/R * D.d(i) * sI * cI * log(nold(q+1).d(i)/nold(q-1).d(i))/(2*dphi)
		u_phi = -1/R * D.d(i) * sI * cI * (nold(q+1).d(i)-nold(q-1).d(i))/(nold(q+1).d(i)+nold(q-1).d(i))/dphi
		S.d(i, 1) = ((-D.d(i-1)/(hmid.d(i) * h.d(i-1)) + u.d(i-1)/(h.d(i) + h.d(i-1))) * sI**2 + (abs(u_phi)-u_phi)/(2*h.d(i-1)))*tau
		
		S.d(i, 2) = 1 + (k.d(i) + (D.d(i-1)/h.d(i-1) + D.d(i)/h.d(i)) / hmid.d(i) * sI**2 + &
				(abs(u_phi)+u_phi)/(2*h.d(i)) - (abs(u_phi)-u_phi)/(2*h.d(i-1)))*tau
		
		S.d(i, 3) = ((-D.d( i )/(hmid.d(i) * h.d( i )) - u.d(i+1)/(h.d(i) + h.d(i-1))) * sI**2 - &
				(abs(u_phi)+u_phi)/(2*h.d(i)))*tau

	else if (nonlinear_scheme_type .eq. 5) then
	!directed difference scheme with u_phi_{i+-1/2}

		!u_phi_{i-1/2}		
		u_phi_mh = -1/R * D.d(i) * sI * cI * (nold(q+1).interp(z, m.d(i-1))-nold(q-1).interp(z, m.d(i-1)))/(nold(q+1).interp(z, m.d(i-1))+nold(q-1).interp(z, m.d(i-1)))/dphi
		!u_phi_{i+1/2}
		u_phi_ph = -1/R * D.d(i) * sI * cI * (nold(q+1).interp(z, m.d(i))-nold(q-1).interp(z, m.d(i)))/(nold(q+1).interp(z, m.d(i))+nold(q-1).interp(z, m.d(i)))/dphi

		S.d(i, 1) = ((-D.d(i-1)/(hmid.d(i) * h.d(i-1)) + u.d(i-1)/(h.d(i) + h.d(i-1))) * sI**2 + (abs(u_phi_mh)-u_phi_mh)/(2*h.d(i-1)))*tau
		S.d(i, 2) = 1 + (k.d(i) + (D.d(i-1)/h.d(i-1) + D.d(i)/h.d(i)) / hmid.d(i) * sI**2 + &
				(abs(u_phi_ph)+u_phi_ph)/(2*h.d(i)) - (abs(u_phi_mh)-u_phi_mh)/(2*h.d(i-1)))*tau	
		S.d(i, 3) = ((-D.d( i )/(hmid.d(i) * h.d( i )) - u.d(i+1)/(h.d(i) + h.d(i-1))) * sI**2 - &
				(abs(u_phi_ph)+u_phi_ph)/(2*h.d(i)))*tau

	else if (nonlinear_scheme_type .eq. 6) then
	!KIR scheme
		!u_phi_{i-1/2}		
		u_phi_mh = -1/R * D.d(i) * sI * cI * (nold(q+1).interp(z, m.d(i-1))-nold(q-1).interp(z, m.d(i-1)))/(nold(q+1).interp(z, m.d(i-1))+nold(q-1).interp(z, m.d(i-1)))/dphi
		!u_phi_{i+1/2}
		u_phi_ph = -1/R * D.d(i) * sI * cI * (nold(q+1).interp(z, m.d(i))-nold(q-1).interp(z, m.d(i)))/(nold(q+1).interp(z, m.d(i))+nold(q-1).interp(z, m.d(i)))/dphi
		!u_phi(i-1)
		u_phi_m1 = -1/R * D.d(i) * sI * cI * (nold(q+1).d(i-1)-nold(q-1).d(i-1))/(nold(q+1).d(i-1)+nold(q-1).d(i-1))/dphi
		!u_phi(i+1)
		u_phi_p1 = -1/R * D.d(i) * sI * cI * (nold(q+1).d(i+1)-nold(q-1).d(i+1))/(nold(q+1).d(i+1)+nold(q-1).d(i+1))/dphi

		S.d(i, 1) = (-D.d(i-1)*tau/(hmid.d(i) * h.d(i-1)) + u.d(i-1)*tau/(h.d(i) + h.d(i-1))) * sI**2 + u_phi_m1*tau/(h.d(i) + h.d(i-1)) - abs(u_phi_mh)*tau/(h.d(i) + h.d(i-1))
		S.d(i, 2) = 1 + k.d(i)*tau + (D.d(i-1)/h.d(i-1) + D.d(i)/h.d(i)) * tau / hmid.d(i) * sI**2 + 0.5*(abs(u_phi_mh) + abs(u_phi_ph)) * tau / h0
		S.d(i, 3) = (-D.d( i )*tau/(hmid.d(i) * h.d( i )) - u.d(i+1)*tau/(h.d(i) + h.d(i-1))) * sI**2 - u_phi_p1*tau/(h.d(i) + h.d(i-1)) - abs(u_phi_ph)*tau/(h.d(i) + h.d(i-1))

	else if (nonlinear_scheme_type .eq. 7) then
	!New flux scheme
	if(sI .ge. 0) then
		!u_phi_{i-1/2}		
		u_phi_mh = -2/R * D.d(i) * sI * cI * (nold(q+1).interp(z, m.d(i-1))-nold(q).interp(z, m.d(i-1)))/(nold(q+1).interp(z, m.d(i-1))+nold(q-1).interp(z, m.d(i-1)))/dphi
		!u_phi_{i+1/2}
		u_phi_ph = -2/R * D.d(i) * sI * cI * (nold(q+1).interp(z, m.d(i))-nold(q).interp(z, m.d(i)))/(nold(q+1).interp(z, m.d(i))+nold(q-1).interp(z, m.d(i)))/dphi
		!u_phi(i-1)
!		u_phi_m1 = -2/R * D.d(i) * sI * cI * (nold(q+1).d(i-1)-nold(q).d(i-1))/(nold(q+1).d(i-1)+nold(q-1).d(i-1))/dphi
		!u_phi(i+1)
!		u_phi_p1 = -2/R * D.d(i) * sI * cI * (nold(q+1).d(i+1)-nold(q).d(i+1))/(nold(q+1).d(i+1)+nold(q-1).d(i+1))/dphi
	else
		!u_phi_{i-1/2}		
		u_phi_mh = -2/R * D.d(i) * sI * cI * (nold(q).interp(z, m.d(i-1))-nold(q-1).interp(z, m.d(i-1)))/(nold(q+1).interp(z, m.d(i-1))+nold(q-1).interp(z, m.d(i-1)))/dphi
		!u_phi_{i+1/2}
		u_phi_ph = -2/R * D.d(i) * sI * cI * (nold(q).interp(z, m.d(i))-nold(q-1).interp(z, m.d(i)))/(nold(q+1).interp(z, m.d(i))+nold(q-1).interp(z, m.d(i)))/dphi
		!u_phi(i-1)
!		u_phi_m1 = -2/R * D.d(i) * sI * cI * (nold(q).d(i-1)-nold(q-1).d(i-1))/(nold(q+1).d(i-1)+nold(q-1).d(i-1))/dphi
		!u_phi(i+1)
!		u_phi_p1 = -2/R * D.d(i) * sI * cI * (nold(q).d(i+1)-nold(q-1).d(i+1))/(nold(q+1).d(i+1)+nold(q-1).d(i+1))/dphi
	end if

		S.d(i, 1) = (-D.d(i-1)*tau/(hmid.d(i) * h.d(i-1)) + u.d(i-1)*tau/(h.d(i) + h.d(i-1))) * sI**2 + (u_phi_mh+abs(u_phi_mh))*tau/(2*h0)
		S.d(i, 2) = 1 + k.d(i)*tau + (D.d(i-1)/h.d(i-1) + D.d(i)/h.d(i)) * tau / hmid.d(i) * sI**2 +(u_phi_mh-abs(u_phi_mh))*tau/(2*h0)-(u_phi_ph+abs(u_phi_ph))*tau/(2*h0)
		S.d(i, 3) = (-D.d( i )*tau/(hmid.d(i) * h.d( i )) - u.d(i+1)*tau/(h.d(i) + h.d(i-1))) * sI**2 - (u_phi_ph-abs(u_phi_ph))*tau/(2*h0)



	end if	
	end do

!	print *
!	call S.print()
!	print *

	nnew(q) = tridiagonal_matrix_algorithm(S, b)
!	call nnew(q).print_stdout()

!block to output the diurnal evolution
!	if (q .eq. 44 .and. j .ge. 288) then
!		call nnew(q).print(10)
!		do i = 1, z.n
!			write(11,*) (j-288)*5, 100+400/(z.n-1)*(i-1), nnew(q).d(i)
!		end do
!		write(11, *)
!	end if

!block to output the stationary solution
!	if (j .eq. 288/2) then
!		do i = 1, z.n
!			write(11,*) (q-1)*(1E+0)*180/(Nphi-1), 100+400/(z.n-1)*(i-1), nnew(q).d(i)
!		end do
!		write (11, *)
!	end if
!
!	if (j .eq. 288*1) then
!		do i = 1, z.n
!			write(12,*) (q-1)*(1E+0)*180/(Nphi-1), 100+400/(z.n-1)*(i-1), nnew(q).d(i)
!		end do
!		write (12, *)
!	end if

!	if (j .eq. 288*2) then
!		do i = 1, z.n
!			write(13,*) (q-1)*(1E+0)*180/(Nphi-1), 100+400/(z.n-1)*(i-1), nnew(q).d(i)
!		end do
!		write (13, *)
!	end if
!
!	if (j .eq. 288*5) then
!		do i = 1, z.n
!			write(20,*) (q-1)*(1E+0)*180/(Nphi-1), 100+400/(z.n-1)*(i-1), nnew(q).d(i)
!		end do
!		write (20, *)
!	end if
!

	if (j*tau .eq. 86400/2) then
		call nnew(q).print(105)
	end if

	if (j*tau .eq. 86400) then
		call nnew(q).print(11)
	end if

	if (j*tau .eq. 86400*2) then
		call nnew(q).print(12)
	end if

	if (j*tau .eq. 86400*5) then
		call nnew(q).print(15)
	end if

	if (j*tau .eq. 86400*10) then
		call nnew(q).print(110)
	end if

	if (j*tau .eq. 86400*15) then
		call nnew(q).print(115)
	end if

	if (j*tau .eq. 86400*20) then
		call nnew(q).print(20)
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
