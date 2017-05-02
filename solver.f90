program mat

use tridi_matrix
use vector

implicit none

type (tridiagonal_matrix) A, S
type (vect) b, z, h, hmid, nO, nO2, nN2, k, p, pcurr, m, njold, njnew, delta, D, u, Tn, Ti, Te, Tr, Tp, gradTp, nday, tau0, nnew(91), nold(91)
integer i, j, q, Te0, Tn0, Ti0, day, nonlinear_scheme_type, diurnal_on
real (8) tau, h0, Fub, delta_norm, eps, tgdelta, sindelta, cosdelta, dphi, phi, coschi, pi, omega, sigma_O2, sigma_N2, sigma_O, sI, cI, R, u_phi, u_phi_N, u_phi_Nm1

!opening the file res.txt for writing the output
open(unit=10, name='res.txt')
open(unit=11, name='res_gnp.txt')

pi = 3.141592653589793238462643

!nonlinear_scheme_type variable switches the u_phi-approximation. 
!If n_s_t = 1, u_phi = 1/n d(ln n)/dphi
!If n_s_t = 2, u_phi = (n(phi+dphi)-n(phi-dphi))/(n(phi+dphi)+n(phi-dphi)) * 1/dphi
nonlinear_scheme_type = 2


!latitude
dphi = 2 * pi / 180
!angle velocity of the earth 
omega = 2*pi/24/60/60
!magnetic inclination sin I
sI = 1
!Earth radius
R = 637100000



!Time step (in seconds) 5 min
tau = 300

!Vector of altitudes. Step h_i = z(i) - z(i - 1). Counting from 100 km to 500 km. z.d(i) is in metres.
call z.init(81)

!Space step (in cm) 5 km
h0 = 400E+5 / (z.n - 1)

do i = 1, z.n
	z.d(i) = 100E+3 + h0 * (i-1)/100
end do

!Vector of middles of [z(i), z(i+1)].  m.d(i) = z_{i+1/2} (in metres)
call m.init(z.n - 1)
do i = 1, z.n - 1
	m.d(i) = 0.5 * (z.d(i) + z.d(i + 1))
end do

call h.init(z.n - 1)
do i = 1, z.n - 1
	h.d(i) = 100 * (z.d(i + 1) - z.d(i))
end do

call hmid.init(z.n - 1)
do i = 1, z.n - 1
	hmid.d(i) = 100 * (m.d(i) - m.d(i - 1))
end do

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

!print *, "Temperatures (in K):"
!print *

!print *, "Printing T_i:"
!call Ti.print()
!print *, "Printing T_n:"
!call Tn.print()
!print *, "Printing T_e:"
!call Te.print()
!print *, "Printing T_r:"
!call Tr.print()
!print *, "Printing T_p:"
!call Tp.print()


!dT/dz = s(T_\infty - T_0)exp(-s(z-z0)) = (T_\infty - T)*s - for Te, Tn, Ti. This derivative is in [K/cm].
call gradTp.init(z.n)
do i = 1, z.n
	gradTp.d(i) = 0.5 * ((Te0 - Te.d(i)) * 9.8 / (287 * Te0) + (Ti0 - Ti.d(i)) * 9.8 / (287 * Ti0) ) * 1E-2
end do

!print *
!print *, "Printing dTp/dz (in K/cm) :"
!call gradTp.print_long_frac()

call nO.init(z.n)
call nO2.init(z.n)
call nN2.init(z.n)

do i = 1, z.n
	nO.d(i)  = 2.8E+10 * exp(-9.8 * 16E-3 / (8.31 * Tn.d(i)) * (z.d(i) - 140000))
	nO2.d(i) = 5.6E+9  * exp(-9.8 * 32E-3 / (8.31 * Tn.d(i)) * (z.d(i) - 140000))
	nN2.d(i) = 5.2E+10 * exp(-9.8 * 28E-3 / (8.31 * Tn.d(i)) * (z.d(i) - 140000))
end do

!print *
!print *, "Concentrations (in cm^-3):"
!print *

!print *, "Printing n_O:"
!call nO.print()
!print *, "Printing n_O2:"
!call nO2.print() 
!print *, "Printing n_N2:"
!call nN2.print()

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

!print *, "Printing tau0:"
!call tau0.print_column(10)

call p.init(z.n)
call k.init(z.n)
do i = 1, z.n
	p.d(i) = 4E-7 * nO.d(i)
	k.d(i) = 1.2E-12 * nN2.d(i) + 2.1E-11 * nO2.d(i)
end do

!print *
!print *, "Printing P (in cm^{-3}*s^{-1}):"
!call p.print()

!print *
!print *, "Printing k (in s^{-1}):"
!call k.print_long_frac()


!Diffusion coefficients vector. D.d(i) = D_{i+1/2}
call D.init(z.n - 1)
do i = 1, z.n - 1
	D.d(i) = 3E+17 * Tp.interp(z, m.d(i)) / (nO.interp(z, m.d(i)) * sqrt(Tr.interp(z, m.d(i))))
end do

!print *
!print *, "Printing D_{i+1/2} (in cm^2*s^{-1}):"
!call D.print()


!Effective velocity vector
!u.d(i) = u_i = D_i * (dTp/dz + mg/2k)/Tp, mg/2k ~ 5.6*10^{-5} [K/cm]
call u.init(z.n)
do i = 1, z.n
	u.d(i) = 3E+17 / (nO.d(i) * sqrt(Tr.d(i))) * (56E-6 + gradTp.d(i)) 
end do

!print *
!print *, "Printing effective velocity vector (in cm/s):"
!call u.print()

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
	
		j = j + 1
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

do j = 0, 90
	call nnew(j).init(z.n)
	call nold(j).init(z.n)
	call nold(j).gen()
end do

	nold( 0) = nday
	nold(90) = nday

do j = 0, 288+288 
do q = 1, 89
! angles phi from -90 to 90; conditions in -90 and 90 are set, the step is dphi = 2

diurnal_on = 1 !switcher for the diurnal evolution mode. 0 corresponds to stationary solutions.

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
	sI = sin(atan(2*tan(-pi/2+q*dphi)))
	cI = cos(atan(2*tan(-pi/2+q*dphi)))

	!lower boundary condition: n_1 = P_1/k_1
	S.d(1, 2) = 1
	!new upper boundary condition:
	if (nold(q+1).d(z.n) .eq. 0 .or. nold(q-1).d(z.n) .eq. 0 .or. nold(q+1).d(z.n-1) .eq. 0 .or. nold(q+1).d(z.n-1) .eq. 0) then
		S.d(z.n, 1) =    (- D.d(z.n-1)*tau/(h.d(z.n-1)**2) + 0.5 * u.d(z.n-1)*tau/h.d(z.n-1)) * sI**2
		S.d(z.n, 2) = +1 + (D.d(z.n-1)*tau/(h.d(z.n-1)**2) + 0.5 * u.d( z.n )*tau/h.d(z.n-1)) * sI**2
	else if (nonlinear_scheme_type .eq. 1) then
		u_phi = -1/R * D.d(i) * sI * cI * log(nold(q+1).d(z.n-1)/nold(q-1).d(z.n-1))/(2*dphi)
		S.d(z.n, 1) =    (- D.d(z.n-1)*tau/(h.d(z.n-1)**2) + 0.5 * u.d(z.n-1)*tau/h.d(z.n-1)) * sI**2 + 0.5*u_phi*tau/h.d(z.n-1)
		u_phi = -1/R * D.d(i) * sI * cI * log(nold(q+1).d(z.n)/nold(q-1).d(z.n))/(2*dphi)
		S.d(z.n, 2) = +1 + (D.d(z.n-1)*tau/(h.d(z.n-1)**2) + 0.5 * u.d( z.n )*tau/h.d(z.n-1)) * sI**2 + 0.5*u_phi*tau/h.d(z.n-1)
	else if (nonlinear_scheme_type .eq. 2) then
		u_phi = -1/R * D.d(i) * sI * cI * (nold(q+1).d(z.n-1)-nold(q-1).d(z.n-1))/(nold(q+1).d(z.n-1)+nold(q-1).d(z.n-1))/dphi
		S.d(z.n, 1) =    (- D.d(z.n-1)*tau/(h.d(z.n-1)**2) + 0.5 * u.d(z.n-1)*tau/h.d(z.n-1)) * sI**2 + 0.5*u_phi*tau/h.d(z.n-1)
		u_phi = -1/R * D.d(i) * sI * cI * (nold(q+1).d(z.n)-nold(q-1).d(z.n))/(nold(q+1).d(i)+nold(q-1).d(i))/dphi
		S.d(z.n, 2) = +1 + (D.d(z.n-1)*tau/(h.d(z.n-1)**2) + 0.5 * u.d( z.n )*tau/h.d(z.n-1)) * sI**2 + 0.5*u_phi*tau/h.d(z.n-1)
	end if

	do i = 2, z.n - 1
	if (nold(q+1).d(i) .eq. 0 .or. nold(q-1).d(i) .eq. 0) then

	! symmetric scheme without u_{phi} = 1/a D sinIcosI ln(n_{q+1}/n_{q-1}) / 2dphi
		S.d(i, 1) = (-D.d(i-1)*tau/(hmid.d(i) * h.d(i-1)) + u.d(i-1)*tau/(h.d(i) + h.d(i-1))) * sI**2
		S.d(i, 2) = 1 + k.d(i)*tau + (D.d(i-1)/h.d(i-1) + D.d(i)/h.d(i)) * tau / hmid.d(i) * sI**2
		S.d(i, 3) = (-D.d( i )*tau/(hmid.d(i) * h.d( i )) - u.d(i+1)*tau/(h.d(i) + h.d(i-1))) * sI**2
	else if (nonlinear_scheme_type .eq. 1) then
	! symmetric scheme with u_{phi} ~ 1/n d(ln n)/dphi
		u_phi = -1/R * D.d(i) * sI * cI * log(nold(q+1).d(i)/nold(q-1).d(i))/(2*dphi)

		S.d(i, 1) = (-D.d(i-1)*tau/(hmid.d(i) * h.d(i-1)) + u.d(i-1)*tau/(h.d(i) + h.d(i-1))) * sI**2 + u_phi*tau/(h.d(i)+h.d(i-1))
		S.d(i, 2) = 1 + k.d(i)*tau + (D.d(i-1)/h.d(i-1) + D.d(i)/h.d(i)) * tau / hmid.d(i) * sI**2 
		S.d(i, 3) = (-D.d( i )*tau/(hmid.d(i) * h.d( i )) - u.d(i+1)*tau/(h.d(i) + h.d(i-1))) * sI**2 - u_phi*tau/(h.d(i)+h.d(i-1)) 
	else if (nonlinear_scheme_type .eq. 2) then
	! symmetric scheme with u_{phi} ~ 1/n d(ln n)/dphi
		u_phi = -1/R * D.d(i) * sI * cI * (nold(q+1).d(i)-nold(q-1).d(i))/(nold(q+1).d(i)+nold(q-1).d(i))/dphi

		S.d(i, 1) = (-D.d(i-1)*tau/(hmid.d(i) * h.d(i-1)) + u.d(i-1)*tau/(h.d(i) + h.d(i-1))) * sI**2 + u_phi*tau/(h.d(i)+h.d(i-1))
		S.d(i, 2) = 1 + k.d(i)*tau + (D.d(i-1)/h.d(i-1) + D.d(i)/h.d(i)) * tau / hmid.d(i) * sI**2 
		S.d(i, 3) = (-D.d( i )*tau/(hmid.d(i) * h.d( i )) - u.d(i+1)*tau/(h.d(i) + h.d(i-1))) * sI**2 - u_phi*tau/(h.d(i)+h.d(i-1)) 
	end if	
	end do

!if (q .eq. 10 .and. j .eq. 200) then
!	print *
!	call S.print()
!	print *
!end if

	nnew(q) = tridiagonal_matrix_algorithm(S, b)


!	if (q .eq. 20 ) then
!		call nnew(q).print(10)
!		do i = 1, z.n
!			write(11,*) (j-288)*5, 100+400/(z.n-1)*(i-1), nnew(q).d(i)
!		end do
!		write(11, *)
!	end if

	if (q .eq. 44 .and. j.eq. 288) then
		do i = 1, z.n
			write(11,*) 100+400/(z.n-1)*(i-1), nnew(q).d(i)
		end do
	end if


end do	

	do i = 1, 89
		nold(i) = nnew(i)
	end do

end do


 close(unit=10)
do j = 0, 90
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
