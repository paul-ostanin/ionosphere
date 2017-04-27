program mat

use tridi_matrix
use vector

implicit none

type (tridiagonal_matrix) A, S
type (vect) b, z, h, hmid, nO, nO2, nN2, k, p, pcurr, m, njold, njnew, delta, D, u, Tn, Ti, Te, Tr, Tp, gradTp, nday, tau0
integer i, j, out_file, scheme_type, Te0, Tn0, Ti0, day
real (8) tau, h0, Fub, delta_norm, eps, tgdelta, sindelta, cosdelta, phi, coschi, pi, omega, sigma_O2, sigma_N2, sigma_O

!opening the file res.txt for writing the output
open(unit=10, name='res.txt')

pi = 3.141592653589793238462643

!scheme type
!0 - old monotone scheme
!1 - symmetric scheme with old boundary condition approximation
!2 - symmetric scheme with new boundary condition approximation
scheme_type = 0

!Time step (in seconds) 5 min
tau = 300

!Vector of altitudes. Step h_i = z(i) - z(i - 1). Counting from 100 km to 500 km. z.d(i) is in metres.
call z.init(401)

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

Ti0 = 1.0 * 950 
Tn0 = 1.0 * 800
Te0 = 1.0 * 2200

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
	nO.d(i)  = 1.0 * 2.8E+10 * exp(-9.8 * 16E-3 / (8.31 * Tn.d(i)) * (z.d(i) - 140000))
	nO2.d(i) = 1.0 * 5.6E+9  * exp(-9.8 * 32E-3 / (8.31 * Tn.d(i)) * (z.d(i) - 140000))
	nN2.d(i) = 1.0 * 5.2E+10 * exp(-9.8 * 28E-3 / (8.31 * Tn.d(i)) * (z.d(i) - 140000))
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

!P = P_0 exp(tau_0(z)*(1-sec chi)
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
	p.d(i) = 1.0 * 4E-7 * nO.d(i)
	k.d(i) = 1.0 * (1.2E-12 * nN2.d(i) + 2.1E-11 * nO2.d(i))
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

!System matrix. n(first) and n(last) will be known from boundary conditions.
call S.init(z.n)
!lower boundary condition: n_1 = P_1/k_1
S.d(1, 2) = 1

if (scheme_type == 1 .OR. scheme_type == 0) then
!old upper boundary condition: dn/dz+u_N/D_N*n_N=Fub -> 0 <=> n_N-n_{N-1}+u_N/D_N*h*n_N = Fub*h <=> -n_{N-1}+(1+u_N/D_N*h)*n_N=Fub*h/D_N
S.d(z.n, 1) = -D.d(z.n - 1)/h0**2
S.d(z.n, 2) = D.d(z.n - 1)/h0**2 +  u.d(z.n)/h0

else if (scheme_type == 2) then
!new upper boundary condition:
S.d(z.n, 1) =    - D.d(z.n-1)*tau/(h.d(z.n-1)**2) + 0.5 * u.d(z.n-1)*tau/h.d(z.n-1)
S.d(z.n, 2) = +1 + D.d(z.n-1)*tau/(h.d(z.n-1)**2) + 0.5 * u.d( z.n )*tau/h.d(z.n-1)
end if

!middle matrix rows are similar
do i = 2, z.n - 1

if (scheme_type == 0) then
! old scheme
	S.d(i, 1) = -D.d(i-1)*tau/h0**2
	S.d(i, 2) = (1 + tau * k.d(i)) + u.d(i)*tau/h0 + (D.d(i-1) + D.d(i))*tau/h0**2
	S.d(i, 3) = -D.d( i )*tau/h0**2 - u.d(i+1)*tau/h0

else if (scheme_type == 1 .OR. scheme_type == 2) then
! symmetric scheme
	S.d(i, 1) = -D.d(i-1)*tau/(hmid.d(i) * h.d(i-1)) + u.d(i-1)*tau/(h.d(i) + h.d(i-1))
	S.d(i, 2) = 1 + k.d(i)*tau + (D.d(i-1)/h.d(i-1) + D.d(i)/h.d(i)) * tau / hmid.d(i)
	S.d(i, 3) = -D.d( i )*tau/(hmid.d(i) * h.d( i )) - u.d(i+1)*tau/(h.d(i) + h.d(i-1))
end if

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

j=0
Fub = 0
delta_norm = 1

	do while (delta_norm > 1E-5)
	
		j = j + 1
		!Setting RHS in the middle
	        do i = 2, z.n - 1
	                b.d(i) = njold.d(i) + tau * p.d(i)
	        end do

		!Setting boundary conditions for the RHS
		if(scheme_type == 2) then 
			b.d(z.n) = +tau/h.d(z.n-1) * Fub + njold.d(z.n)
		else if(scheme_type == 0 .OR. scheme_type == 1) then 
			b.d(z.n) = Fub/h0
		end if
		b.d(1) = P.d(1)/ k.d(1)
!		b.d(1) = 0

		!Solving the system
	        njnew = tridiagonal_matrix_algorithm(S, b)

		delta = njnew - njold
		delta_norm = delta.norm()
		njold = njnew

!		print *
!		print *, "time is", tau * j
!		print *
!		call njold.print_old()

	end do

	nday = njold


	call nday.print_column(10)


!A block to calculate the solution evolution after setting P = 0
!do i = 1, z.n
!	p.d(i) = 0
!end do

!do j = 1, 86400/tau
!        do i = 2, z.n - 1
!                b.d(i) = njold.d(i)
!        end do

	!Setting boundary conditions for the RHS
!	if(scheme_type == 2) then 
!		b.d(z.n) = tau/h.d(z.n-1) * Fub + njold.d(z.n)
!	else if(scheme_type == 1) then 
!		b.d(z.n) = Fub/h0
!	end if
!	b.d(1) = 0
!
!	!Solving the system
!        njnew = tridiagonal_matrix_algorithm(S, b)
!
!	njold = njnew
!
!	write(10, *)
!	call njnew.print(10)
!end do
!	close(unit=10)

!A block to calculate day-night evolution


phi = -65 * pi / 180
omega = 2*pi/24/60/60

do j = 1, 86400/tau

	day = (j*tau)/86400 + 1/2 !starting from the afternoon of the 1-st day, day starts from 1/2
	tgdelta = tan(pi/180*23.5) * sin(2*pi/365 * (day - 80))
	sindelta = tgdelta/sqrt(1+tgdelta**2)
	cosdelta = sqrt(1-sindelta**2) !cos of the zenith angle is > 0

        do i = 2, z.n - 1
		coschi = sin(phi)*sindelta - cos(phi)*cosdelta*cos(omega * (tau * j + 86400/2)) !starting from the middle of the day
		if(coschi > 1E-1) then 
		        b.d(i) = njold.d(i) + tau * p.d(i) * exp(tau0.d(i) * (1-1/coschi))
        	else 
			b.d(i) = njold.d(i)
        	end if
        end do

!	print *, j
!	print *, coschi

	!Setting boundary conditions for the RHS
	if(scheme_type == 2) then 
		b.d(z.n) = tau/h.d(z.n-1) * Fub + njold.d(z.n)
	else if(scheme_type == 1) then 
		b.d(z.n) = Fub/h0
	end if
	if(coschi > 1E-2) then 
	        b.d(1) = p.d(1)/k.d(1) * exp(tau0.d(1) * (1-1/coschi))
       	else 
		b.d(1) = 0
       	end if

	!Solving the system
        njnew = tridiagonal_matrix_algorithm(S, b)

	do i = 1, z.n
		if(coschi > 1E-2) then 
			pcurr.d(i) = p.d(i) * exp(tau0.d(i) * (1-1/cos(65*pi/180)))
	       	else 
			pcurr.d(i) = 0
	       	end if
	end do

	njold = njnew
!	write(10, '(e14.7)') p.d
!	call njnew.print(10)

end do
!	call pcurr.print_column(10)


	close(unit=10)




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
