program mat

use tridi_matrix
use vector

implicit none

type (tridiagonal_matrix) S_y, S_z
type (vect) rhs_z, rhs_y, z, h, hmid, nO, nO2, nN2, k, p, pcurr, m, n_old, n_new, delta, D, D_node, u, Tn, Ti, Te, Tr, Tp, gradTp, n_day, tau0, n_new_z(1441), n_old_z(1441), n_new_y(401), n_old_y(401)
integer i, j, t, Te0, Tn0, Ti0, day, nonlinear_scheme_type, diurnal_on, Nphi, pk_switch, mixed_z_switch, mixed_y_switch, transf_yz_switch, transf_y_switch
real (8) tau, h0, F_z, delta_norm, eps, tgdelta, sindelta, cosdelta, dphi, phi, coschi, pi, omega, sigma_O2, sigma_N2, sigma_O, sI, cI, R, u_phi, u_phi_1, u_phi_2, u_phi_3, u_phi_4, u_z, u_z_mh, u_z_ph, u_z_m1, u_z_p1, x, A, B, u_phi_ph, u_phi_mh, Ndays, Niter

!opening files for writing the output
open(unit=10, name='res.txt')
open(unit=11, name='res_gnp.txt')

pi = 3.141592653589793238462643

!nonlinear_scheme_type variable switches the u_phi-approximation. 
nonlinear_scheme_type = 7

!number of nodes in phi
Nphi = 180
!latitude
dphi = pi / Nphi
!angle velocity of the Earth 
omega = 2*pi/24/60/60
!magnetic inclination sin I
sI = 1
!Earth radius
R = 637100000
!number of calculation days
Ndays = 1
Niter = 800

!Time step (in seconds) 5 min
tau = 300

!switches for physical processes and terms of the equation
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
	p.d(i) = ( 4E-7 * nO.d(i) ) * pk_switch
	k.d(i) = ( 1.2E-12 * nN2.d(i) + 2.1E-11 * nO2.d(i) ) * pk_switch
end do

!print *
!print *, "Printing P (in cm^{-3}*s^{-1}):"
!call p.print()

!print *
!print *, "Printing k (in s^{-1}):"
!call k.print_long_frac()


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

!print *
!print *, "Printing D_{i+1/2} (in cm^2*s^{-1}):"
!call D.print()


!Effective velocity vector
!u.d(i) = u_i = D_i * (dTp/dz + mg/2k)/Tp, mg/2k ~ 5.6*10^{-5} [K/cm]
call u.init(z.n)
do i = 1, z.n
	u.d(i) = ( 3E+17 / (nO.d(i) * sqrt(Tr.d(i))) * (56E-6 + gradTp.d(i)) ) * transf_yz_switch
end do

!print *
!print *, "Printing effective velocity vector (in cm/s):"
!call u.print()

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
call pcurr.init(z.n)

F_z = 0

do j = 1, Nphi
	call n_new_z(j).init(z.n)
	call n_old_z(j).init(z.n)
	call n_old_z(j).gen()
end do

do i = 1, z.n
	call n_new_y(i).init(Nphi)
	call n_old_y(i).init(Nphi)
end do


!do i = 1, z.n
!	write(11,*) -90, 100+400/(z.n-1)*(i-1), n_day.d(i)
!end do
!	write (11, *)

do t = 0, Ndays*86400/tau !Niter + 1!
print *, t
do j = 1, Nphi
! angles phi from -90 to 90; conditions in -90 and 90 are set


        do i = 2, z.n - 1
		rhs_z.d(i) = n_old_z(j).d(i) + tau * p.d(i)
	end do

	rhs_z.d(z.n) = +tau/h.d(z.n-1) * F_z + n_old_z(j).d(z.n)
	rhs_z.d(1) = p.d(1)/k.d(1)*pk_switch + (pk_switch - 1) * (-1)

	!sinus and cosinus of magnetic inclination angle I
	sI = sin(atan(2*tan(-pi/2+(j-0.5)*dphi)))
	cI = cos(atan(2*tan(-pi/2+(j-0.5)*dphi)))

	!lower boundary condition: n_1 = P_1/k_1
	S_z.d(1, 2) = 1
	!upper boundary condition:
		S_z.d(z.n, 1) =    (- D.d(z.n-1)*tau/(h.d(z.n-1)**2) + 0.5 * u.d(z.n-1)*tau/h.d(z.n-1)) * sI**2
		S_z.d(z.n, 2) = +1 + (D.d(z.n-1)*tau/(h.d(z.n-1)**2) + 0.5 * u.d( z.n )*tau/h.d(z.n-1)) * sI**2
	

	do i = 2, z.n - 1
	! symmetric scheme without u_{phi} = 1/a D sinIcosI ln(n_{j+1}/n_{j-1}) / 2dphi

		S_z.d(i, 1) = (-D.d(i-1)*tau/(hmid.d(i) * h.d(i-1)) + u.d(i-1)*tau/(h.d(i) + h.d(i-1))) * sI**2

		S_z.d(i, 2) = 1 + k.d(i)*tau + (D.d(i-1)/h.d(i-1) + D.d(i)/h.d(i)) * tau / hmid.d(i) * sI**2

		S_z.d(i, 3) = (-D.d( i )*tau/(hmid.d(i) * h.d( i )) - u.d(i+1)*tau/(h.d(i) + h.d(i-1))) * sI**2


	end do

	n_new_z(j) = tridiagonal_matrix_algorithm(S_z, rhs_z)

end do

!sending 1st step result to the 2nd step as the initial condition
do i = 1, z.n
	do j = 1, Nphi
		n_old_y(i).d(j) = n_new_z(j).d(i)
	end do

end do

	!boundary conditions at z = 100 and z = 500
	do j = 1, Nphi
		n_new_y(1).d(j) = n_new_z(j).d(1)
		n_new_y(z.n).d(j) = n_new_z(j).d(z.n)
	end do

!simulating block for the 2nd step
!do i = 2, z.n-1
!do j = 2, Nphi
!	n_new_y(i).d(j) = n_old_y(i).d(j)
!end do
!end do


do i = 2, z.n-1

	do j = 1, Nphi
		rhs_y.d(j) = n_old_y(i).d(j) + 0*tau*p.d(i)/2
	end do


	j = 1
	phi = (j-0.5)*dphi-pi/2


		S_y.d(j, 1) =  0
		S_y.d(j, 2) = 1 + (A(phi+dphi/2)/(dphi**2)-transf_y_switch*(+u.d(i)/2)*B(phi)/(2*dphi))*tau/(R)/cos(phi)
		S_y.d(j, 3) =  (-D_node.d(i)/R * A(phi+dphi/2)/(dphi**2) - &
				transf_y_switch*(-u.d(i)/2)*B(phi+dphi)/(2*dphi))*tau/(R)/cos(phi)



	do j = 2, Nphi-1	
		phi = (j-0.5)*dphi-pi/2

		S_y.d(j, 1) =  (-D_node.d(i)/R * A(phi-dphi/2)/(dphi**2) + &
				transf_y_switch*(-u.d(i)/2)*B(phi-dphi)/(2*dphi))*tau/(R)/cos(phi)
		S_y.d(j, 2) = 1 + (D_node.d(i)/R * (A(phi-dphi/2) + A(phi+dphi/2))/(dphi**2)) * tau/(R)/cos(phi)
		S_y.d(j, 3) =  (-D_node.d(i)/R * A(phi+dphi/2)/(dphi**2) - &
				transf_y_switch*(-u.d(i)/2)*B(phi+dphi)/(2*dphi))*tau/(R)/cos(phi)


	end do

	j = Nphi
	phi = (j-0.5)*dphi-pi/2

		S_y.d(j, 1) =  (-D_node.d(i)/R * A(phi-dphi/2)/(dphi**2) + &
				 transf_y_switch*(-u.d(i)/2)*B(phi-dphi)/(2*dphi))*tau/(R)/cos(phi)
		S_y.d(j, 2) = 1 + (D_node.d(i)/R * A(phi-dphi/2)/(dphi**2) + &
				 transf_y_switch*(+u.d(i)/2)*B(   phi  )/(2*dphi))*tau/(R)/cos(phi)
		S_y.d(j, 3) = 0


	n_new_y(i) = tridiagonal_matrix_algorithm(S_y, rhs_y)

end do


!sending the result back to the 1-st steps
do i = 1, z.n
	do j = 1, Nphi
		n_old_z(j).d(i) = n_new_y(i).d(j)
	end do
end do

!block to output the stationary solution
	if (t*tau .eq. 86400*Ndays) then 
		do j = 1, Nphi
		do i = 1, z.n
			write(11,*) (j - 5E-1)*180/(Nphi)-90, 100+400/(z.n-1)*(i-1), n_new_z(j).d(i)

		end do
		write (11, *)
		end do

		do j = 1, Nphi
			call n_new_z(j).print(10)
		end do
	end if

end do



 close(unit=10)
 close(unit=11)
do j = 1, Nphi
 call n_new_z(j).destroy()
 call n_old_z(j).destroy()
end do
call rhs_z.destroy()
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
call n_old.destroy()
call delta.destroy()
call n_new.destroy()
call S_z.destroy()
end

real*8 function A(phi)
implicit none
real*8 phi, pi
real*8 smallphi 

pi = 3.141592653589793238462643
!smallphi = -1

!if(abs(phi-pi/2) .le. smallphi) then
!	A = -0.25 * (phi - pi/2)**3
!else if (abs(phi+pi/2) .le. smallphi) then
!	A = 0.25 * (phi + pi/2)**3
!else
	A = cos(phi) / (1+4*(tan(phi)**2))
!end if

return
end

real*8 function B(phi)
implicit none
real*8 phi, pi
real*8 smallphi 

pi = 3.141592653589793238462643
!smallphi = -1

!if(abs(phi-pi/2) .le. smallphi) then
!	B = (phi - pi/2)**2
!else if (abs(phi+pi/2) .le. smallphi) then
!	B = -(phi + pi/2)**2
!else
	B = 4*sin(phi) / (1+4*(tan(phi)**2))
!end if

return
end

