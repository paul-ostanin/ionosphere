program mat

use tridi_matrix
use vector

implicit none

type (tridiagonal_matrix) A, S
type (vect) b, z, nO, nO2, nN2, res, k, p, m, njold, njnew, nday, nnight, delta, D, u, Tn, Ti, Te, Tr, Tp, gradTp, Hp, t_char, dif_mul
integer ios, i, j
real (8) tau, h, z0, Fub, delta_norm, eps

!Time step (in seconds) 5 min
tau = 300

!Space step (in cm) 5 km
h = 50E+5

!Vector of altitudes. Step h_i = z(i) - z(i - 1). Counting from 100 km to 500 km. z.d(i) is in metres.
call z.init(9)
do i = 1, z.n
	z.d(i) = 50E+3 + 50E+3 * i
end do

!Vector of middles of [z(i), z(i+1)].  m.d(i) = z_{i+1/2} (in metres)
call m.init(z.n - 1)
do i = 1, z.n - 1
	m.d(i) = 0.5 * (z.d(i) + z.d(i + 1))
end do

call Ti.init(z.n)
call Tn.init(z.n)
call Te.init(z.n)
call Tr.init(z.n)
call Tp.init(z.n)

do i = 1, z.n
        Ti.d(i) = 1.0 * 950  - (1.0 * 950  - 200) * exp(-9.8 / (287 * 1.0 * 950 ) * (z.d(i) - 100000))
        Tn.d(i) = 1.0 * 800  - (1.0 * 800  - 200) * exp(-9.8 / (287 * 1.0 * 800 ) * (z.d(i) - 100000))
        Te.d(i) = 1.0 * 2200 - (1.0 * 2200 - 200) * exp(-9.8 / (287 * 1.0 * 2200) * (z.d(i) - 100000))
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
	gradTp.d(i) = 0.5 * ((1.0 * 2200 - Te.d(i)) * 9.8 / (287 * 1.0 * 2200) + (1.0 * 950 - Ti.d(i)) * 9.8 / (287 * 1.0 * 950) ) * 1E-2
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

call Hp.init(z.n - 1)
do i = 1, z.n - 1
	Hp.d(i) = Tp.interp(z, m.d(i))/ ( (56E-6 + gradTp.interp(z, m.d(i))) )
end do

print *
print *, "Printing Hp (in сm):"
call Hp.print()

call t_char.init(z.n - 1)
do i = 1, z.n - 1
	t_char.d(i) = Hp.d(i) * Hp.d(i) / (3E+17 * Tp.d(i) / (nO.d(i) * sqrt(Tr.d(i))))
end do

!print *
!print *, "Printing t_character (in s):"
!call t_char.print()

call dif_mul.init(z.n - 1)
do i = 1, z.n - 1
	dif_mul.d(i) = 1+h/Hp.d(i)
end do

print *
print *, "Printing diffusion multiplyer:"
call dif_mul.print()


!Diffusion coefficients vector. D.d(i) = D_{i+1/2}
call D.init(z.n - 1)
do i = 1, z.n - 1
	D.d(i) = 3E+17 * Tp.interp(z, m.d(i)) / (nO.interp(z, m.d(i)) * sqrt(Tr.interp(z, m.d(i))))
!	D.d(i) = D.d(i) * (1+h/Hp.d(i))
end do

print *
print *, "Printing D_{i+1/2} (in cm^2*s^{-1}):"
call D.print()



!Effective velocity vector
!u.d(i) = u_i = D_i * (dTp/dz + mg/2k)/Tp, mg/2k ~ 5.6*10^{-5} [K/cm]
call u.init(z.n - 1)
do i = 1, z.n - 1
	u.d(i) = 3E+17 / (nO.d(i) * sqrt(Tr.d(i))) * (56E-6 + gradTp.d(i)) 
end do

!print *
!print *, "Printing effective velocity vector (in cm/s):"
!call u.print()



!System matrix. n(1) and n(11) will be known from boundary conditions.
call S.init(z.n)
!lower boundary condition: n_1 = P_1/k_1
S.d(1, 2) = 1
!upper boundary condition: dn/dz+u_N/D_N*n_N=Fub -> 0 <=> n_N-n_{N-1}+u_N/D_N*h*n_N = Fub*h <=> -n_{N-1}+(1+u_N/D_N*h)*n_N=Fub*h/D_N
S.d(z.n, 1) = -1
S.d(z.n, 2) = 1 + h * (56E-6 + 0.5 * ((1.0 * 2200 - Te.d(z.n)) * 9.8 / (287 * 1.0 * 2200) + (1.0 * 950 - Ti.d(z.n)) * 9.8 / (287 * 1.0 * 950) ) * 1E-2)/Tp.d(z.n)
!middle matrix rows are similar
do i = 2, z.n - 1

! old scheme

!	S.d(i, 1) = -D.d(i - 1) * tau / (h * h)
!	S.d(i, 2) = (1 + tau * k.d(i)) + u.d(i) * tau / h + (D.d(i - 1) + D.d(i))* tau / (h * h) 
!	S.d(i, 3) = -D.d(i)* tau / (h * h) - u.d(i + 1) * tau / h


! symmetric scheme

	S.d(i, 1) = -D.d(i - 1) * tau / (100 * (m.d(i)-m.d(i-1)) * 100 * (z.d(i)-z.d(i-1))) + tau * u.d(i-1)/(100 * (z.d(i+1)-z.d(i-1)))
	S.d(i, 2) = 1 + k.d(i) * tau + ( D.d(i - 1)/(100 * (z.d(i)-z.d(i-1))) + D.d(i)/(100 * (z.d(i+1)-z.d(i))) ) * tau / (100 * (m.d(i)-m.d(i-1)))
	S.d(i, 3) = -D.d(i) * tau / (100 * (m.d(i)-m.d(i-1)) * 100 * (z.d(i+1)-z.d(i))) - tau * u.d(i+1)/(100 * (z.d(i+1)-z.d(i-1)))

end do

	print *
	call S.print()
	print *

!One iteration: getting njold = (1, 1, ... , 1) and calculating the next vector njnew
call njold.init(z.n)
call njnew.init(z.n)
call delta.init(z.n)
call njold.gen() !njold = (1, 1, ..., 1)
call nday.init(z.n)
call nnight.init(z.n)
call b.init(z.n) !generating RHS for S * njnew = b

!These are the right boundary conditions
!b.d(z.n) = -Fub * h
!b.d(1) = P.d(1) / k.d(1)

Fub = 10


	print *
	print *, "Daytime solution", tau * j
	print *
	
	do i = 2, z.n - 1
		b.d(i) = (njold.d(i) + tau * p.d(i))
	end do

	b.d(z.n) = Fub * h / D.d(z.n - 1)
	b.d(1) = P.d(1) / k.d(1)

	njnew = tridiagonal_matrix_algorithm(S, b)
	
	delta = njnew - njold
	delta_norm = delta.norm()
        njold = njnew

!	print *
!	call njold.print()
	
	
	do while (delta_norm > 1E-5)
	
		j = j + 1
	        do i = 2, z.n - 1
	                b.d(i) = (njold.d(i) + tau * p.d(i)) 
	        end do
	
	        b.d(z.n) = Fub * h / D.d(z.n - 1)
	
		b.d(1) = P.d(1)/ k.d(1)
	
	        njnew = tridiagonal_matrix_algorithm(S, b)

		delta = njnew - njold
		delta_norm = delta.norm()
		njold = njnew

!		print *
!		print *, "Daytime solution", tau * j
!		print *
!		call njold.print

	end do


	nday = njold

	print *
	print *, "Time: 0.00"
	print *
	call nday.print()


!temperature sensitivity

!diurnal evolution

!do j = 1, 288
!	
!	print *
!	print *, "Time:", tau * j
!	print *
!	
!	do i = 2, z.n - 1
!		b.d(i) = (njold.d(i) + tau * p.d(i) * cos(3.141591/288 * j) ** 2) / (1 + tau * k.d(i))
!	end do

!	b.d(z.n) = -Fub * h / D.d(z.n - 1)
!	b.d(1) = P.d(1) * cos(3.141591/288 * j) ** 2 / k.d(1)

!	njnew = tridiagonal_matrix_algorithm(S, b)
!	
!	print *
!	call njnew.print()
!	
!end do




call z.destroy()
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
call njnew.destroy()
call S.destroy()
end
