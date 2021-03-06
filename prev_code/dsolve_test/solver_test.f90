program mat

use tridi_matrix
use vector

implicit none

type (tridiagonal_matrix) S
type (vect) b, z, h, hmid, res, m, njold, njnew, nday, nnight, delta
integer i, j
real (8) tau, h0, Fub, delta_norm, eps, alpha

!Time step (in seconds) 5 min
tau = 1E-1

!Space step (in cm) 5 km
h0 = 1E+0/2000

!Expected solution: u = sin (alpha z) = sin (pi z)
alpha = 3.14159

!Vector of altitudes. Step h_i = z(i) - z(i - 1)
call z.init(2001)
do i = 1, z.n
	z.d(i) = h0 * (i - 1)
end do

!call z.print()

!Vector of middles of [z(i), z(i+1)].  m.d(i) = z_{i+1/2} 
call m.init(z.n - 1)
do i = 1, z.n - 1
	m.d(i) = 0.5 * (z.d(i) + z.d(i + 1))
end do

call h.init(z.n - 1)
do i = 1, z.n - 1
	h.d(i) = z.d(i + 1) - z.d(i)
end do

call hmid.init(z.n - 1)
do i = 1, z.n - 1
	hmid.d(i) = m.d(i) - m.d(i - 1)
end do


!System matrix. 
call S.init(z.n)
!lower boundary condition: n_1 = P_1/k_1
S.d(1, 2) = 1




!new upper boundary condition:
S.d(z.n, 1) =    - 1/(h.d(z.n-1)) - 0.5/h.d(z.n-1)
S.d(z.n, 2) = h.d(z.n-1)/tau + 1/h.d(z.n-1) - 0.5 - alpha**2 * h.d(z.n-1)

!middle matrix rows are similar
do i = 2, z.n - 1

! old scheme

!	S.d(i, 1) = -D.d(i - 1) * tau / (h0 * h0)
!	S.d(i, 2) = (1 + tau * k.d(i)) + u.d(i) * tau / h0 + (D.d(i - 1) + D.d(i))* tau / (h0 * h0) 
!	S.d(i, 3) = -D.d(i)* tau / (h0 * h0) - u.d(i + 1) * tau / h0


! symmetric scheme

	S.d(i, 1) = -1/(hmid.d(i) * h.d(i-1)) + 1/(h.d(i) + h.d(i-1))
	S.d(i, 2) = 1/tau - alpha**2 + (1/h.d(i-1)+1/h.d(i))/ hmid.d(i)
	S.d(i, 3) = -1/(hmid.d(i) * h.d( i )) - 1/(h.d(i) + h.d(i-1))

end do

!	print *
!	call S.print()
!	print *


call njold.init(z.n)
call njnew.init(z.n)
call delta.init(z.n)
call njold.gen() !njold = (1, 1, ..., 1)
call b.init(z.n) !generating RHS for S * njnew = b

Fub = alpha


!One iteration: setting njold = (1, 1, ... , 1) and calculating the next vector njnew
!One iteration: setting njold = (1, 1, ... , 1) and calculating the next vector njnew
	print *
	print *, "Daytime solution", tau * j
	print *

	!Setting RHS in the middle	
	do i = 2, z.n - 1
		b.d(i) = njold.d(i)/tau - alpha * cos(alpha * z.d(i))
	end do

	!Setting boundary conditions for the RHS
	b.d(z.n) = njold.d(z.n)/tau*h.d(z.n-1) - alpha + alpha*h.d(z.n-1)/2 
	b.d(1) = 0

	!Solving the system
	njnew = tridiagonal_matrix_algorithm(S, b)
	
	delta = njnew - njold
	delta_norm = delta.norm()
        njold = njnew

!	print *
!	call njold.print()

j = 0

!Next iterations, until we get the stationarity point	
	do while (delta_norm > 1E-5)
	
		j = j + 1
		!Setting RHS in the middle
	        do i = 2, z.n - 1
	                b.d(i) = njold.d(i)/tau - alpha * cos(alpha * z.d(i))
	        end do

		!Setting boundary conditions for the RHS
	        b.d(z.n) = njold.d(z.n)/tau*h.d(z.n-1) - alpha + alpha*h.d(z.n-1)/2 
		b.d(1) = 0

		!Solving the system
	        njnew = tridiagonal_matrix_algorithm(S, b)

		delta = njnew - njold
		delta_norm = delta.norm()
		njold = njnew

!		print *
!		print *, "Daytime solution", j
!		print *
!		call njold.print

	end do


        open(unit=10, name='res.txt')
	call njold.print_column(10)
	close(unit=10)




call z.destroy()
call njold.destroy()
call njnew.destroy()
call S.destroy()
end
