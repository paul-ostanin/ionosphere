program prog2
implicit none

integer, parameter :: N = 400
integer, parameter :: max_size = 1000000
integer, parameter :: max_nonzero = 10000000
integer, parameter :: maxwr = max_nonzero + 8 * max_size
integer, parameter :: maxwi = max_nonzero + 2 * max_size + 1
real*8 , parameter :: DX = 1d0
real*8 , parameter :: DY = 1d0
real*8 , parameter :: DT = 25d-3
real*8 , parameter :: RT = 50d0
real*8 , parameter :: C0 = 1
real*8 , parameter :: SZL = 200
real*8 , parameter :: PI = 3.14159265359

integer ia(max_size + 1), ja(max_nonzero)
real*8 a(max_nonzero), f(max_size), u(max_size), pu(max_size)
real*8 rw(maxwr)
integer iw(maxwi)
real*8 v(max_size)

external matvec, prevec0
integer ITER, INFO, NUNIT
real*8 RESID

integer ierr, ipalu, ipjlu, ipju, ipiw

real*8 ddot
external ddot, dcopy

integer imatvec(1), iprevec(1), ipbcg
real*8 resinit
integer matrix_size, i, j, nonzero

real*8 h
integer sz
integer counter_a, counter_ia
integer counter_rhs
real*8 alpha1, alpha2, alpha3, alpha4, alpha5, beta1
real*8 currt
real*8 C_norm, L2_norm

real*8 prev, analytical_solution, cubature_formula

h = SZL / N
sz = N - 1
pu = 0d0
currt = 0d0

do while (currt .lt. RT)

currt = currt + DT
counter_a = 1
counter_ia = 1
counter_rhs = 1

do i = 1, N - 1
    do j = 1, N - 1

        alpha1 = -1 / (2 * h) - DX / (h * h)
        alpha2 = 1 / DT  + 2 * (DX / (h * h) + DY / (h * h))
        alpha3 = +1 / (2 * h) - DX / (h * h) 
        alpha4 = -DY / (h * h)
        alpha5 = -DY / (h * h)
        beta1 = prev(pu, i, j, N, SZL, C0) / DT

        ia(counter_ia) = counter_a

        if (i .ne. 1) then
            a(counter_a) = alpha1
            ja(counter_a) = (i - 2) * sz + (j - 1) + 1
            counter_a = counter_a + 1
        else if (h * j .le. SZL / 2) then
            beta1 = beta1 - alpha1 * C0
        end if
            
        a(counter_a) = alpha2
        ja(counter_a) = (i - 1) * sz + (j - 1) + 1
        counter_a = counter_a + 1

        if (i .ne. N - 1) then
            a(counter_a) = alpha3
            ja(counter_a) = (i) * sz + (j - 1) + 1
            counter_a = counter_a + 1
        end if
        
        if (j .ne. 1) then
            a(counter_a) = alpha4
            ja(counter_a) = (i - 1) * sz + (j - 2) + 1
            counter_a = counter_a + 1
        end if

        if (j .ne. N - 1) then
            a(counter_a) = alpha5
            ja(counter_a) = (i - 1) * sz + (j) + 1
            counter_a = counter_a + 1
        end if

        f(counter_rhs) = beta1
        counter_ia = counter_ia + 1
        counter_rhs = counter_rhs + 1

    end do
end do

ia(counter_ia) = counter_a

matrix_size = sz * sz
write(*, *) matrix_size, counter_rhs - 1
nonzero = counter_a - 1

ipalu = 1
ipbcg = ipalu + nonzero
ipju = 1
ipjlu = ipju + matrix_size + 1
ipiw = ipjlu + nonzero

call ilu0(matrix_size, a, ja, ia, rw(ipalu), iw(ipjlu), iw(ipju), iw(ipiw), ierr)
write(*, *) ierr

call dcopy(matrix_size, 0d0, 0, u, 1)
resinit = dsqrt(ddot(matrix_size, f, 1, f, 1))

write(*, *) "resinit", resinit

ITER = 1000
RESID = 1d-6 * resinit
INFO = 0
NUNIT = 6
iprevec(1) = matrix_size
imatvec(1) = matrix_size

call slpbcgs(prevec0, iprevec, iw, rw,   &
             matvec, imatvec, ia, ja, a, &
             rw(ipbcg), matrix_size, 8, matrix_size, f, u,   &
             ITER, RESID, INFO, NUNIT)
write(*, *) INFO
pu = u

end do

open(UNIT = 1, FILE = "out.txt")
do j = 1, sz
    do i = 1, sz
        write(1, *) i * h, j * h - SZL / 2, u((i - 1) * sz + (j - 1) + 1)
    end do
	write (1, *)
end do
close(UNIT = 1)

!do i = 1, sz
!    do j = 1, sz
!        v((i - 1) * sz + (j - 1) + 1) = analytical_solution(i * h, j * h - SZL / 2, RT, DX, DY, C0)
!    end do
!end do

!open(UNIT = 2, FILE = "out2.txt")
!do j = 1, sz
!    do i = 1, sz
!        write(2, *) i * h, j * h - SZL / 2, v((i - 1) * sz + (j - 1) + 1)
!    end do
!	write (2, *)
!end do
!close(UNIT = 2)

C_norm = 0d0
do i = 1, sz
    do j = 1, sz
!        C_norm = max(C_norm, abs(u((i - 1) * sz + (j - 1) + 1) - analytical_solution(i * h, j * h - SZL / 2, RT, DX, DY, C0)))
    end do
end do
write(*, *) C_norm

L2_norm = 0
do i = 1, N - 2
    do j = 1, N - 2
!        L2_norm = L2_norm + cubature_formula(i * h, j * h - SZL / 2, &
!                        u((i - 1) * sz + (j - 1) + 1), &
!                        i * h, (j + 1) * h - SZL / 2, &
!                        u((i - 1) * sz + (j) + 1), &
!                        (i + 1) * h, j * h - SZL / 2, &
!                        u((i) * sz + (j - 1) + 1) &
!                        , RT, DX, DY, C0) &
!                        * h * h / 2
!        L2_norm = L2_norm + cubature_formula((i + 1) * h, (j + 1) * h - SZL / 2, &
!                        u((i) * sz + (j) + 1), &
!                        i * h, (j + 1) * h - SZL / 2, &
!                        u((i - 1) * sz + (j) + 1), &
!                        (i + 1) * h, j * h - SZL / 2, &
!                        u((i) * sz + (j - 1) + 1) &
!                        , RT, DX, DY, C0) &
!                        * h * h / 2
    end do
end do
write(*, *) sqrt(L2_norm)

end program prog2

real*8 function prev(pu, i, j, N, SZL, C0)
implicit none
integer :: i, j, N
real*8 SZL, C0
real*8 h
real*8 :: pu(*)
h = SZL / N
if (i .lt. 1 .or. i .gt. N - 1 .or. j .lt. 1 .or. j .gt. N - 1) then
    if (i .eq. 0 .and. h * j .le. SZL / 2) then
        prev = C0
    else
        prev = 0d0
    end if
else
    prev = pu((i - 1) * (N - 1) + (j - 1) + 1)
end if
return
end

real*8 function analytical_solution(x, y, t, DX, DY, C0)
implicit none
real*8 :: x, y, t, DX, DY, C0
real*8 , parameter :: PI = 3.14159265359
real*8 , parameter :: EPS = 1d-8
real*8 integr
real*8 curr, prev, eta
integer N, j

curr = -1d0
prev = 0d0
N = 128
do while (abs(curr - prev) > EPS)
    prev = curr
    curr = 0d0
    N = 2 * N
    do j = 1, N
        eta = cos((2 * j - 1) * PI / (2 * n))
        curr = curr + (1 - eta ** 2) ** 0.5 * integr(x, y, t * (eta + 1) / 2, DX, DY)
    end do
    curr = curr * PI * t / (2 * n)
end do

analytical_solution = x * C0 / (16 * PI * DX) ** 0.5 * curr
end

real*8 function integr(x, y, tau, DX, DY)
implicit none
real*8 x, y, tau, DX, DY
integr = tau ** (-1.5) * (erfc(y / (4 * DY * tau) ** 0.5))* exp(-((x - tau) / (4 * DX * tau) ** 0.5) ** 2)
return
end

real*8 function cubature_formula(x1, y1, val1, x2, y2, val2, x3, y3, val3, t, DX, DY, C0)
implicit none
real*8 t, DX, DY, C0
real*8 x1, x2, x3, y1, y2, y3, val1, val2, val3, wa, wb, a1, a2, a3, b1, b2, b3, pi, x, y, res
real*8 analytical_solution
	
wa = 0.205950504760887
wb = 0.063691414286223
a1 = 0.124949503233232
a2 = 0.437525248383384
a3 = 0.437525248383384
b1 = 0.797112651860071 
b2 = 0.165409927389841
b3 = 0.037477420750088
res = 0

x = a1*x1 + a2*x2 + a3*x3
y = a1*y1 + a2*y2 + a3*y3
res = res + wa*(a1*val1 + a2*val2 + a3*val3 - analytical_solution(x, y, t, DX, DY, C0))**2

x = a2*x1 + a1*x2 + a3*x3
y = a2*y1 + a1*y2 + a3*y3
res = res + wa*(a2*val1 + a1*val2 + a3*val3 - analytical_solution(x, y, t, DX, DY, C0))**2

x = a3*x1 + a2*x2 + a1*x3
y = a3*y1 + a2*y2 + a1*y3
res = res + wa*(a3*val1 + a2*val2 + a1*val3 - analytical_solution(x, y, t, DX, DY, C0))**2

x = b1*x1 + b2*x2 + b3*x3
y = b1*y1 + b2*y2 + b3*y3
res = res + wb*(b1*val1 + b2*val2 + b3*val3 - analytical_solution(x, y, t, DX, DY, C0))**2

x = b1*x1 + b3*x2 + b2*x3
y = b1*y1 + b3*y2 + b2*y3
res = res + wb*(b1*val1 + b3*val2 + b2*val3 - analytical_solution(x, y, t, DX, DY, C0))**2

x = b2*x1 + b1*x2 + b3*x3
y = b2*y1 + b1*y2 + b3*y3
res = res + wb*(b2*val1 + b1*val2 + b3*val3 - analytical_solution(x, y, t, DX, DY, C0))**2

x = b2*x1 + b3*x2 + b1*x3
y = b2*y1 + b3*y2 + b1*y3
res = res + wb*(b2*val1 + b3*val2 + b1*val3 - analytical_solution(x, y, t, DX, DY, C0))**2

x = b3*x1 + b1*x2 + b2*x3
y = b3*y1 + b1*y2 + b2*y3
res = res + wb*(b3*val1 + b1*val2 + b2*val3 - analytical_solution(x, y, t, DX, DY, C0))**2

x = b3*x1 + b2*x2 + b1*x3
y = b3*y1 + b2*y2 + b1*y3
res = res + wb*(b3*val1 + b2*val2 + b1*val3 - analytical_solution(x, y, t, DX, DY, C0))**2


cubature_formula = res
return
end

