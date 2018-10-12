module vector

	type vect 
		integer(4) n
		double precision, allocatable::d(:)
	contains
		procedure::init => vec_init
		procedure::destroy => vec_destroy
		procedure::print => vec_print
		procedure::print_stdout => vec_print_old
		procedure::print_column => vec_print_column
		procedure::print_column_with_z => vec_print_column_with_z
		procedure::gen => generate_vector
		procedure::genzero => generate_zero
		procedure::norm => norm
		procedure::interp => interpolation
	end type

	interface operator (-)
		module procedure vec_diff
	end interface


	interface operator (+)
		module procedure vec_sum
	end interface

        interface operator (*)
                module procedure vec_mul_num
        end interface


contains

subroutine vec_init (this, n)
	class(vect)::this
	integer(4) n
	this.n = n
	allocate (this.d(n))
	this.d = 0.0_8
end subroutine vec_init

subroutine vec_print_old(this)
	class(vect)::this
	integer(4) j
	write(*, '(1000(e13.3))') (this.d(j), j = 1, this.n)
end subroutine vec_print_old

subroutine vec_print(this, descriptor)
	class(vect)::this
	integer(4) j, descriptor
!	write(*, '("Printing vector:")')
	write(descriptor, '(1000(e13.3))') (this.d(j), j = 1, this.n)
end subroutine vec_print

subroutine vec_print_column_with_z(this, descriptor)
        class(vect)::this
        integer(4) j, descriptor
!       write(*, '("Printing vector:")')
        write(descriptor, '(e14.7, e14.7)') ((100+4E2*j/this.n), this.d(j), j = 1, this.n)
end subroutine vec_print_column_with_z

subroutine vec_print_column(this, descriptor)
        class(vect)::this
        integer(4) j, descriptor
!       write(*, '("Printing vector:")')
        write(descriptor, '(e14.7)') (this.d(j), j = 1, this.n)
end subroutine vec_print_column

subroutine generate_vector(this)
	class(vect)::this
	integer(4) n, i
	n = this.n
	do i = 1, n
		this.d(i) = 1
	end do
end subroutine generate_vector

subroutine generate_zero(this)
        class(vect)::this
        integer(4) n, i
        n = this.n
        do i = 1, n
                this.d(i) = 0
        end do
end subroutine generate_zero


function norm(this) result(res)
	class(vect)::this
	integer(4) n, i
	double precision:: res
	res = 0
	n = this.n
	do i = 1, n
		res = res + this.d(i) * this.d(i)
	end do
	res = sqrt(res)
end function

function vec_sum(this, r) result (res)
   class(vect), intent(in) :: this
   type (vect), intent(in) :: r
   type (vect) :: res
   call res.init(r.n)
   res.d = this.d + r.d
end function vec_sum

function vec_diff(this, r) result (res)
   class(vect), intent(in) :: this
   type (vect), intent(in) :: r
   type (vect) :: res
   call res.init(r.n)
   res.d = this.d - r.d
end function vec_diff

function vec_mul_num(this, r) result (res)
   class(vect), intent(in) :: this
   double precision, intent(in) :: r
   type (vect) :: res
   call res.init(this.n)
   res.d = this.d * r
end function vec_mul_num

!Vector contains a set of values f(z1), f(z2),...,f(zn).
!Altitudes z1, z2,...,zn are taken from vector z
!This function returns a value at z=z0 by calculating a linear interpolation f(z0) = f(z_{i-1}) + (z_i-z0)/(z_i-z_{i-1})(f(z_i)-f(z_{i-1}) 
function interpolation(this, z, z0) result (res)
   class(vect), intent(in) :: this
   type (vect), intent(in) :: z
   double precision z0, res
   integer(4) i
	i = 1
	do while (z.d(i) < z0)
		i = i + 1
	end do
	res = this.d(i - 1) + (z.d(i) - z0)/(z.d(i) - z.d(i - 1)) * (this.d(i) - this.d(i - 1))
end function interpolation

subroutine vec_destroy(this)
	class(vect)::this
	deallocate(this.d)
end subroutine vec_destroy

end module
