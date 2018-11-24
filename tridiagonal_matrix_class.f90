module tridi_matrix
use vector

!This class contains a tridiagonal matrix, keeping in memory n * (three nonzero elements in the i-th row)
!The first and the last rows of a tridiagonal matrix contain only 2 nonzero elements => d(1, 1) = d(n, 3) = 0, 
!as they are over the matrix's dimensions


	type tridiagonal_matrix
		integer(4) n
		double precision, allocatable::d(:, :)
	contains
		procedure::init => matrix_init
		procedure::destroy => matrix_destroy
		procedure::print => matrix_print
		procedure::gen => generate_matrix
	end type
	
	interface operator (+)
		module procedure matrix_add
	end interface
	
	interface operator (-)
		module procedure matrix_diff
	end interface

	interface operator  (*)
		module procedure matrix_mul_vec
	end interface
	

contains

subroutine matrix_init (this, n)
!Initializes a tridiagonal matrix (allocates the array). Fills the array with zeros.
	class(tridiagonal_matrix)::this
	integer(4) n
	this.n = n
	allocate (this.d(n, 3))	
	this.d = 0
end subroutine matrix_init

subroutine matrix_print(this)
!Prints out the nonzero elements row by row (except the first and the last rows). 
!The first row starts with a zero element d(1,1), as the last row ends with zero d(n, 3).
	class(tridiagonal_matrix)::this
	integer(4) j
	!write(*, '("Printing matrix:")')
	do j = 1, this.n
		write(*, '(150(e15.7))') this.d(j, :)
	end do
end subroutine matrix_print

function matrix_add(this, m2) result (res)
!Redefined operator "+" for two matrixes
   class(tridiagonal_matrix), intent(in) :: this
   type (tridiagonal_matrix), intent(in) :: m2
   type (tridiagonal_matrix) :: res
   call res.init(m2.n)
   res.d = this.d + m2.d
end function matrix_add


function matrix_diff(this, m2) result (res)
!Redefined operator "+" for two matrixes
   class(tridiagonal_matrix), intent(in) :: this
   type (tridiagonal_matrix), intent(in) :: m2
   type (tridiagonal_matrix) :: res
   call res.init(m2.n)
   res.d = this.d - m2.d
end function matrix_diff


subroutine generate_matrix(this)
!Generates a discrete laplasian matrix with Neumann bounary conditions
	class(tridiagonal_matrix)::this
	integer(4) n, i, j
	n = this.n
	do i = 1, n
		this.d(i, 1) = -1
		this.d(i, 2) = 2
		this.d(i, 3) = -1
	end do
	this.d(1, 1) = 0
	this.d(n, 3) = 0
end subroutine generate_matrix

function matrix_mul_vec(this, r) result (res)
!Redefined operator "*" for a matrix and a vector
	class(tridiagonal_matrix), intent(in) :: this
	type (vect), intent(in) :: r
	type (vect)::res
	integer(4) i
	
	n = this.n
	call res.init(r.n)
	do i = 1, n
		res.d(i) = this.d(i, 1) * r.d(i - 1) + this.d(i, 2) * r.d(i) + this.d(i, 3) * r.d(i + 1)
	end do
end function matrix_mul_vec

function tridiagonal_matrix_algorithm(this, b) result (res)
	class(tridiagonal_matrix), intent(in)::this
	type (vect), intent(in)::b
	type (vect)::res, p, q
	integer(4) i

	call res.init(b.n)
	call p.init(b.n)
	call q.init(b.n)

	p.d(1) = - this.d(1, 3) / this.d(1, 2)
	q.d(1) = b.d(1) / this.d(1, 2)

	do i = 2, b.n - 1
		p.d(i) = this.d(i, 3) / (- this.d(i, 1) * p.d(i - 1) - this.d(i, 2))
		q.d(i) = (- b.d(i) + this.d(i, 1) * q.d(i - 1)) / (- this.d(i, 1) * p.d(i - 1) - this.d(i, 2))
	end do

!	call p.print() 
!	call q.print()
!	do i = 1, b.n
!		write(11,*) l, i, p.d(i)
!	end do
!	write(11,*)

!	do i = 1, b.n
!		write(*,*) l, i, q.d(i)
!	end do
!	write(*,*)
       
	res.d(b.n) = (- b.d(b.n) + this.d(b.n, 1) * q.d(b.n - 1)) / (- this.d(b.n, 1) * p.d(b.n - 1) - this.d(b.n, 2))
	
	i = b.n
	do while (i > 1)
		res.d(i - 1) = res.d(i) * p.d(i - 1) + q.d(i - 1)
		i = i - 1
	end do
end function

subroutine matrix_destroy(this)
	class(tridiagonal_matrix)::this
	deallocate(this.d)
end subroutine matrix_destroy

end module
