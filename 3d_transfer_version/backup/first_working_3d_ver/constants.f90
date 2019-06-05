module constants

integer, parameter::&
&                   Nx = 145, & !37, & !**********
&                   Ny = 91, & !91, & !**********
&                   Nz = 81, & !41, & !**********
&                   Nx1 = Nx-1, &
&                   Ny1 = Ny-1, &
&                   Nz1 = Nz-1
real, parameter::pi = 3.1415926, &
&                a = 6.37e6, &
&                eps = 1e-5, &
&                dx=2*pi/Nx1, &
&                dy=pi/Ny1 

end module
