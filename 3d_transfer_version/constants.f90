module constants

integer, parameter::&
&                   Nx = 145, & !37, & !**********
&                   Ny = 91, & !91, & !**********
&                   Nz = 79, & !41, & !**********
&                   Nx1 = Nx-1, &
&                   Ny1 = Ny-1, &
&                   Nz1 = Nz-1, &
&                   Nadv=21 ! number of advection timesteps inside one diffusion timestep
real, parameter::pi = 3.1415926, &
&                a = 6.37e6, &
&                eps = 1e-5, &
&                dx=2*pi/Nx1, &
&                dy=pi/Ny1 
real*8, parameter:: Hmax = 460000

end module
