module constants

character(len = *), parameter::maska = 'ion_upb_trd2'
integer, parameter::&
&                   sp = 3, &
&                   ic = 20, & 
&                   sf = 2**(sp-1), &
&                   Nx = 145, &
&                   Ny = 91, &
&                   Nz = 81, &
&                   Nout = 24, &
&                   Nt = Nout*sf*ic, &
&                   Nx1 = Nx-1, &
&                   Ny1 = Ny-1, &
&                   Nz1 = Nz-1, &
&                   Nnorms = 5
real, parameter::pi = 3.1415926, &
&                tend = 1.*86400, &
&                h_atm = 1.0e5, &
&                a = 6.37e6, &
&                dx = 2*pi/Nx1, &
&                dy = pi/Ny1, &
&                dz = h_atm/Nz1, &
&                dt = tend/Nt, &
&                u0 = 2*pi*a/tend, &
&                N_vert_kol = 1., &
&                Omega = 2*pi*N_vert_kol/tend, &
&                hraz = 0.1*h_atm, &
&                w0 = hraz*Omega, &
&                lambda_c = -80.*pi/180., &
&                phi_c = 0.*pi/180., &
&                h0 = 0.5*h_atm, &!0.1*h_atm, &
&                hb = h_atm, &
&                z_c = 2.5*h_atm, &
&                eps = 1e-5, &
&                alpha = 0.*pi/180

end module