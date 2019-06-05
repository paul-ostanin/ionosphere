module functions
implicit none

contains

real*8 function dA_dphi(phi)
    implicit none
    real*8 phi 
    dA_dphi = -8*tan(phi) / cos(phi) / (1+4*(tan(phi)**2))**2 - sin(phi) / (1+4*(tan(phi)**2))
    return
end

real*8 function dB_dphi(phi)
    implicit none
    real*8 phi
    dB_dphi =  -32*tan(phi)**2 / cos(phi) / (1+4*(tan(phi)**2))**2 + 4*cos(phi) / (1+4*(tan(phi)**2))
    return
end

real*8 function f_m(z, phi)
    implicit none
    real*8 phi, z
    f_m =  3E+6 * ((z/1000)-100)/133 * exp((100-(z/1000))/133) * cos(phi/2)**2
    return
end

real*8 function df_dz(z, phi)
    implicit none
    real*8 phi, z
    df_dz =  3E+6 * exp((100-(z/1000))/133) * cos(phi/2)**2 * (1 - ((z/1000)-100)/133) / 133  * 1E-3
    return
end

real*8 function df_dphi(z, phi)
    implicit none
    real*8 phi, z
    df_dphi =  3E+6 * ((z/1000)-100)/133 * exp((100-(z/1000))/133) * cos(phi/2)**2 * (-tan(phi/2))
    return
end

real*8 function ddf_dzdz(z, phi)
    implicit none
    real*8 phi, z
    ddf_dzdz =  1E-6 * 3E+6 * cos(phi/2)**2 * &
            (-exp((100-(z/1000))/133) - &
            &exp((100-(z/1000))/133) * (1 - ((z/1000)-100)/133) )/133/133
    return
end

real*8 function ddf_dphidphi(z, phi)
    implicit none
    real*8 phi, z
    ddf_dphidphi = 3E+6 * ((z/1000)-100)/133 * exp((100-(z/1000))/133) * (-cos(phi)/2)
    !ddf_dphidphi = -2*cos(z)**2 * sin(2*phi)
    return
end

real*8 function ddf_dzdphi(z, phi)
    implicit none
    real*8 phi, z
    ddf_dzdphi = 3E+6 * exp((100-(z/1000))/133) * cos(phi/2)**2 * (1 - ((z/1000)-100)/133) * (-tan(phi/2)) / 133 * 1E-3
    !ddf_dzdphi = sin(2*z)*sin(2*phi)
    return
end


real*8 function vel_z(z, phi) !v_z = d(wind_d)/dz
    implicit none
    real*8 phi, z, pi
    pi = 3.141592653589793238462643
    vel_z =  -0.7 * 1E-1*cos((z/1000-100)/400 * pi - pi/2) * 50000

    return
end

real*8 function dvelz_dz(z, phi) !dv_z/dz
    implicit none
    real*8 phi, z, pi
    pi = 3.141592653589793238462643
    dvelz_dz =  +0.7 * 1E-1*sin((z/1000-100)/400 * pi - pi/2) * pi/400 * 1E-3 * 50000
    return
end

real*8 function vel_phi(z, phi) !v_phi = d(wind_d)/dphi
    implicit none
    real*8 phi, z, pi
    pi = 3.141592653589793238462643
    vel_phi = 0
    return
end

real*8 function dvelphi_dphi(z, phi) !dv_phi/dphi
    implicit none
    real*8 phi, z, pi
    pi = 3.141592653589793238462643
    dvelphi_dphi =  0
    return
end

end module

