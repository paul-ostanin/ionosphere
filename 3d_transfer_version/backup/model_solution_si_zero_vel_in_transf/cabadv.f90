module cabadv_mod
implicit none

contains

subroutine cab_adv(Nx,Ny,Nz,Nx1,Ny1,Nz1,c,u,fu,ux,uy,uz,u1,u2,dlam,dphi,dlev,phi_m,cphi,cphi_m,lambda,lambda_m,z_m,a,pi,dt)

! Solve mass continiuty equation on one time step with help of two-level "cabaret" scheme
! predictor-corrector variant of the scheme,
! 3d spherical layer case, C-grid,
! uniform grid in longitude, latitude and nonuniform in height
! limiter for monotonicity based on maximum principle 

! Written by S.Kostrykin (25/02/2010)

! DESCRIPTION OF INPUT-OUTPUT VARIABLES:

! Nx,Ny,Nz (input) - dimensions of velocity components array in x,y,z-directions, 
!    suggesting following boundary conditions: periodicity in x  - u,v,w(1,:,:)=u,v,w(Nx,:,:), on y - v(:,1,:)-South Pole, v(:,Ny,:)-North Pole, on z - w(:,:1)=w(:,:,Nz)=0
! Nx1=Nx-1,Ny1=Ny-1,Nz1=Nz-1 (input)
! c (input) - array with velocity components (in m/s) defined on C-grid at the centers of cell boundaries: c(:,:,:,1,:) - u - component, c(:,:,:,2,:) - v - component, c(:,:,:,3,:) - w - component, 
!    c(:,:,:,:,1) - on t=n*dt , c(:,:,:,:,2) - on t=(n+1)*dt 
! u - array with tracer density (in kg/m^3)  defined on C-grid in the cell centers: u(:,:,:,1) - on t=n*dt (input), u(:,:,:,2) - on t=(n+1)*dt (output)
! ux,uy,uz - arrays with tracer fluxes defined at the positions of u,v,w - velocity components:  ux,uy,uz(:,:,:,1) - on t=n*dt (input), ux,ux,uz(:,:,:,2) - on t=(n+1)*dt (output),
!    (at the first time-step these arrays are defined by interpolation from u(:,:,:,1))
! fu (input) - array with non-advective tendencies from mass continuity equation
! dlam,dphi,dlev (input) - grid sizes in x,y,z - dimensions: dlam,dphi - in radians, dlev - in meters
! lambda_m, phi_m, z_m (input) - coordinates of grid boundaries at x,y,z - dimensions: lambda_m, phi_m - in radians, z_m - in meters
! cphi_m (input) - cosine of latitudes of grid boundaries
! lambda (input) - longitudes of cell centers (in radians)
! cphi (input) - cosine of latitudes of cell centers
! pi (input) - pi-number
! dt (input) - time step (in seconds)
! a (input) - sphere radius (in meters)
! u1,u2 (input) - arrays of boundary values at the upper and lower boundaries

implicit none

logical,parameter::fcorr=.true.
integer Nx,Ny,Nz,Nx1,Ny1,Nz1,i,j,k,i1,p,q,r,ip1,ip2,jq1,jq2,kr1,kr2,itr,itrr(1),i_0(1:2,1:Nz1)
real pi,a,uxmin,uxmax,uymin,uymax,uzmin,uzmax,dt,up1,up2,ftx,fty,ftz
real c(1:Nx,1:Ny,1:Nz,1:3,1:2),u(1:Nx1,1:Ny1,1:Nz1,1:2),fu(1:Nx1,1:Ny1,1:Nz1),ux(1:Nx,1:Ny1,1:Nz1,1:2),uy(1:Nx,1:Ny,1:Nz1,1:2),uz(1:Nx,1:Ny1,1:Nz,1:2),&
&      dphi(1:Ny1),dlev(1:Nz1),dlam(1:Nx),lambda_0(2),lambda(1:Nx),lambda_m(1:Nx+1),phi_m(1:Ny),z_m(1:Nz),cphi_m(1:Ny),cphi(1:Ny1),&
&      u1(1:Nx1,1:Ny1),u2(1:Nx1,1:Ny1)

! predictor step

do k=1,Nz1
  do i=1,Nx1 
    i1=modulo(i,Nx1)+1
    do j=1,Ny1     
      u(i,j,k,2)=u(i,j,k,1)-0.5*dt*(((c(i1,j,k,1,1)*ux(i1,j,k,1)-c(i,j,k,1,1)*ux(i,j,k,1))/dlam(1)&
&                                  +(c(i,j+1,k,2,1)*uy(i,j+1,k,1)*cphi_m(j+1)-c(i,j,k,2,1)*uy(i,j,k,1)*cphi_m(j))/dphi(1))/(a*cphi(j))&
&                                  +((c(i,j,k+1,3,1)*uz(i,j,k+1,1)-c(i,j,k,3,1)*uz(i,j,k,1)))/dlev(k))+0.5*fu(i,j,k)*dt
    end do
  end do
end do

! flux variables update

do k=1,Nz1
  do i=1,Nx1
    do j=1,Ny1
      p=(sign(1.,c(i,j,k,1,2))-1)/2
      ip1=modulo(i-2-p,Nx1)+1
      ip2=modulo(i-2-2*p,Nx1)+1   
      ux(i,j,k,2)=2*u(ip1,j,k,2)-ux(ip2,j,k,1)
     if (fcorr) then
         ftx=2*u(ip1,j,k,2)-2*u(ip1,j,k,1)+0.5*dt/a*(c(i,j,k,1,1)+c(ip2,j,k,1,1))*(ux(i,j,k,1)-ux(ip2,j,k,1))/(cphi(j)*dlam(1)*(2*p+1))
      else
        ftx=0.
      end if
        uxmin=min(ux(i,j,k,1)+ftx,ux(ip2,j,k,1)+ftx,u(ip1,j,k,2)+ftx/2)
        uxmax=max(ux(i,j,k,1)+ftx,ux(ip2,j,k,1)+ftx,u(ip1,j,k,2)+ftx/2)
        ux(i,j,k,2)=max(0.,min(max(ux(i,j,k,2),uxmin),uxmax))
    end do
    do j=2,Ny1
      q=(sign(1.,c(i,j,k,2,2))-1)/2
      jq1=j-q-1
      jq2=j-2*q-1
      uy(i,j,k,2)=2*u(i,jq1,k,2)-uy(i,jq2,k,1)
      if (fcorr) then
         fty=2*u(i,jq1,k,2)-2*u(i,jq1,k,1)+0.5*dt/a*(c(i,j,k,2,1)+c(i,jq2,k,2,1))*(uy(i,j,k,1)-uy(i,jq2,k,1))/(phi_m(j)-phi_m(jq2))
      else
        fty=0.
      endif
        uymin=min(uy(i,j,k,1)+fty,uy(i,jq2,k,1)+fty,u(i,jq1,k,2)+fty/2)
        uymax=max(uy(i,j,k,1)+fty,uy(i,jq2,k,1)+fty,u(i,jq1,k,2)+fty/2)
        uy(i,j,k,2)=max(0.,min(max(uy(i,j,k,2),uymin),uymax))
    end  do
  end do

  ! to define flux variables at the pole points
  ! North Pole

  j=Ny
 itrr(:)=maxloc(c(1:Nx1,j,k,2,1))
 i=itrr(1)
 q=0
  jq1=j-q-1
  jq2=j-2*q-1
  up1=2*u(i,j-1,k,2)-uy(i,j-1,k,1)  
  if (fcorr) then
      fty=2*u(i,jq1,k,2)-2*u(i,jq1,k,1)+0.5*dt/a*(c(i,j,k,2,1)+c(i,jq2,k,2,1))*(uy(i,j,k,1)-uy(i,jq2,k,1))/(phi_m(j)-phi_m(jq2))
  else
    fty=0.
  end if
    uymin=min(uy(i,j,k,1)+fty,uy(i,j-1,k,1)+fty,u(i,j-1,k,2)+fty/2)
    uymax=max(uy(i,j,k,1)+fty,uy(i,j-1,k,1)+fty,u(i,j-1,k,2)+fty/2)
    uy(1:Nx1,j,k,2)=max(0.,min(max(up1,uymin),uymax))

 ! South Pole
  j=1
  itrr(:)=minloc(c(1:Nx1,j,k,2,1))
  i=itrr(1)
  q=-1
  jq1=j-q-1
  jq2=j-2*q-1
  up2=2*u(i,j,k,2)-uy(i,j+1,k,1)
  if (fcorr) then
      fty=2*u(i,jq1,k,2)-2*u(i,jq1,k,1)+0.5*dt/a*(c(i,j,k,2,1)+c(i,jq2,k,2,1))*(uy(i,j,k,1)-uy(i,jq2,k,1))/(phi_m(j)-phi_m(jq2))
  else
    fty=0.
  end if
    uymin=min(uy(i,j,k,1)+fty,uy(i,j+1,k,1)+fty,u(i,j,k,2)+fty/2)
    uymax=max(uy(i,j,k,1)+fty,uy(i,j+1,k,1)+fty,u(i,j,k,2)+fty/2)
    uy(1:Nx1,j,k,2)=max(0.,min(max(up2,uymin),uymax))
end do

do k=2,Nz1
  do i=1,Nx1
    do j=1,Ny1
      r=(sign(1.,-c(i,j,k,3,2))-1)/2
      kr1=k-r-1
      kr2=k-2*r-1
      uz(i,j,k,2)=2*u(i,j,kr1,2)-uz(i,j,kr2,1)
      if (fcorr) then
         ftz=2*u(i,j,kr1,2)-2*u(i,j,kr1,1)+0.5*dt*(c(i,j,k,3,1)+c(i,j,kr2,3,1))*(uz(i,j,k,1)-uz(i,j,kr2,1))/(z_m(k)-z_m(kr2))
      else
        ftz=0.
       end if
        uzmin=min(uz(i,j,k,1)+ftz,uz(i,j,kr2,1)+ftz,u(i,j,kr1,2)+ftz/2)
        uzmax=max(uz(i,j,k,1)+ftz,uz(i,j,kr2,1)+ftz,u(i,j,kr1,2)+ftz/2)
        uz(i,j,k,2)=max(0.,min(max(uz(i,j,k,2),uzmin),uzmax))
    end do
  end do
end do
! to define flux variables at upper and lower boundaries
! upper boundary
k=1
  do i=1,Nx1
    do j=1,Ny1
      r=(sign(1.,-c(i,j,k,3,2))-1)/2
      if (r==-1) then
        kr1=k-r-1
        kr2=k-2*r-1
        uz(i,j,k,2)=2*u(i,j,kr1,2)-uz(i,j,kr2,1)
        if (fcorr) then
           ftz=2*u(i,j,kr1,2)-2*u(i,j,kr1,1)+0.5*dt*(c(i,j,k,3,1)+c(i,j,kr2,3,1))*(uz(i,j,k,1)-uz(i,j,kr2,1))/(z_m(k)-z_m(kr2))
        else
          ftz=0.
        end if
        uzmin=min(uz(i,j,k,1)+ftz,uz(i,j,kr2,1)+ftz,u(i,j,kr1,2)+ftz/2)
        uzmax=max(uz(i,j,k,1)+ftz,uz(i,j,kr2,1)+ftz,u(i,j,kr1,2)+ftz/2)
        uz(i,j,k,2)=max(0.,min(max(uz(i,j,k,2),uzmin),uzmax))
      else
!        uz(i,j,k,2)=u(i,j,k,2)
        uz(i,j,k,2)=u1(i,j)
      end if
    end do
  end do
! lower boundary
k=Nz
  do i=1,Nx1
    do j=1,Ny1
      r=(sign(1.,-c(i,j,k,3,2))-1)/2
      if (r==0) then
        kr1=k-r-1
        kr2=k-2*r-1
        uz(i,j,k,2)=2*u(i,j,kr1,2)-uz(i,j,kr2,1)
        if (fcorr) then
           ftz=2*u(i,j,kr1,2)-2*u(i,j,kr1,1)+0.5*dt*(c(i,j,k,3,1)+c(i,j,kr2,3,1))*(uz(i,j,k,1)-uz(i,j,kr2,1))/(z_m(k)-z_m(kr2))
        else
          ftz=0.
        end if
        uzmin=min(uz(i,j,k,1)+ftz,uz(i,j,kr2,1)+ftz,u(i,j,kr1,2)+ftz/2)
        uzmax=max(uz(i,j,k,1)+ftz,uz(i,j,kr2,1)+ftz,u(i,j,kr1,2)+ftz/2)
        uz(i,j,k,2)=max(0.,min(max(uz(i,j,k,2),uzmin),uzmax))
      else
        uz(i,j,k,2)=u2(i,j)
      end if
    end do
  end do


! corrector step

do k=1,Nz1
  do i=1,Nx1
    do j=1,Ny1
      i1=modulo(i,Nx1)+1
      u(i,j,k,2)=u(i,j,k,2)-0.5*dt*(((c(i1,j,k,1,2)*ux(i1,j,k,2)-c(i,j,k,1,2)*ux(i,j,k,2))/dlam(1)&
&                                  +(c(i,j+1,k,2,2)*uy(i,j+1,k,2)*cphi_m(j+1)-c(i,j,k,2,2)*uy(i,j,k,2)*cphi_m(j))/dphi(1))/(a*cphi(j))&
&                                  +(c(i,j,k+1,3,2)*uz(i,j,k+1,2)-c(i,j,k,3,2)*uz(i,j,k,2))/dlev(k))+0.5*fu(i,j,k)*dt
    end do
  end do
end do

end subroutine

end module
