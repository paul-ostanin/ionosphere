module restore_pwind_mod
implicit none

contains

subroutine restore_pwind(Nx,Ny,Nz,Nx1,Ny1,Nz1,c,lambda_m,eps,pi)

implicit none

integer i,j,k,i1,Nx,Ny,Nz,Nx1,Ny1,Nz1
real eps,pi,s1,s2,lambda_0(1:2),b(1:2),lambda_m(1:Nx+1)
real c(1:Nx,1:Ny,1:Nz,1:3)

!   to define v-component of wind at the pole point 

do k=1,Nz1

!   North Pole

    j=Ny
    s1=sum(c(1:Nx1,j-1,k,2)*cos(lambda_m(1:Nx1)))/pi
    s2=sum(c(1:Nx1,j-1,k,2)*sin(lambda_m(1:Nx1)))/pi
   if (abs(s1)>eps) then
      lambda_0(1)=atan(s2/s1)+pi*(1.-sign(1.,s1))/2.
      b(1)=sqrt(s1*s1+s2*s2)
    else
      lambda_0(1)=pi/2*sign(1.,s2)
      b(1)=abs(s2)
    end if
   c(1:Nx,j,k,2)=b(1)*cos(lambda_m(1:Nx)-lambda_0(1))  
 !  ind(:)=minloc(modulo(abs(lambda_m(1:Nx1)-lambda_0(1)),2*pi))
 !  i_0(1,k)=ind(1)

!  South Pole

   j=1
   s1=sum(c(1:Nx1,j+1,k,2)*cos(lambda_m(1:Nx1)))/pi
   s2=sum(c(1:Nx1,j+1,k,2)*sin(lambda_m(1:Nx1)))/pi
   if (abs(s1)>eps) then
      lambda_0(2)=atan(s2/s1)+pi*(1.-sign(1.,s1))/2.
      b(2)=sqrt(s1*s1+s2*s2) 
    else
      lambda_0(2)=pi/2*sign(1.,s2)
      b(2)=abs(s2)
    end if
   c(1:Nx,j,k,2)=b(2)*cos(lambda_m(1:Nx)-lambda_0(2))
 !   ind(:)=minloc(modulo(abs(lambda_m(1:Nx1)-lambda_0(2)+pi),2*pi))
 !   i_0(2,k)=ind(1)   
end do

end subroutine

end module
