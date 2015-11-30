
!======================================================================================
subroutine limit
implicit none

!call limiter_VKN
!call limiter_VKN1
call limiter_min_max
!call limiter_vanAlbada


end subroutine limit

!======================================================================================

subroutine limiter_VKN1 
use grid 
use commons
implicit none
integer(kind=i4) :: i,j,k,ie,in,out,p1,p2,c
real(kind=dp) :: q2, pv(nvar),dist(ndim)
real(kind=dp)    :: x1,y1,x2,y2,dx,dy

real(kind=dp) :: ql(nvar),qr(nvar)
real(kind=dp) :: nr, dr, phi, q_min, q_max,kappa
real(kind=dp) :: TOL,alfa,D_L,D_L0

do i=1,noc
      cell(i)%DUmin(:)=0.d0
      cell(i)%DUmax(:)=0.d0 
      pv(:)=cell(i)%qp(:) 
     
   do j=1,nvar
      q_min=pv(j)
      q_max=pv(j)
      do k=1,cell(i)%nc2c
         c=cell(i)%c2c(k)
         q_min=dmin1(cell(c)%qp(j),q_min)
         q_max=dmax1(cell(c)%qp(j),q_max)
      enddo
      cell(i)%dumax(j)=q_max-pv(j)
      cell(i)%dumin(j)=q_min-pv(j)
   enddo  
enddo

kappa=0.01
do i=1,noc
   
   !TOL=(kappa*cell(in)%ds)**2
   TOL=(kappa*dsqrt(cell(i)%cv))**3

   do j=1,nvar 

      do k=1,cell(i)%nc2f
         c=cell(i)%c2f(k)
         dist(:)=fc(c)%cen(:)-cell(i)%cen(:)
         D_L=sum(cell(i)%grad(:,j)*dist(:))
         D_L0=D_L

         if(D_L>eps) then
            alfa=cell(i)%dumax(j)
            if(dabs(alfa)<eps) alfa=0.d0
            nr=(alfa*alfa+TOL)*D_L+2.d0*D_L*D_L*alfa
            dr=alfa*alfa+2.d0*D_L*D_L+alfa*D_L+TOL
            phi=dmin1(1.d0,nr/dr/D_L0)
         elseif(D_L<-eps) then
            alfa=cell(i)%dumin(j)
            if(dabs(alfa)<eps) alfa=0.d0
            nr=(alfa*alfa+TOL)*D_L+2.d0*D_L*D_L*alfa
            dr=alfa*alfa+2.d0*D_L*D_L+alfa*D_L+TOL
            phi=dmin1(1.d0,nr/dr/D_L0)
         else
            phi=1.d0
         endif 
      
         cell(i)%phi(j)=dmin1(cell(i)%phi(j),phi)

      enddo  
   enddo  
enddo  

do i=1,noc
   do j=1,nvar 
      cell(i)%grad(1:ndim,j)=cell(i)%grad(1:ndim,j)*cell(i)%phi(j)
   enddo
enddo

end subroutine limiter_VKN1 

!======================================================================================

subroutine limiter_VKN 
use grid 
use commons
implicit none
integer(kind=i4):: i,j,k,c
real(kind=dp)   :: nr,dr,phi,q_min,q_max,kappa
real(kind=dp)   :: TOL,alfa,D_L,var,dist(ndim),D_L0

do i=1,noc
   cell(i)%DUmin(:)=0.d0
   cell(i)%DUmax(:)=0.d0 
   do j=1,nvar
      var=cell(i)%qp(j) 
      q_min=var   
      q_max=-var  
      do k=1,cell(i)%nc2c
         c=cell(i)%c2c(k)
         q_min=dmin1(cell(c)%qp(j),q_min)
         q_max=dmax1(cell(c)%qp(j),q_max)
      enddo
      cell(i)%dumax(j)=q_max-var   
      cell(i)%dumin(j)=q_min-var   
   enddo  
enddo

kappa=0.3
do i=1,noc
   TOL=(kappa*dsqrt(cell(i)%cv))**3
   do j=1,nvar 
      do k=1,cell(i)%nc2c
         c=cell(i)%c2c(k)
         dist(:)=fc(c)%cen(:)-cell(i)%cen(:)
         D_L=sum(cell(i)%grad(:,j)*dist(:))
         !D_L0=dsign(D_L,1.d0)*(dabs(D_L)+1e-12)
         D_L0=D_L
 
         if(D_L>eps ) then
            alfa=cell(i)%dumax(j)
            nr=(alfa*alfa+TOL)*D_L+2.d0*D_L*D_L*alfa
            dr=alfa*alfa+2.d0*D_L*D_L+D_L*alfa+TOL    
            phi=nr/dr/D_L0
         elseif(D_L<-eps) then
            alfa=cell(i)%dumin(j)
            nr=(alfa*alfa+TOL)*D_L+2.d0*D_L*D_L*alfa
            dr=alfa*alfa+2.d0*D_L*D_L+D_L*alfa+TOL      
            phi=nr/dr/D_L0
         else
            phi=1.d0   
         endif 
         cell(i)%phi(j)=dmin1(cell(i)%phi(j),phi)
      enddo  
   enddo  
enddo  

do i=1,noc
   do j=1,nvar 
      cell(i)%grad(1:ndim,j)=cell(i)%grad(1:ndim,j)*cell(i)%phi(j)
   enddo
enddo

end subroutine limiter_VKN 

!======================================================================================

subroutine limiter_vanAlbada
use grid 
use commons
implicit none
integer(kind=i4) :: i,j,k,ie,in,out,p1,p2,c
real(kind=dp) :: q2, pv(nvar)
real(kind=dp)    :: x1,y1,x2,y2,dx,dy

real(kind=dp) :: ql(nvar),qr(nvar),dist(ndim)
real(kind=dp) :: nr, dr, phi, q_min, q_max
real(kind=dp) :: alfa,D_L,D_L0,rk

do i=1,noc
      cell(i)%DUmin(:)=0.d0
      cell(i)%DUmax(:)=0.d0 
      pv(:)=cell(i)%qp(:) 
     
   do j=1,nvar
      q_min=pv(j)
      q_max=pv(j)
      do k=1,cell(i)%nc2c
         c=cell(i)%c2c(k)
         q_min=dmin1(cell(c)%qp(j),q_min)
         q_max=dmax1(cell(c)%qp(j),q_max)
      enddo
      cell(i)%dumax(j)=q_max-pv(j)
      cell(i)%dumin(j)=q_min-pv(j)
   enddo  
enddo

do i=1,noc
   
   do j=1,nvar 

      do k=1,cell(i)%nc2c
         c=cell(i)%c2c(k)
         dist(:)=fc(c)%cen(:)-cell(i)%cen(:)
         D_L=sum(cell(i)%grad(:,j)*dist(:))
         D_L0=dsign(D_L,1.d0)*(dabs(D_L)+eps) 
         !D_L0=D_L

         if(D_L>eps) then
            !rk=cell(i)%dumax(j)/D_L0
            !phi=VanAlbada(rk)
            phi=Albada(cell(i)%dumax(j),D_L0)
         elseif(D_L<-eps) then
            !rk=cell(i)%dumin(j)/D_L0
            !phi=VanAlbada(rk)
            phi=Albada(cell(i)%dumin(j),D_L0)
         else
            phi=1.d0
         endif 
      
         cell(i)%phi(j)=dmin1(cell(i)%phi(j),phi)

      enddo  
   enddo  
enddo  

do i=1,noc
   do j=1,nvar 
      cell(i)%grad(1:ndim,j)=cell(i)%grad(1:ndim,j)*cell(i)%phi(j)
   enddo
enddo

contains
!-----------------------------------------------------------------------------
!real(kind=dp) function VanAlbada(rk)
real(kind=dp) function Albada(a,b)
implicit none 
real(kind=dp) :: rk,a,b


!VanAlbada=(rk*rk+rk)/(1.d0+rk*rk)
!VanAlbada=(rk*rk+2.d0*rk)/(rk*rk+rk+2.d0)
Albada=dmax1(0.d0,(2.d0*a*b+eps*eps)/(a*a+b*b+eps*eps) )


end function Albada

end subroutine limiter_vanAlbada

!-----------------------------------------------------------------------------

subroutine limiter_min_max
use grid 
use commons
implicit none
integer(kind=i4) :: i,j,k,c

real(kind=dp) :: var,dist(ndim)
real(kind=dp) :: phi, q_min, q_max,kappa
real(kind=dp) :: TOL,alfa,D_L

do i=1,noc
      cell(i)%DUmin(:)=0.d0
      cell(i)%DUmax(:)=0.d0 
     
   do j=1,nvar
      var=cell(i)%qp(j)
      q_min=var
      q_max=var
      do k=1,cell(i)%nc2c
         c=cell(i)%c2c(k)
         q_min=dmin1(cell(c)%qp(j),q_min)
         q_max=dmax1(cell(c)%qp(j),q_max)
      enddo
      cell(i)%dumax(j)=q_max-var   
      cell(i)%dumin(j)=q_min-var
      !cell(i)%phi(j)=1e20
   enddo  
enddo

do i=1,noc
   
   do j=1,nvar 
      q_max=cell(i)%dumax(j)
      q_min=cell(i)%dumin(j)

      do k=1,cell(i)%nc2c
         c=cell(i)%c2c(k)
         dist(:)=fc(c)%cen(:)-cell(i)%cen(:)
         D_L=sum(cell(i)%grad(:,j)*dist(:))

         if(D_L>eps ) then
            phi=dmin1(1.d0,q_max/D_L)
         elseif(D_L< -eps) then
            phi=dmin1(1.d0,q_min/D_L)
         else
            phi=1.d0
         endif 
      
         cell(i)%phi(j)=dmin1(cell(i)%phi(j),phi)

      enddo  
   enddo  
enddo  

do i=1,noc
   do j=1,nvar 
      cell(i)%grad(1:ndim,j)=cell(i)%grad(1:ndim,j)*cell(i)%phi(j)
   enddo
enddo
end subroutine limiter_min_max
