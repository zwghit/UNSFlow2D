
!======================================================================================
subroutine limit
implicit none

!call limiter_VKN1
call limiter_VKN
!call limiter_vanAlbada
!call limiter_min_max


end subroutine limit

!======================================================================================

subroutine limiter_VKN1 
use grid 
use commons
implicit none
integer(kind=i4) :: i,j,k,ie,p1,p2,c,cc
real(kind=dp) :: q2, pv(nvar),dist(ndim)
real(kind=dp)    :: x1,y1,x2,y2,dx,dy

real(kind=dp) :: ql(nvar),qr(nvar)
real(kind=dp) :: nr, dr, phi, q_min, q_max,kappa
real(kind=dp) :: TOL,alfa,D_L,D_L0

do i=1,nop
      pt(i)%DUmin(:)=0.d0
      pt(i)%DUmax(:)=0.d0 
      pv(:)=pt(i)%qp(:) 
     
   do j=1,nvar
      q_min=pv(j)
      q_max=pv(j)
      do k=1,pt(i)%nv2v
         c=pt(i)%v2v(k)
         q_min=dmin1(pt(c)%qp(j),q_min)
         q_max=dmax1(pt(c)%qp(j),q_max)
      enddo
      pt(i)%dumax(j)=q_max-pv(j)
      pt(i)%dumin(j)=q_min-pv(j)
   enddo  
enddo

kappa=0.3d0
do i=1,nop
   
   !TOL=(kappa*pt(i)%ds)**2
   TOL=(kappa*dsqrt(pt(i)%cv))**3

   do j=1,nvar 

      do k=1,pt(i)%nv2f
         c=pt(i)%v2f(k)
         p1=fc(c)%pt(1)
         p2=fc(c)%pt(2)
         if(p2==i) then 
            print*,i,p1,p2
            p1=fc(c)%pt(2)
            p2=fc(c)%pt(1)
         endif  
         dist(1)=pt(p2)%x-pt(p1)%x
         dist(2)=pt(p2)%y-pt(p1)%y 
         D_L=sum(pt(i)%grad(:,j)*dist(:))
         !D_L0=dsign(D_L,1.d0)*(dabs(D_L)+eps)
         D_L0=D_L

         if(D_L>eps) then
            alfa=pt(i)%dumax(j)
            if(dabs(alfa)<eps) alfa=0.d0
            nr=(alfa*alfa+TOL)*D_L+2.d0*D_L*D_L*alfa
            dr=alfa*alfa+2.d0*D_L*D_L+alfa*D_L+TOL
            phi=dmin1(1.d0,nr/dr)
         elseif(D_L<-eps) then
            alfa=pt(i)%dumin(j)
            if(dabs(alfa)<eps) alfa=0.d0
            nr=(alfa*alfa+TOL)*D_L+2.d0*D_L*D_L*alfa
            dr=alfa*alfa+2.d0*D_L*D_L+alfa*D_L+TOL
            phi=dmin1(1.d0,nr/dr)
         elseif(D_L==eps) then
         !else
            print*,i,j,D_L
            phi=0.d0
         endif 
      
         pt(i)%phi(j)=dmin1(pt(i)%phi(j),phi)

      enddo  
   enddo  
enddo  

do i=1,nop
   do j=1,nvar 
      pt(i)%grad(1:ndim,j)=pt(i)%grad(1:ndim,j)*pt(i)%phi(j)
      if(pt(i)%grad(1,j)== eps)  stop
      if(pt(i)%grad(2,j)== eps)  stop
   enddo
enddo

end subroutine limiter_VKN1 

!======================================================================================

subroutine limiter_VKN 
use grid 
use commons
implicit none
integer(kind=i4):: i,j,k,c,p1,p2
real(kind=dp)   :: nr,dr,phi,q_min,q_max,kappa
real(kind=dp)   :: TOL,alfa,D_L,var,D_L0
real(kind=dp)   :: pv(nvar),dist(ndim)

do i=1,nop
      pt(i)%DUmin(:)=0.d0
      pt(i)%DUmax(:)=0.d0 
     
   do j=1,nvar
      q_min=1e20
      q_max=-1e20
      do k=1,pt(i)%nv2v
         c=pt(i)%v2v(k)
         q_min=dmin1(pt(c)%qp(j),q_min)
         q_max=dmax1(pt(c)%qp(j),q_max)
      enddo
      q_min=dmin1(pt(i)%qp(j),q_min)
      q_max=dmax1(pt(i)%qp(j),q_max)

      pt(i)%dumax(j)=q_max-pt(i)%qp(j)
      pt(i)%dumin(j)=q_min-pt(i)%qp(j)
   enddo  
enddo

kappa=0.1
do i=1,nop
   !TOL=kappa*pt(i)%cv
   !TOL=(kappa*pt(i)%ds)**2
   TOL=(kappa*dsqrt(pt(i)%cv))**3
   
   do j=1,nvar 
      do k=1,pt(i)%nv2f
         c=pt(i)%v2f(k)
         p1=fc(c)%pt(1)
         p2=fc(c)%pt(2)
         if(p2==i) then 
            p1=fc(c)%pt(2)
            p2=fc(c)%pt(1)
         endif  
         dist(1)=pt(p2)%x-pt(p1)%x
         dist(2)=pt(p2)%y-pt(p1)%y 
         D_L=sum(pt(i)%grad(:,j)*dist(:))
         !D_L=dsign(D_L,1.d0)*(dabs(D_L)+eps)

         !TOL=(kappa*(pt(i)%dumax(j)-pt(i)%dumin(j)))**2
 
         alfa=pt(i)%dumax(j)
         if(D_L<0.d0 ) alfa=pt(i)%dumin(j)
         phi=1.d0 
         if(D_L/=0.d0) then
            if(dabs(alfa)<eps) alfa=0.d0
            nr=alfa*alfa+2.d0*D_L*alfa+TOL
            dr=alfa*alfa+2.d0*D_L*D_L+D_L*alfa+TOL    
            phi=dmin1(1.d0,nr/dr)
            !phi=nr/dr
         endif

         pt(i)%phi(j)=dmin1(pt(i)%phi(j),phi)
         !pt(i)%phi(j)=phi
      enddo  
   enddo  
enddo  

do i=1,nop
   do j=1,nvar 
      pt(i)%grad(1:ndim,j)=pt(i)%grad(1:ndim,j)*pt(i)%phi(j)
   enddo
enddo

end subroutine limiter_VKN 

!======================================================================================

subroutine limiter_vanAlbada
use grid 
use commons
implicit none
integer(kind=i4) :: i,j,k,ie,p1,p2,c
real(kind=dp) :: q2, pv(nvar)
real(kind=dp)    :: x1,y1,x2,y2,dx,dy

real(kind=dp) :: ql(nvar),qr(nvar),dist(ndim)
real(kind=dp) :: nr, dr, phi, q_min, q_max
real(kind=dp) :: alfa,D_L,D_L0,rk,tol,kappa

do i=1,nop
      pt(i)%DUmin(:)=0.d0
      pt(i)%DUmax(:)=0.d0 
     
   do j=1,nvar
      q_min=1e20
      q_max=-1e20
      do k=1,pt(i)%nv2v
         c=pt(i)%v2v(k)
         q_min=dmin1(pt(c)%qp(j),q_min)
         q_max=dmax1(pt(c)%qp(j),q_max)
      enddo
      q_min=dmin1(pt(i)%qp(j),q_min)
      q_max=dmax1(pt(i)%qp(j),q_max)

      pt(i)%dumax(j)=q_max-pt(i)%qp(j)
      pt(i)%dumin(j)=q_min-pt(i)%qp(j)
   enddo  
enddo

kappa=0.3
do i=1,nop
   TOL=(kappa*dsqrt(pt(i)%cv))**3
   
   do j=1,nvar 

      do k=1,pt(i)%nv2f
         c=pt(i)%v2f(k)
         p1=fc(c)%pt(1)
         p2=fc(c)%pt(2)
         if(p2==i) then 
            p1=fc(c)%pt(2)
            p2=fc(c)%pt(1)
         endif  
         dist(1)=pt(p2)%x-pt(p1)%x
         dist(2)=pt(p2)%y-pt(p1)%y 
         D_L=sum(pt(i)%grad(:,j)*dist(:))
         !D_L0=dsign(D_L,1.d0)*(dabs(D_L)+tol) 
         D_L0=D_L

         if(D_L>eps) then
            !rk=pt(i)%dumax(j)/D_L0
            !phi=VanAlbada(rk)
            phi=Albada(pt(i)%dumax(j),D_L0,tol)
         elseif(D_L<-eps) then
            !rk=pt(i)%dumin(j)/D_L0
            !phi=VanAlbada(rk)
            phi=Albada(pt(i)%dumin(j),D_L0,tol)
         else
            phi=1.d0
         endif 
      
         pt(i)%phi(j)=dmin1(pt(i)%phi(j),phi)

      enddo  
   enddo  
enddo  

do i=1,nop
   do j=1,nvar 
      pt(i)%grad(1:ndim,j)=pt(i)%grad(1:ndim,j)*pt(i)%phi(j)
   enddo
enddo

contains
!-----------------------------------------------------------------------------
!real(kind=dp) function VanAlbada(rk)
real(kind=dp) function Albada(a,b,tol)
implicit none 
real(kind=dp) :: rk,a,b,tol


!VanAlbada=(rk*rk+rk)/(1.d0+rk*rk)
!VanAlbada=(rk*rk+2.d0*rk)/(rk*rk+rk+2.d0)
Albada=dmax1(0.d0,(2.d0*a*b+tol*tol)/(a*a+b*b+tol*tol) )


end function Albada

end subroutine limiter_vanAlbada

!-----------------------------------------------------------------------------

subroutine limiter_min_max
use grid 
use commons
implicit none
integer(kind=i4) :: i,j,k,c,p1,p2

real(kind=dp) :: var,dist(ndim)
real(kind=dp) :: phi, q_min, q_max,kappa
real(kind=dp) :: TOL,alfa,D_L

do i=1,nop
      pt(i)%DUmin(:)=0.d0
      pt(i)%DUmax(:)=0.d0 
     
   do j=1,nvar
      q_min=1e20
      q_max=-1e20
      do k=1,pt(i)%nv2v
         c=pt(i)%v2v(k)
         q_min=dmin1(pt(c)%qp(j),q_min)
         q_max=dmax1(pt(c)%qp(j),q_max)
      enddo
      q_min=dmin1(pt(i)%qp(j),q_min)
      q_max=dmax1(pt(i)%qp(j),q_max)

      pt(i)%dumax(j)=q_max-pt(i)%qp(j)
      pt(i)%dumin(j)=q_min-pt(i)%qp(j)
   enddo  
enddo


do i=1,nop
   
   do j=1,nvar 
      q_max=pt(i)%dumax(j)
      q_min=pt(i)%dumin(j)

      do k=1,pt(i)%nv2f
         c=pt(i)%v2f(k)
         p1=fc(c)%pt(1)
         p2=fc(c)%pt(2)
         if(p2==i) then 
            p1=fc(c)%pt(2)
            p2=fc(c)%pt(1)
         endif  
         dist(1)=pt(p2)%x-pt(p1)%x
         dist(2)=pt(p2)%y-pt(p1)%y 
         D_L=sum(pt(i)%grad(:,j)*dist(:))

         if(D_L>0.d0 ) then
            phi=dmin1(1.d0,q_max/D_L)
         elseif(D_L< 0.d0) then
            phi=dmin1(1.d0,q_min/D_L)
         else
            phi=1.d0
         endif 
      
         pt(i)%phi(j)=dmin1(pt(i)%phi(j),phi)

      enddo  
   enddo  
enddo  

do i=1,nop
   do j=1,nvar 
      pt(i)%grad(1:ndim,j)=pt(i)%grad(1:ndim,j)*pt(i)%phi(j)
   enddo
enddo
end subroutine limiter_min_max
