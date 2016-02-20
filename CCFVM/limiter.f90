
!======================================================================================
subroutine limit
implicit none

call limiter_VKN
!call limiter_vanAlbada
!call limiter_min_max


end subroutine limit

!======================================================================================
subroutine limiter_VKN 
use grid 
use commons
implicit none
integer(kind=i4):: i,j,k,c1,c2
real(kind=dp)   :: nr,dr,phi,kappa,du
real(kind=dp)   :: TOL,alfa,D_L,var,dist(ndim),D_L0

do i=1,noc
   cell(i)%DUmin(:)=eps
   cell(i)%DUmax(:)=-eps
enddo

do i=1,nof
   c1 = fc(i)%in
   c2 = fc(i)%out
   if(c1/=0.and.c2/=0) then 
   do j=1,nvar
      du=cell(c2)%qp(j)-cell(c1)%qp(j)
      cell(c1)%dumax(j)=dmax1(cell(c1)%dumax(j), du)
      cell(c1)%dumin(j)=dmin1(cell(c1)%dumin(j), du)
      cell(c2)%dumax(j)=dmax1(cell(c2)%dumax(j),-du)
      cell(c2)%dumin(j)=dmin1(cell(c2)%dumin(j),-du)
   enddo
   endif
enddo

kappa=0.3

do i=1,nof
   c1 = fc(i)%in
   c2 = fc(i)%out


   if(c1/=0.and.c2/=0) then 
   do j=1,nvar

         !---------------i
         dist(:)=fc(i)%cen(:)-cell(c1)%cen(:)
         D_L=sum(cell(c1)%grad(:,j)*dist(:))

         !TOL=(kappa*(cell(c1)%dumax(j)-cell(c1)%dumin(j)))**2
         TOL=(kappa*dsqrt(cell(c1)%cv))**3

         if(D_L>0.d0) then
            alfa=cell(c1)%dumax(j)
         else
            alfa=cell(c1)%dumin(j)
         endif

         if(dabs(alfa)<eps) alfa=0.d0

         nr=alfa*alfa+2.d0*D_L*alfa+TOL
         dr=alfa*alfa+2.d0*D_L*D_L+D_L*alfa+TOL
         !phi=nr/dr
         phi=dmin1(1.d0,nr/dr)
         cell(c1)%phi(j)=dmin1(cell(c1)%phi(j),phi)

         !---------------j
         dist(:)=fc(i)%cen(:)-cell(c2)%cen(:)
         D_L=sum(cell(c2)%grad(:,j)*dist(:))

         !TOL=(kappa*(cell(c2)%dumax(j)-cell(c2)%dumin(j)))**2
         TOL=(kappa*dsqrt(cell(c2)%cv))**3

         if(D_L>0.d0) then
            alfa=cell(c2)%dumax(j)
         else
            alfa=cell(c2)%dumin(j)
         endif

         if(dabs(alfa)<eps) alfa=0.d0

         nr=alfa*alfa+2.d0*D_L*alfa+TOL
         dr=alfa*alfa+2.d0*D_L*D_L+D_L*alfa+TOL
         !phi=nr/dr
         phi=dmin1(1.d0,nr/dr)
         cell(c2)%phi(j)=dmin1(cell(c2)%phi(j),phi)

   enddo
   endif
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
integer(kind=i4) :: i,j,c1,c2
real(kind=dp) :: q2, pv(nvar)
real(kind=dp)    :: x1,y1,x2,y2,dx,dy

real(kind=dp) :: ql(nvar),qr(nvar),dist(ndim)
real(kind=dp) :: phi,du
real(kind=dp) :: alfa,D_L,D_L0,rk


do i=1,noc
   cell(i)%DUmin(:)=eps
   cell(i)%DUmax(:)=-eps
enddo

do i=1,nof
   c1 = fc(i)%in
   c2 = fc(i)%out
   if(c1/=0.and.c2/=0) then 
   do j=1,nvar
      du=cell(c2)%qp(j)-cell(c1)%qp(j)
      cell(c1)%dumax(j)=dmax1(cell(c1)%dumax(j), du)
      cell(c1)%dumin(j)=dmin1(cell(c1)%dumin(j), du)
      cell(c2)%dumax(j)=dmax1(cell(c2)%dumax(j),-du)
      cell(c2)%dumin(j)=dmin1(cell(c2)%dumin(j),-du)
   enddo
   endif
enddo

do i=1,nof
   c1 = fc(i)%in
   c2 = fc(i)%out

   if(c1/=0.and.c2/=0) then 
   do j=1,nvar

         !---------------i
         dist(:)=fc(i)%cen(:)-cell(c1)%cen(:)
         D_L=sum(cell(c1)%grad(:,j)*dist(:))
         if(D_L>0.d0) then
            phi=Albada(cell(c1)%dumax(j),D_L)
         else
            phi=Albada(cell(c1)%dumin(j),D_L)
         endif

         cell(c1)%phi(j)=dmin1(cell(c1)%phi(j),phi)

         !---------------j
         dist(:)=fc(i)%cen(:)-cell(c2)%cen(:)
         D_L=sum(cell(c2)%grad(:,j)*dist(:))
         !D_L0=dsign(D_L,1.d0)*(dabs(D_L)+eps) 
         if(D_L>0.d0) then
            phi=Albada(cell(c2)%dumax(j),D_L)
         else
            phi=Albada(cell(c2)%dumin(j),D_L)
         endif

         cell(c2)%phi(j)=dmin1(cell(c2)%phi(j),phi)

   enddo
   endif
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
integer(kind=i4) :: i,j,k,c1,c2

real(kind=dp) :: var,dist(ndim)
real(kind=dp) :: phi, kappa
real(kind=dp) :: TOL,alfa,D_L,du

do i=1,noc
   cell(i)%DUmin(:)=eps
   cell(i)%DUmax(:)=-eps
enddo

do i=1,nof
   c1 = fc(i)%in
   c2 = fc(i)%out
   if(c1/=0.and.c2/=0) then 
   do j=1,nvar
      du=cell(c2)%qp(j)-cell(c1)%qp(j)
      cell(c1)%dumax(j)=dmax1(cell(c1)%dumax(j), du)
      cell(c1)%dumin(j)=dmin1(cell(c1)%dumin(j), du)
      cell(c2)%dumax(j)=dmax1(cell(c2)%dumax(j),-du)
      cell(c2)%dumin(j)=dmin1(cell(c2)%dumin(j),-du)
   enddo
   endif
enddo

do i=1,nof
   c1 = fc(i)%in
   c2 = fc(i)%out

   if(c1/=0.and.c2/=0) then 
   do j=1,nvar

         !---------------i
         dist(:)=fc(i)%cen(:)-cell(c1)%cen(:)
         D_L=sum(cell(c1)%grad(:,j)*dist(:))
         D_L0=dsign(D_L,1.d0)*(dabs(D_L)+eps) 
         if(D_L>0.d0) then
            phi=dmin1(1.d0,cell(c1)%dumax(j) /D_L0)
         else
            phi=dmin1(1.d0,cell(c1)%dumin(j) /D_L0)
         endif

         cell(c1)%phi(j)=dmin1(cell(c1)%phi(j),phi)

         !---------------j
         dist(:)=fc(i)%cen(:)-cell(c2)%cen(:)
         D_L=sum(cell(c2)%grad(:,j)*dist(:))
         D_L0=dsign(D_L,1.d0)*(dabs(D_L)+eps) 
         if(D_L>0.d0) then
            phi=dmin1(1.d0,cell(c2)%dumax(j) /D_L0)
         else
            phi=dmin1(1.d0,cell(c2)%dumin(j) /D_L0)
         endif

         cell(c2)%phi(j)=dmin1(cell(c2)%phi(j),phi)

   enddo
   endif
enddo


do i=1,noc
   do j=1,nvar 
      cell(i)%grad(1:ndim,j)=cell(i)%grad(1:ndim,j)*cell(i)%phi(j)
   enddo
enddo
end subroutine limiter_min_max
