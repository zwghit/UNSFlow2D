!======================================================================================
subroutine Gradient_LSQR
use grid 
use commons
implicit none
integer(kind=i4) :: i,j,k,c
real(kind=dp)    :: xc,yc,dx,dy,wt
real(kind=dp)    :: wx,wy     
real(kind=dp)    :: r11,r12,r22,alfa1,alfa2


do i=1,nop
   pt(i)%grad(:,:)=0.d0
enddo

do i=1,nop
   xc=pt(i)%x
   yc=pt(i)%y

   do k=1,nvar
      do j=1,pt(i)%nv2v
         c=pt(i)%v2v(j)
         pt(i)%grad(1,k)=pt(i)%grad(1,k)+(pt(c)%qp(k)-pt(i)%qp(k))*pt(i)%wx(j)
         pt(i)%grad(2,k)=pt(i)%grad(2,k)+(pt(c)%qp(k)-pt(i)%qp(k))*pt(i)%wy(j)
      enddo
   enddo

enddo

call limit

end subroutine Gradient_LSQR

!======================================================================================

!======================================================================================
subroutine Gradient_GG
use grid 
use commons
implicit none
integer(kind=i4) :: i,j,k,in,out,c,p1,p2,ps
real(kind=dp) :: var,dx,dy,xc,yc,ds 
real(kind=dp) :: ql(nvar),qr(nvar)
real(kind=dp) :: nr, dr, phi, q_min, q_max,kappa
real(kind=dp) :: TOL,alfa,D_L,D_L0,check,dotprod

do i=1,nop
   pt(i)%grad(:,:)=0.d0
enddo

do i=1,nop

   do j=1,pt(i)%nv2f
      c=pt(i)%v2f(j)
      p1=fc(c)%pt(1)
      p2=fc(c)%pt(2)
      ds=fc(c)%area
      dx=fc(c)%nx*ds
      dy=fc(c)%ny*ds
      do k=1,nvar
         var=0.5d0*(pt(p1)%qp(k)+pt(p2)%qp(k))
         pt(i)%grad(1,k)=pt(i)%grad(1,k)+var*dx
         pt(i)%grad(2,k)=pt(i)%grad(2,k)+var*dy
      enddo
   enddo

enddo

do i=1,nop
   var=pt(i)%cv
   do j=1,nvar
   pt(i)%grad(1,j)=pt(i)%grad(1,j)/var
   pt(i)%grad(2,j)=pt(i)%grad(2,j)/var
   enddo
enddo

!return

call limit

end subroutine Gradient_GG

!======================================================================================
