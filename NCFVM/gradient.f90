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
   r11=pt(i)%r11
   r12=pt(i)%r12
   r22=pt(i)%r22

  do k=1,nvar
   do j=1,pt(i)%nv2v
      c=pt(i)%v2v(j)
      dx=pt(c)%x-xc
      dy=pt(c)%y-yc

      alfa1=dx/r11/r11
      alfa2=(dy-dx*r12/r11)/r22/r22
      wx=alfa1-alfa2*r12/r11
      wy=alfa2

      pt(i)%grad(1,k)=pt(i)%grad(1,k)+(pt(c)%qp(k)-pt(i)%qp(k))*wx
      pt(i)%grad(2,k)=pt(i)%grad(2,k)+(pt(c)%qp(k)-pt(i)%qp(k))*wy
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

do k=1,nop

   do i=1,pt(k)%nv2f
      c=pt(k)%v2f(i)
      p1=fc(c)%pt(1)
      p2=fc(c)%pt(2)
      ds=fc(c)%area
      dx=fc(c)%nx*ds
      dy=fc(c)%ny*ds
      do j=1,nvar
         var=0.5d0*(pt(p1)%qp(j)+pt(p2)%qp(j))
         pt(k)%grad(1,j)=pt(k)%grad(1,j)+var*dx
         pt(k)%grad(2,j)=pt(k)%grad(2,j)+var*dy
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
