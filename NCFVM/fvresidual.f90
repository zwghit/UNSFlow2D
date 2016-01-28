!------------------------------------------------------------------------------
! Computes finite volume residual
!------------------------------------------------------------------------------
subroutine fvresidual
use grid 
implicit none
integer(kind=i4)  :: i, j, ie, c1, c2
real(kind=dp) :: nx, ny,x,y


!call Gradient_LSQR
call Gradient_GG


do i=1,nop
   do j=1,nvar
      pt(i)%res(j) = 0.0d0
   enddo
enddo
pt(:)%la=0.d0
fc(:)%la=0.d0

!     Compute flux for interior edges
do ie=1,nof
   c1 = fc(ie)%pt(1)
   c2 = fc(ie)%pt(2)
   call vanleer_flux(ie,c1,c2)
   !if(fc(ie)%bc==0) call vanleer_flux(ie,c1,c2)
   !if(c1 /=0 .and.c2 /=0) call vanleer_flux(ie,c2,c1)
   !if(c1 /=0 .and.c2 /=0) call roe_flux(ie,c1,c2)
   !if(c1 /=0 .and.c2 /=0) call ausmPlus_flux(ie,c1,c2)
   if(fc(ie)%bc==1001)  call solid_flux(ie,c1,c2)
   if(fc(ie)%bc==2001)  call farfield_flux1(ie,c1,c2)
enddo

call time_step2
!call time_step02

end
