!------------------------------------------------------------------------------
! Computes finite volume residual
!------------------------------------------------------------------------------
subroutine fvresidual
use grid 
implicit none
integer(kind=i4)  :: i, j, ie, c1, c2

do i=1,noc
   do j=1,nvar
      cell(i)%res(j) = 0.0d0
   enddo
enddo

!     Compute flux for interior edges
do ie=1,nof
   c1 = fc(ie)%in
   c2 = fc(ie)%out
   if(c1 /=0 .and.c2 /=0) call vanleer_flux(ie,c1,c2)
   !if(c1 /=0 .and.c2 /=0) call roe_flux(ie,c1,c2)
enddo


!     Compute flux for solid wall edges

do ie=1,nof
   c1 = fc(ie)%in
   c2 = fc(ie)%out
   if(fc(ie)%bc==1001.and.c1/=0)  call solid_flux(ie,c1)
   if(fc(ie)%bc==1001.and.c2/=0)  call solid_flux(ie,c2)
enddo

!     Flux for far-field points

do ie=1,nof
   c1 = fc(ie)%in
   c2 = fc(ie)%out
   if(fc(ie)%bc==2001.and.c1/=0)  call farfield_flux(ie,c1)
   if(fc(ie)%bc==2001.and.c2/=0)  call farfield_flux(ie,c2)
enddo

end
