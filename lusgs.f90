!-----------------------------------------------------------------------------
! LUSGS implicit scheme - matrix free
!-----------------------------------------------------------------------------
subroutine lusgs
use commons
use pri
use grid
implicit none
!!------------------------------------------------------------------------------
integer(kind=i4):: i,j,it,iv,in,out,c,itr,inner_itr
real(kind=dp):: omega,sign 
real(kind=dp):: D(noc),cres(nvar), flux1(nvar), flux2(nvar) 
real(kind=dp):: sx,sy,ds
type fvms
       real(kind=dp),dimension(nvar):: dqc
!       real(kind=dp),dimension(nvar):: qflux,qold
end type fvms
type(fvms),allocatable,dimension(:)::fvm



allocate(fvm(noc))
! Over-relaxation factor: higher value improves stability but retards
! convergence. Needs to be tuned.
inner_itr=1
omega = 1.5d0

! Compute residual vector
call fvresidual



! Compute diagonal term
do it=1,noc
   D(it)=cell(it)%cv/cell(it)%dt+0.5d0*omega*cell(it)%la
   fvm(it)%dqc=0.d0
!   fvm(it)%qflux(:)=0.d0
!!   fvm(it)%qold(:)=cell(it)%qold(:)
enddo


do itr=1,inner_itr

! Forward loop
do it=1,noc
   cres(:) = 0.0d0

   do j=1,cell(it)%nc2f

      c=cell(it)%c2f(j) 
      sx = fc(c)%sx
      sy = fc(c)%sy
      ds = dsqrt(sx*sx + sy*sy)
      in = fc(c)%in
      out = fc(c)%out
      sign=1.d0
      if(out==it) then  
         out=in  
         sign=-1.d0 
      endif

      if(out .lt. it .and. out .gt. 0)then
         flux1(:)=cell(out)%qold(:)
         call normalflux(flux1)

         flux2(:)=cell(out)%qc(:)
         call normalflux(flux2)

         do iv=1,nvar
            cres(iv)=cres(iv)+ 0.5d0*( sign*(flux2(iv)-flux1(iv))-omega*fc(c)%la*fvm(out)%dqc(iv))
         enddo
      endif

   enddo

   do iv=1,nvar
      fvm(it)%dqc(iv)=(-cell(it)%res(iv)-cres(iv))/D(it)
      cell(it)%qc(iv) = cell(it)%qold(iv) + fvm(it)%dqc(iv)
   enddo

enddo

! Reverse loop
do it=noc,1,-1
   cres(:) = 0.0d0

   do j=1,cell(it)%nc2f

      c=cell(it)%c2f(j) 
      sx = fc(c)%sx
      sy = fc(c)%sy
      ds = dsqrt(sx*sx + sy*sy)
      in = fc(c)%in
      out = fc(c)%out
      sign=1.d0
      if(out==it) then  
         out=in  
         sign=-1.d0
      endif

      if(out .gt. it .and. out .gt. 0)then
         flux1(:)=cell(out)%qold(:)
         call normalflux(flux1)

         flux2(:)=cell(out)%qc(:)
         call normalflux(flux2)

         do iv=1,nvar
            cres(iv)=cres(iv)+ 0.5d0*( sign*(flux2(iv)-flux1(iv))-omega*fc(c)%la*fvm(out)%dqc(iv))
         enddo
      endif

   enddo

   do iv=1,nvar
      fvm(it)%dqc(iv)=fvm(it)%dqc(iv)-cres(iv)/D(it)
      cell(it)%qc(iv) = cell(it)%qold(iv) + fvm(it)%dqc(iv)
   enddo

enddo

enddo

do it=1,noc
   do iv=1,nvar
      cell(it)%qc(iv) = cell(it)%qold(iv) + fvm(it)%dqc(iv)
   enddo
enddo


deallocate(fvm)

contains

!-----------------------------------------------------------------------------
! Computes flux along (sx,sy)
!-----------------------------------------------------------------------------
subroutine normalflux(flux)
use param
use pri
implicit none
real(kind=dp) :: flux(nvar)
real(kind=dp) :: un,et

call con2prim(flux)
!e       = p/GAMMA1 + 0.5d0*rho*(u*u + v*v)

et=e*rho
un      = u*sx + v*sy
flux(1) = (Et+ p)*un
flux(2) = rho*un
flux(3) = flux(2)*u + p*sx
flux(4) = flux(2)*v + p*sy

end subroutine normalflux

end
