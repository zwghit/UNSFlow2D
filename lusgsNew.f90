!-----------------------------------------------------------------------------
! LUSGS implicit scheme - matrix free
!-----------------------------------------------------------------------------
subroutine lusgs
use commons
use pri
use grid
implicit none
!!------------------------------------------------------------------------------
integer(kind=i4):: j,in,out,c,itr,inner_itr
integer  :: i, it, iv
real(dp) :: lam, sx,sy, &
            ds,D(noc), dqc(nvar,noc), &
           cres(nvar), flux1(nvar), flux2(nvar), &
           lam1, omega
real(dp) :: qa(nvar)

!type fvms
!       real(kind=dp),dimension(nvar):: dqc
!       real(kind=dp),dimension(nvar):: qflux,qold
!end type fvms
!type(fvms),allocatable,dimension(:)::fvm
! Over-relaxation factor: higher value improves stability but retards
! convergence. Needs to be tuned.
omega = 1.5d0

! Compute residual vector
call fvresidual



! Compute diagonal term
do it=1,noc
   D(it) = cell(it)%cv/cell(it)%dt
   dqc(:,it)=0.d0
enddo

do i=1,nof

   in = fc(i)%in 
   out= fc(i)%out

   if(in/=0.and.out/=0) then
       call con2prim(cell(in)%qc(:))
       flux1(:)=prim(:)
       call con2prim(cell(out)%qc(:))
       flux2(:)=prim(:)
       
       do iv=1,nvar
          qa(iv) = 0.5d0*(flux1(iv)+flux2(iv))
       enddo
       rho = qa(1)
       u   = qa(2)
       v   = qa(3)
       p   = qa(4)
       a   = dsqrt(GAMMA*p/rho)
    
       sx = fc(i)%sx 
       sy = fc(i)%sy 
       ds =  dsqrt(sx*sx + sy*sy)
       lam =  dabs(u*sx + v*sy) + a*ds
       lam = 0.5d0*omega*lam
       D(in) = D(in) + lam
       D(out) = D(out) + lam
   elseif(in/=0.and.out==0) then
       call con2prim(cell(in)%qc(:))
       sx = fc(i)%sx 
       sy = fc(i)%sy 
       ds =  dsqrt(sx*sx + sy*sy)
       lam =  dabs(u*sx + v*sy) + a*ds
       lam = 0.5d0*omega*lam
       D(in) = D(in) + lam
   elseif(in==0.and.out/=0) then
       call con2prim(cell(out)%qc(:))
       sx = fc(i)%sx 
       sy = fc(i)%sy 
       ds =  dsqrt(sx*sx + sy*sy)
       lam =  dabs(u*sx + v*sy) + a*ds
       lam = 0.5d0*omega*lam
       D(out) = D(out) + lam
   endif
enddo

! Forward loop
do it=1,noc

   do iv=1,nvar
      cres(iv) = 0.0d0
   enddo
  
   do i=1,cell(it)%nc2f
      c=cell(it)%c2f(i)
      in  = fc(c)%in
      out = fc(c)%out
      sx =fc(c)%sx 
      sy =fc(c)%sy 
      if(out==it) then
         out=in
         sx =-fc(c)%sx 
         sy =-fc(c)%sy 
      endif

      ds =  dsqrt(sx*sx + sy*sy)

      if(out .lt. it .and. out .gt. 0)then
         call normalflux(sx, sy,cell(out)%qold(1:nvar),  flux1)
         call normalflux(sx, sy,cell(out)%qc(1:nvar), flux2)
         call maxeig(sx,sy, ds,cell(out)%qold(1:nvar), lam1)
         lam1 = omega*lam1
         do iv=1,nvar
            cres(iv) = cres(iv) + (flux2(iv) - flux1(iv)) -  &
                       lam1*dqc(iv,out)
         enddo
      endif
   enddo

   do iv=1,nvar
      dqc(iv,it) = ( -cell(it)%res(iv) - 0.5d0*cres(iv) )/D(it)
      cell(it)%qc(iv) = cell(it)%qold(iv) + dqc(iv,it)
   enddo

enddo

! Reverse loop
do it=noc,1,-1

   do iv=1,nvar
      cres(iv) = 0.0d0
   enddo

   do i=1,cell(it)%nc2f
      c=cell(it)%c2f(i)
      in = fc(c)%in
      out = fc(c)%out
      sx = fc(c)%sx 
      sy = fc(c)%sy 
      if(out==it) then
         out=in
         sx = -fc(c)%sx 
         sy = -fc(c)%sy 
      endif
      ds =  dsqrt(sx*sx + sy*sy)

      if(out.gt. it .and.out.gt. 0)then
         call normalflux(sx, sy, cell(out)%qold(1:nvar),  flux1)
         call normalflux(sx, sy, cell(out)%qc(1:nvar), flux2)
         call maxeig(sx,sy, ds, cell(out)%qold(1:nvar), lam1)
         lam1 = omega*lam1
         do iv=1,nvar
            cres(iv) = cres(iv) + (flux2(iv) - flux1(iv)) - &
                       lam1*dqc(iv,out)
         enddo
      endif

   enddo

   do iv=1,nvar
      dqc(iv,it) = dqc(iv,it) - 0.5d0*cres(iv)/D(it)
      cell(it)%qc(iv)  = cell(it)%qold(iv) + dqc(iv,it)
   enddo

enddo

end
!-----------------------------------------------------------------------------
! Computes flux along (sx,sy)
!-----------------------------------------------------------------------------
subroutine normalflux(sx, sy, qc, flux)
use commons
use pri
use grid
implicit none
real(dp) :: sx, sy, qc(nvar), flux(nvar)

real(dp) ::  un,et

call con2prim(qc)

!et       = p/GAMMA1 + 0.5d0*rho*(u*u + v*v)
et=e*rho
un      = u*sx + v*sy
flux(1) = (et + p)*un
flux(2) = rho*un
flux(3) = p*sx + u*flux(2)
flux(4) = p*sy + v*flux(2)

end
!-----------------------------------------------------------------------------
! Computes maximum eigenvalue normal to a face with normal (sx, sy)
! and ds = dsqrt(sx*sx + sy*sy) is face length
!-----------------------------------------------------------------------------
subroutine maxeig(sx, sy, ds, qc, lam)
use commons
use pri
use grid
implicit none
real(dp) :: sx, sy, ds, qc(nvar), lam

call con2prim(qc)
!a   = dsqrt(GAMMA*p/d)
lam = dabs(u*sx + v*sy) + a*ds

end
