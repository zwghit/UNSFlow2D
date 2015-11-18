!-----------------------------------------------------------------------------
!.....Definition of some constants
!-----------------------------------------------------------------------------
subroutine math
use param
implicit none


GAMMA        = 1.4d0
GAMMA1       = GAMMA-1.0d0
GAS_CONST    = 1.0d0
PI         = 4.0d0*datan(1.0d0)

end
!-----------------------------------------------------------------------------
!.....Variables stored are primitive - density, u, v, pressure
!.....Initialize primitive variables to free stream values
!-----------------------------------------------------------------------------
subroutine initialize
use param
use grid
implicit none
integer(kind=i4)  :: j

q_inf = 1.d0
aoa     = aoa_deg*PI/180.0d0
qinf(2) = 1.d0
qinf(3) = dcos(aoa)*q_inf
qinf(4) = dsin(aoa)*q_inf
p_inf    = 1.d0/(gamma*m_inf*m_inf)
qinf(1) = p_inf/(gamma-1.d0) + 0.5d0
r_inf = qinf(2)
ent_inf = p_inf/(r_inf**gamma)
a_inf   = dsqrt(gamma*p_inf/r_inf)

u_inf=qinf(3)/qinf(2)
v_inf=qinf(4)/qinf(2)

!     Print some useful info
if(flow_type == 'inviscid') print*,'Euler computation'
if(flow_type == 'laminar')  print*,'Laminar Navier-Stokes computation'
if(flow_type == 'rans')print*,'Turbulent Navier-Stokes computation'
print*,'Free-stream values:'
write(*,'(5x, " Mach number =", f8.4)')m_inf
write(*,'(5x, " AOA         =", f8.4)')aoa_deg
write(*,'(5x, " u velocity  =", f8.4)')u_inf
write(*,'(5x, " v velocity  =", f8.4)')v_inf
write(*,'(5x, " Pressure    =", f15.6)')p_inf

do j=1,noc
   cell(j)%qc(:) = qinf(:)
enddo

end
!-----------------------------------------------------------------------------
! Save flow solution into a file. Conserved variables are saved.
!-----------------------------------------------------------------------------
subroutine save_flow
use grid
implicit none
integer(kind=i4)  :: i, j

open(unit=50, file='FLO.DAT')
do i=1,noc
   write(50,'(4e20.10)') (cell(i)%qc(j), j=1,nvar)
enddo
close(50)

end
!-----------------------------------------------------------------------------
! Read flow solution from file
!-----------------------------------------------------------------------------
subroutine read_flow
use grid
implicit none
integer(kind=i4)  :: i, j

open(unit=50, file='FLO.DAT', status='OLD')
do i=1,noc
   read(50,*)(cell(i)%qc(j), j=1,nvar)
enddo
close(50)

end
!-----------------------------------------------------------------------------
! Save old solution
!-----------------------------------------------------------------------------
subroutine save_old
use param
use grid
use pri
implicit none
integer(kind=i4)  :: i, j

do i=1,noc
do j=1,nvar
      cell(i)%qold(j) = cell(i)%qc(j)
enddo
enddo


end
!-----------------------------------------------------------------------------
! L2 and Linf norm of the finite volume residual
!-----------------------------------------------------------------------------
subroutine residue
use param
use grid
implicit none
integer(kind=i4)  :: i, j
real(kind=dp) :: fr,dt

! Store current residual norm into fres_old for gmres
fres_old = fres

fres  = 0.0d0
fresi = 0.0d0
iresi = 0

!do i=1,noc
!      fresi = cell(i)%qc(2)/cell(i)%qold(2) - 1.d0
!      dt=cell(i)%dt
!      fres = fres + fresi*fresi/(dt*dt)
!enddo
!     fres=dsqrt(fres)/noc
!return

do i=1,noc
   fr = 0.0d0
   do j=1,nvar
      fr = fr + cell(i)%res(j)**2
   enddo
   fres = fres + fr
   fr   = dsqrt(fr)
   if(fr .gt. fresi)then
      fresi = fr
      iresi = i
   endif
enddo

fres = dsqrt(fres/noc)

if(iter .eq. 1)then
   fres1 = fres
   print*,'Residue in first iteration =',fres1
endif

if(fres1 .ne. 0.0d0) fres = fres/fres1



end


!-----------------------------------------------------------------------------
! Time-step from cfl condition: refined version
!-----------------------------------------------------------------------------
!======================================================================================
subroutine time_step2
use param
use pri
use grid
!use commons
implicit none
integer(kind=i4) :: i

call spectral

do i=1,noc
cell(i)%dt=0.d0
enddo

dtglobal = 1.0d20
do i=1,noc
   cell(i)%dt  = cfl*cell(i)%cv/cell(i)%la ! local  timestep
   dtglobal = dmin1(dtglobal, cell(i)%dt)   ! global timestep
   !if(cell(i)%dt<=0.0) print*,i,cell(i)%dt
enddo
end subroutine time_step2

!======================================================================================
!  time  step
subroutine time_step02
use param
use grid
use pri
use commons
implicit none
integer(kind=i4) :: i
real(kind=dp)    :: ll,Ui,Vi,pr,rhoi,ai
!real(kind=dp)    :: con(nvar) 

call spectral

cell(:)%dt=0.d0
dtglobal = 1.0d20
do i=1,noc
   rhoi=cell(i)%qp(1)
   Ui  =cell(i)%qp(2)
   Vi  =cell(i)%qp(3)
   pr  =cell(i)%qp(4)
   ai=dsqrt(gamma*pr/rhoi)
   ll  = (dabs(Ui) + ai)*cell(i)%ds+(dabs(Vi) + ai)*cell(i)%ds
   cell(i)%dt  = cfl*cell(i)%cv/ll          ! local  timestep
   dtglobal = dmin1(dtglobal, cell(i)%dt)   ! global timestep
   if(cell(i)%dt<=0.0) print*,i,cell(i)%dt
enddo



end
!======================================================================================

subroutine spectral 
use param
use pri
use grid
!use commons
implicit none
integer(kind=i4) :: i,j,in,out
real(kind=dp)    :: r0, u0, v0, p0, a0, nx, ny,un,nl,ll
real(kind=dp)    :: con(nvar),prim1(nvar),prim2(nvar)

do i=1,noc
cell(i)%la=0.d0
enddo

do i=1,nof

   in=fc(i)%in
   out=fc(i)%out
   nx  = fc(i)%sx
   ny  = fc(i)%sy
   nl  = dsqrt(nx*nx + ny*ny)
    
   if(in/=0.and.out/=0) then
      do j=1,nvar
         con(j) = 0.5d0*( cell(in)%qp(j) + cell(out)%qp(j) )
      enddo

      r0   = con(1)
      u0   = con(2)
      v0   = con(3)
      p0   = con(4)
      a0   = dsqrt(gamma*p0/r0)
      un  = u0*nx + v0*ny
      ll  = dabs(un) + a0*nl
      cell(in)%la=cell(in)%la+ll
      cell(out)%la=cell(out)%la+ll
      fc(i)%la=ll
   endif
   if(in/=0.and.out==0) then
      con(:)=cell(in)%qc(:)
      call con2prim(con)
      un  = u*nx + v*ny
      ll  = dabs(un) + a*nl
      cell(in)%la=cell(in)%la+ll
      fc(i)%la=ll
   endif

   if(in==0.and.out/=0) then
      con(:)=cell(out)%qc(:)
      call con2prim(con)
      un  = u*nx + v*ny
      ll  = dabs(un) + a*nl
      cell(out)%la=cell(out)%la+ll
      fc(i)%la=ll
   endif

enddo


end subroutine spectral 
!======================================================================================
subroutine Gradient_GG
use grid 
use commons
implicit none
integer(kind=i4) :: i,j,k,in,out,c,p1,p2,ps
real(kind=dp) :: var,dx,dy,xc,yc 
real(kind=dp) :: ql(nvar),qr(nvar)
real(kind=dp) :: nr, dr, phi, q_min, q_max,kappa
real(kind=dp) :: TOL,alfa,D_L,D_L0,check,dotprod

do i=1,noc
   cell(i)%qx(:)=0.d0
   cell(i)%qy(:)=0.d0
enddo

do k=1,noc

   do i=1,cell(k)%nc2v-1
      c=cell(k)%c2v(i)
      p1=cell(k)%c2v(i)
      p2=cell(k)%c2v(i+1)
      dx=pt(p2)%y-pt(p1)%y 
      dy=-(pt(p2)%x-pt(p1)%x)
      do j=1,nvar
         var=0.5d0*(pt(p1)%prim(j)+pt(p2)%prim(j))
         cell(k)%qx(j)=cell(k)%qx(j)+var*dx
         cell(k)%qy(j)=cell(k)%qy(j)+var*dy
      enddo
   enddo

enddo

do i=1,noc
   var=cell(i)%cv
   do j=1,nvar
   cell(i)%qx(j)=cell(i)%qx(j)/var
   cell(i)%qy(j)=cell(i)%qy(j)/var
   enddo
enddo

!return

call limiter_VKN

end subroutine Gradient_GG

!======================================================================================

subroutine Gradient_GG_FC
use grid 
use commons
implicit none
integer(kind=i4) :: i,j,k,ie,in,out,p1,p2,c
real(kind=dp) :: q2, pv(nvar)
real(kind=dp)    :: x1,y1,x2,y2,dx,dy

real(kind=dp) :: ql(nvar),qr(nvar)
real(kind=dp) :: nr, dr, phi, q_min, q_max,kappa
real(kind=dp) :: TOL,alfa,D_L,D_L0


do i=1,noc
   cell(i)%qx(:)=0.d0
   cell(i)%qy(:)=0.d0
enddo

do i=1,nof
      in = fc(i)%in
      out = fc(i)%out
      p1=fc(i)%pt(1)       
      p2=fc(i)%pt(2)       

      do j=1,nvar    
          fc(i)%qx(j)=0.d0
          fc(i)%qy(j)=0.d0
          if(in/=0.and.out/=0) then
                  pv(j)=0.5d0*(pt(p1)%prim(j)+cell(out)%qp(j)) 
                  x1 = pt(p1)%x    ; y1 = pt(p1)%y
                  x2 = cell(out)%xc ; y2 = cell(out)%yc
                  dx = y2-y1       ;  dy = -(x2-x1) 
                  fc(i)%qx(j)=fc(i)%qx(j)+pv(j)*dx
                  fc(i)%qy(j)=fc(i)%qy(j)+pv(j)*dy
         
                  pv(j)=0.5d0*(cell(out)%qp(j)+pt(p2)%prim(j)) 
                  x1 = cell(out)%xc ; y1 = cell(out)%yc
                  x2 = pt(p2)%x    ; y2 = pt(p2)%y
                  dx = y2-y1       ;  dy = -(x2-x1) 
                  fc(i)%qx(j)=fc(i)%qx(j)+pv(j)*dx
                  fc(i)%qy(j)=fc(i)%qy(j)+pv(j)*dy
              
                  pv(j)=0.5d0*(pt(p2)%prim(j)+cell(in)%qp(j)) 
                  x1 = pt(p2)%x    ; y1 = pt(p2)%y
                  x2 = cell(in)%xc ; y2 = cell(in)%yc
                  dx = y2-y1       ;  dy = -(x2-x1) 
                  fc(i)%qx(j)=fc(i)%qx(j)+pv(j)*dx
                  fc(i)%qy(j)=fc(i)%qy(j)+pv(j)*dy
              
                  pv(j)=0.5d0*(cell(in)%qp(j)+pt(p1)%prim(j)) 
                  x1 = cell(in)%xc ; y1 = cell(in)%yc
                  x2 = pt(p1)%x    ; y2 = pt(p1)%y
                  dx = y2-y1       ;  dy = -(x2-x1) 
                  fc(i)%qx(j)=fc(i)%qx(j)+pv(j)*dx
                  fc(i)%qy(j)=fc(i)%qy(j)+pv(j)*dy
          endif
    
          if(in==0.and.out/=0) then
                  pv(j)=0.5d0*(pt(p1)%prim(j)+cell(out)%qp(j)) 
                  x1 = pt(p1)%x    ; y1 = pt(p1)%y
                  x2 = cell(out)%xc ; y2 = cell(out)%yc
                  dx = y2-y1       ;  dy = -(x2-x1) 
                  fc(i)%qx(j)=fc(i)%qx(j)+pv(j)*dx
                  fc(i)%qy(j)=fc(i)%qy(j)+pv(j)*dy
         
                  pv(j)=0.5d0*(cell(out)%qp(j)+pt(p2)%prim(j)) 
                  x1 = cell(out)%xc ; y1 = cell(out)%yc
                  x2 = pt(p2)%x    ; y2 = pt(p2)%y
                  dx = y2-y1       ;  dy = -(x2-x1) 
                  fc(i)%qx(j)=fc(i)%qx(j)+pv(j)*dx
                  fc(i)%qy(j)=fc(i)%qy(j)+pv(j)*dy
         
                  pv(j)=0.5d0*(pt(p2)%prim(j)+pt(p1)%prim(j)) 
                  x1 = pt(p2)%x ; y1 = pt(p2)%y
                  x2 = pt(p1)%x ; y2 = pt(p1)%y
                  dx = y2-y1    ;  dy = -(x2-x1) 
                  fc(i)%qx(j)=fc(i)%qx(j)+pv(j)*dx
                  fc(i)%qy(j)=fc(i)%qy(j)+pv(j)*dy
          endif
    
          if(out==0.and.in/=0) then
                 pv(j)=0.5d0*(pt(p1)%prim(j)+pt(p2)%prim(j)) 
                 x1 = pt(p1)%x ; y1 = pt(p1)%y
                 x2 = pt(p2)%x ; y2 = pt(p2)%y
                 dx = y2-y1    ;  dy = -(x2-x1) 
                 fc(i)%qx(j)=fc(i)%qx(j)+pv(j)*dx
                 fc(i)%qy(j)=fc(i)%qy(j)+pv(j)*dy
        
                 pv(j)=0.5d0*(pt(p2)%prim(j)+cell(in)%qp(j)) 
                 x1 = pt(p2)%x    ; y1 = pt(p2)%y
                 x2 = cell(in)%xc ; y2 = cell(in)%yc
                 dx = y2-y1       ;  dy = -(x2-x1) 
                 fc(i)%qx(j)=fc(i)%qx(j)+pv(j)*dx
                 fc(i)%qy(j)=fc(i)%qy(j)+pv(j)*dy
             
                 pv(j)=0.5d0*(cell(in)%qp(j)+pt(p1)%prim(j)) 
                 x1 = cell(in)%xc ; y1 = cell(in)%yc
                 x2 = pt(p1)%x    ; y2 = pt(p1)%y
                 dx = y2-y1       ;  dy = -(x2-x1) 
                 fc(i)%qx(j)=fc(i)%qx(j)+pv(j)*dx
                 fc(i)%qy(j)=fc(i)%qy(j)+pv(j)*dy
          endif
    
          fc(i)%qx(j)=fc(i)%qx(j)/fc(i)%cov
          fc(i)%qy(j)=fc(i)%qy(j)/fc(i)%cov
    
          if(in/=0) then
          cell(in)%qx(j)=cell(in)%qx(j)+fc(i)%qx(j)*fc(i)%cov  
          cell(in)%qy(j)=cell(in)%qy(j)+fc(i)%qy(j)*fc(i)%cov  
          endif
    
          if(out/=0) then
          cell(out)%qx(j)=cell(out)%qx(j)+fc(i)%qx(j)*fc(i)%cov  
          cell(out)%qy(j)=cell(out)%qy(j)+fc(i)%qy(j)*fc(i)%cov  
          endif
      enddo 

enddo


do i=1,noc
   do j=1,nvar
   cell(i)%qx(j)= cell(i)%qx(j)/cell(i)%cov
   cell(i)%qy(j)= cell(i)%qy(j)/cell(i)%cov
   enddo
enddo

call limiter_VKN


100 format(1x,4(f15.8,1x))
end subroutine Gradient_GG_FC

!======================================================================================
subroutine Gradient_LSQR
use grid 
use commons
implicit none
integer(kind=i4) :: i,j,k,c
real(kind=dp)    :: xc,yc,dx,dy,wt
real(kind=dp)    :: wx,wy     
real(kind=dp)    :: r11,r12,r22,alfa1,alfa2


do i=1,noc
   cell(i)%qx(:)=0.d0
   cell(i)%qy(:)=0.d0
enddo

do i=1,noc
   xc=cell(i)%xc
   yc=cell(i)%yc
   r11=cell(i)%r11
   r12=cell(i)%r12
   r22=cell(i)%r22

  do k=1,nvar
   do j=1,cell(i)%nc2c
      c=cell(i)%c2c(j)
      dx=cell(c)%xc-xc
      dy=cell(c)%yc-yc

      alfa1=dx/r11/r11
      alfa2=(dy-dx*r12/r11)/r22/r22
      wx=alfa1-alfa2*r12/r11
      wy=alfa2

      cell(i)%qx(k)=cell(i)%qx(k)+(cell(c)%qp(k)-cell(i)%qp(k))*wx
      cell(i)%qy(k)=cell(i)%qy(k)+(cell(c)%qp(k)-cell(i)%qp(k))*wy
    enddo
   enddo

enddo

call limiter_VKN

end subroutine Gradient_LSQR

!======================================================================================

subroutine limiter_VKN 
use grid 
use commons
implicit none
integer(kind=i4) :: i,j,k,ie,in,out,p1,p2,c
real(kind=dp) :: q2, pv(nvar)
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
      cell(i)%phi(j)=1e20
   enddo  
enddo

kappa=0.10
do i=1,noc
   
   !TOL=(kappa*cell(in)%ds)**2
   TOL=(kappa*dsqrt(cell(i)%cv))**3
   do j=1,nvar 

      do k=1,cell(i)%nc2c
         c=cell(i)%c2c(k)
         D_L=cell(c)%qp(j)-cell(i)%qp(j)
         !D_L0=dsign(D_L,1.d0)*(dabs(D_L)+eps) 
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
            phi=0.d0
         endif 
      
         cell(i)%phi(j)=dmin1(cell(i)%phi(j),phi)

      enddo  
   enddo  
enddo  

do i=1,noc
   do j=1,nvar 
      cell(i)%qx(j)=cell(i)%qx(j)*cell(i)%phi(j)
      cell(i)%qy(j)=cell(i)%qy(j)*cell(i)%phi(j)
   enddo
enddo

end subroutine limiter_VKN 

!-----------------------------------------------------------------------------

!!.....Error function, from Abromovitz and Stegun
!!-----------------------------------------------------------------------------
!double precision function ERRF(X)
!double precision X,ARG,E,VB,T,tmp1,tmp2,tmp3,tmp4,tmp5,tmp6
!
!ARG = X*X
!if(ARG .lt. 20.0d0)then
!      E = exp(-ARG)
!else
!      E = 0.0d0
!endif
!VB = abs(X)
!T = 1.0d0/(1.0d0 + 0.3275911d0*VB)
!tmp1 = 1.061405429d0*T
!tmp2 = (tmp1 - 1.453152027d0)*T
!tmp3 = (tmp2 + 1.421413741d0)*T
!tmp4 = (tmp3 - 0.284496736d0)*T
!tmp5 = (tmp4 + 0.254829592d0)*T
!tmp6 = 1.0d0 - tmp5*E
!if(X .lt. 0.0d0)then
!      ERRF = -tmp6
!else
!      ERRF =  tmp6
!endif
!end
!!-----------------------------------------------------------------------------
!!.....Prints data into a file for visualizing contours using gnuplot
!!.....Taken from NSC2KE of Bijan Mohammadi
!!-----------------------------------------------------------------------------
SUBROUTINE isocont(ifile,F,COOR,NVAL,VAL)
implicit none
integer          ifile, nval
double precision F(3),COOR(2,3),VAL(100)

integer          IP1(3), ival, itr, k
double precision epsi, ff1, ff2, ff3, ffma, d12, d23, val1, fk, &
                 fk1, fmi, fma, dif, eps, hh, x, y
!
epsi   = 1.0d-5
IP1(1) = 2
IP1(2) = 3
IP1(3) = 1
FF1    = F(1)
FF2    = F(2)
FF3    = F(3)
FFMA   = DMAX1(DABS(FF1),DABS(FF2))
FFMA   = DMAX1(ffma,DABS(FF3))
D12    = DABS(FF1-FF2)
D23    = DABS(FF2-FF3)
IF(D12+D23.LT.DMAX1(epsi,epsi*FFMA)) GOTO 1000
!  PAS DE RESTRICTION
!  ******************    
DO 100 IVAL=1,NVAL
      VAL1 = VAL(IVAL)
      ITR  = 0
      DO 110 K=1,3
            FK  = F(K)
            FK1 = F(IP1(K))
            FMI = DMIN1(FK,FK1)
            FMA = DMAX1(FK,FK1)
            DIF = FMA-FMI
            IF(DIF.LT.epsi) GOTO 110
            EPS = epsi*DIF
            IF(VAL1.LT.FMI-EPS .OR. VAL1.GT.FMA+EPS) GOTO 110
            HH  = DABS(FK-VAL1)/DIF
            X   = COOR(1,K) + HH*(COOR(1,IP1(K))-COOR(1,K))
            Y   = COOR(2,K) + HH*(COOR(2,IP1(K))-COOR(2,K))
            IF(ITR.EQ.0) GOTO 115
            write(ifile,*) x,y
            write(ifile,*) 
            GOTO 100
115               ITR = 1
            write(ifile,*) x,y
110         CONTINUE
100   CONTINUE
1000  return      
END



!==============================================================================
real(kind=dp) function m1p(mach)
use data_type
implicit none

real(kind=dp):: mach

m1p=0.50*(mach+dabs(mach))

end function m1p
!====================================
real(kind=dp) function m1m(mach)
use data_type
implicit none

real(kind=dp):: mach

m1m=0.5d0*(mach-dabs(mach))

end function m1m
!====================================
real(kind=dp) function m2p(mach)
use data_type
implicit none

real(kind=dp):: mach

m2p=0.25d0*(mach+1.d0)*(mach+1.d0)

end function m2p
!====================================
real(kind=dp) function m2m(mach)
use data_type
implicit none

real(kind=dp):: mach

m2m=-0.25d0*(mach-1.d0)*(mach-1.d0)

end function m2m
!====================================
real(kind=dp) function m4p(mach,beta)
use data_type
implicit none

real(kind=dp):: mach,beta,m2m,m2p

if(dabs(mach)>=1.d0) then
   m4p=0.5d0*(mach+dabs(mach))
else
   m4p=m2p(mach)*(1.d0-16.d0*beta*m2m(mach))
endif

end function m4p
!====================================
real(kind=dp) function m4m(mach,beta)
use data_type
implicit none

real(kind=dp):: mach,beta,m2m,m2p

if(dabs(mach)>=1.d0) then
   m4m=0.5d0*(mach-dabs(mach))
else
   m4m=m2m(mach)*(1.d0+16.d0*beta*m2p(mach))
endif

end function m4m
!====================================
real(kind=dp) function p5p(mach,alpha)
use data_type
implicit none

real(kind=dp):: mach,alpha,m2m,m2p

if(dabs(mach)>=1.d0) then
   p5p=0.5d0*(mach+dabs(mach))/mach
else
   p5p=m2p(mach)*( (2.d0-mach)-16.d0*alpha*mach*m2m(mach))
endif

end function p5p
!====================================
real(kind=dp) function p5m(mach,alpha)
use data_type
implicit none

real(kind=dp):: mach,alpha,m2m,m2p

if(dabs(mach)>=1.d0) then
   p5m=0.5d0*(mach-dabs(mach))/mach
else
   p5m=m2m(mach)*( (-2.d0-mach)+16.d0*alpha*mach*m2p(mach))
endif

end function p5m
!==============================================================================


subroutine avg_c2v 
use grid 
use pri
use commons
implicit none
integer(kind=i4) :: i,j,c,k
real(kind=dp) :: con(nvar)
real(kind=dp) :: wt,cwt


do i=1,noc
      call con2prim(cell(i)%qc(1:nvar))
      do k=1,nvar
         cell(i)%qp(k)=prim(k)
      enddo
enddo

do i=1,nop

    wt=0.d0
    con(:)=0.d0
    do j=1,pt(i)%nv2c
       c=pt(i)%v2c(j)
       !cwt=1.d0/cell(c)%cv
       cwt=pt(i)%wt(j)
       wt=wt+cwt
       do k=1,nvar
          con(k)=con(k)+cell(c)%qp(k)*cwt 
       enddo
    enddo
    do k=1,nvar
       pt(i)%prim(k)=con(k)/wt    
    enddo
end do


end subroutine avg_c2v

! Gradients
subroutine avg_Grad_c2v 
use grid 
use pri
use commons
implicit none
integer(kind=i4) :: i,j,c
real(kind=dp) :: con(nvar),qc(nvar)
real(kind=dp) :: q2, mach,entropy,mul,mutot,sutherland
real(kind=dp) :: r0,u0,v0,p0,t0,wt,x1,x2,y1,y2,cwt
real(kind=dp) :: r1,u1,v1,p1,t1



do i=1,nop

    u0=0.d0
    v0=0.d0
    p0=0.d0
    r0=0.d0
    u1=0.d0
    v1=0.d0
    p1=0.d0
    r1=0.d0
    wt=0.d0
    do j=1,pt(i)%nv2c
       c=pt(i)%v2c(j)
       !cwt=1.d0/cell(c)%cv
       cwt=pt(i)%wt(j)
       wt=wt+cwt
       r0    = r0 + cell(c)%qx(1)*cwt
       u0    = r0 + cell(c)%qx(2)*cwt
       v0    = r0 + cell(c)%qx(3)*cwt
       p0    = r0 + cell(c)%qx(4)*cwt

       r1    = r1 + cell(c)%qy(1)*cwt
       u1    = u1 + cell(c)%qy(2)*cwt
       v1    = v1 + cell(c)%qy(3)*cwt
       p1    = p1 + cell(c)%qy(4)*cwt
    enddo
    pt(i)%dx(1)=r0/wt    
    pt(i)%dx(2)=u0/wt    
    pt(i)%dx(3)=v0/wt    
    pt(i)%dx(4)=p0/wt    

    pt(i)%dy(1)=r1/wt    
    pt(i)%dy(2)=u1/wt    
    pt(i)%dy(3)=v1/wt    
    pt(i)%dy(4)=p1/wt    

end do

end subroutine avg_Grad_c2v 
