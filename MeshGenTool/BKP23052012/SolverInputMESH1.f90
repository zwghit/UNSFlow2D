module grid
implicit none

integer(kind=4):: np,nf,nc 

type points 
     real(kind=8) :: x,y,z
     integer(kind=4):: bc,flag 
end type points 

type faces  
     integer(kind=4):: pt(2) 
     integer(kind=4):: bc,flag,in,out 
     real(kind=8)   :: sx,sy 
end type faces  

type cells  
     !integer(kind=4):: npt,nfc  
     real(kind=8)   :: xc,yc,zc,cv
     !type(faces),allocatable,dimension(:)::fc 
end type cells  

type(points),allocatable,dimension(:)::pt
type(points),allocatable,dimension(:,:)::mesh
type(faces),allocatable,dimension(:)::fc
type(cells),allocatable,dimension(:)::cell

end module grid
!========================================================
program SolverMesh
use grid
implicit none

integer(kind=4) :: i,j,k,nx,ny,nz,nbk,imax,jmax,kmax,cc,c1
real(kind=8) :: dummy,x0,y0,xi,yi,xj,yj,dx,dy

open(3,file='grid.grd')
read(3,*)nbk
if(nbk>1) stop 'Only Single block is allowed'
read(3,*) nx, ny, nz
print*
print*,'nx,ny,nz :',nx, ny, nz
print*
allocate(mesh(nx,ny))
read(3,*) (((mesh(i,j)%x, i=1,nx), j=1,ny), k=1,nz), &
          (((mesh(i,j)%y, i=1,nx), j=1,ny), k=1,nz), &
          (((dummy      , i=1,nx), j=1,ny), k=1,nz)
close(3)

!For O-grid
!nx=nx-1

np=(nx-1)*ny
nf=(nx-1)*(2*ny-1)
nc=(nx-1)*(ny-1)

print*,'No. of vertices =',np
print*,'No. of faces    =',nf
print*,'No. of cells    =',nc

allocate(pt(np),fc(nf),cell(nc))
cc=0
do j=1,ny
   do i=1,nx-1
      cc=cc+1   
      pt(cc)%x=mesh(i,j)%x      
      pt(cc)%y=mesh(i,j)%y      
      pt(cc)%bc=0
      if(j==1) pt(cc)%bc=1001
      if(j==ny) pt(cc)%bc=2001
      pt(cc)%flag=0
   enddo
enddo
print*
print*,'Node count =',cc

cc=0
do j=1,ny
   do i=1,nx-1
      cc=cc+1   
      fc(cc)%pt(1)=cc
      fc(cc)%pt(2)=cc+1
      if(i==nx-1) then 
        fc(cc)%pt(2)=cc-(nx-2)
      endif 
      fc(cc)%bc=0
      fc(cc)%out=cc-nx+1
      if(j==1) then 
        fc(cc)%bc=1001
        fc(cc)%out=0
      endif
      fc(cc)%in=cc
      if(j==ny) then 
        fc(cc)%bc=2001
        fc(cc)%in=0
      endif
      fc(cc)%flag=0
   enddo
enddo

print*,'Zhi-face count=',cc

c1=0
do j=1,ny-1
   do i=1,nx-1
      cc=cc+1   
      c1=c1+1   
      fc(cc)%pt(1)=c1
      fc(cc)%pt(2)=c1+nx-1
      fc(cc)%bc=0
      fc(cc)%out=c1
      fc(cc)%in=c1-1
      if(i==1) then 
        fc(cc)%bc=2001
        fc(cc)%in=(nx-1)*j
      endif 
      fc(cc)%flag=0
   enddo
enddo

print*
print*,'face count =',cc

call solver_compatible


end program SolverMesh


!=====================================================
subroutine solver_compatible
!=====================================================
use grid
implicit none
integer(kind=4) :: i,j,in,out,count(nc),c1,c2
integer(kind=4) :: n1,n2
real(kind=8)    :: dx,dy,ds,x1,y1,x2,y2,nx,ny
real(kind=8)    :: rx,ry,dotprod

do i=1,nc
   count(i) = 0
   cell(i)%xc = 0.0d0
   cell(i)%yc = 0.0d0
   cell(i)%cv = 0.0d0
enddo

do j=1,nf
   n1  = fc(j)%pt(1)
   n2  = fc(j)%pt(2)
   in  = fc(j)%in
   out = fc(j)%out
   !print*, j,in,out
   x1 = pt(n1)%x ; y1 = pt(n1)%y
   x2 = pt(n2)%x ; y2 = pt(n2)%y

   if(in/=0) then
    cell(in)%xc = cell(in)%xc+x1
    cell(in)%yc = cell(in)%yc+y1
    count(in) = count(in) + 1
    cell(in)%xc = cell(in)%xc+x2
    cell(in)%yc = cell(in)%yc+y2
    count(in) = count(in) + 1
   endif
   if(out/=0) then
    cell(out)%xc = cell(out)%xc+x1
    cell(out)%yc = cell(out)%yc+y1
    count(out) = count(out) + 1
    cell(out)%xc = cell(out)%xc+x2
    cell(out)%yc = cell(out)%yc+y2
    count(out) = count(out) + 1
   endif
enddo

do i = 1,nc
   cell(i)%xc=cell(i)%xc/count(i)
   cell(i)%yc=cell(i)%yc/count(i)
enddo

do i=1,nf
   n1 = fc(i)%pt(1)
   n2 = fc(i)%pt(2)
   c1 = fc(i)%in
   c2 = fc(i)%out
 
   if(n1>np.or.n2>np.or.n1<1.or.n2<1) then
      print*,'fc,n1,n2,in,out:',i,n1,n2,c1,c2
      stop  
   endif    

   x1 = pt(n1)%x ; y1 = pt(n1)%y
   x2 = pt(n2)%x ; y2 = pt(n2)%y
   dx = x2-x1 ; dy = y2-y1
   ds = dsqrt(dx*dx+dy*dy)
   nx = dy/ds ; ny = -dx/ds
   x1 = (x2+x1)/2.0d0 ; y1= (y2+y1)/2.0d0
   if( c1/=0 ) then
      x2 = cell(c1)%xc ; y2 = cell(c1)%yc
      dx = x2-x1 ; dy = y2-y1
      ds = dsqrt(dx*dx+dy*dy)
      rx = dx/ds
      ry = dy/ds
      dotprod = rx*nx+ry*ny
      if(dotprod<0.0d0) then
         fc(i)%in = c1
         fc(i)%out = c2
         print*,'Face Cor :',i
      else
         fc(i)%in = c2
         fc(i)%out = c1
         print*,'Face swap:',i
      endif
   elseif( c1==0 ) then
      x2 = cell(c2)%xc ; y2 = cell(c2)%yc
      dx = x2-x1 ; dy = y2-y1
      ds = dsqrt(dx*dx+dy*dy)
      rx = dx/ds
      ry = dy/ds
      dotprod = rx*nx+ry*ny
      if(dotprod<0.0d0) then
       fc(i)%in = c2
       fc(i)%out = c1
      else
       fc(i)%in = c1
       fc(i)%out = c2
      endif

   endif
enddo





do i=1,nf
   n1 = fc(i)%pt(1)
   n2 = fc(i)%pt(2)
   in = fc(i)%in
   out = fc(i)%out
   x1 = pt(n1)%x ; y1 = pt(n1)%y
   x2 = pt(n2)%x ; y2 = pt(n2)%y
   dx = x2-x1 ; dy = y2-y1
   ds = dsqrt(dx*dx+dy*dy)
   nx = dy/ds ; ny = -dx/ds
   fc(i)%sx = nx*ds
   fc(i)%sy = ny*ds

   if(in/=0) then
    cell(in)%cv = cell(in)%cv + (x1*y2 - x2*y1) + &
                  (x2*cell(in)%yc -cell(in)%xc*y2) + &
                  (cell(in)%xc*y1 -x1*cell(in)%yc)
    cell(in)%cv = 0.5d0* dabs(cell(in)%cv)
   endif

   if(out/=0) then
     cell(out)%cv = cell(out)%cv + (x1*y2 - x2*y1) + &
                    (x2*cell(out)%yc -cell(out)%xc*y2) + &
                    (cell(out)%xc*y1 -x1*cell(out)%yc)
     cell(out)%cv = 0.5d0* dabs(cell(out)%cv)
   endif
enddo

ds=0.0
do i=1,nc
   ds=ds+cell(i)%cv
enddo
print*,'Doamin volume = ',ds


do i=1,nf
   n1=fc(i)%pt(1)
   n2=fc(i)%pt(2)
   x1 = pt(n1)%x ; y1 = pt(n1)%y
   x2 = pt(n2)%x ; y2 = pt(n2)%y
   write(21,100)x1,y1,x2-x1,y2-y1
   !write(21,*)pt(n2)%x,pt(n2)%y
   !write(21,*)
   x1 = (x2+x1)/2.0d0 ; y1= (y2+y1)/2.0d0
   in=fc(i)%in
   out=fc(i)%out
   if(in>0) write(23,100)x1,y1,cell(in)%xc-x1,cell(in)%yc-y1
   if(out>0) write(24,100)x1,y1,cell(out)%xc-x1,cell(out)%yc-y1
enddo

100 format(1x,4(f15.6,1x))

do i=1,nc
   write(22,*)cell(i)%xc,cell(i)%yc
   write(22,*)
enddo


open(13,file='geometry.inp')
write(13,*)np,nc,nf
do i=1,np
 write(13,*)pt(i)%x,pt(i)%y,pt(i)%bc
enddo
do i=1,nc
 write(13,*)cell(i)%xc,cell(i)%yc,cell(i)%cv
enddo
do i=1,nf
 write(13,*)fc(i)%pt(1),fc(i)%pt(2),fc(i)%in,fc(i)%out,&
                       fc(i)%sx,fc(i)%sy,fc(i)%bc
enddo
close(13)
end subroutine solver_compatible

