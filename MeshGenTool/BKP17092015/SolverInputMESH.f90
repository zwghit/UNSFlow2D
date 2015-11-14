module grid
implicit none

integer(kind=4):: np,nf,nc 

type points 
     real(kind=8) :: x,y,z
     integer(kind=4):: bc,flag,nv2c 
     integer(kind=4),dimension(:),pointer::v2c
end type points 

type faces  
     integer(kind=4):: pt(2) 
     integer(kind=4):: bc,flag,in,out 
     real(kind=8)   :: sx,sy 
end type faces  

type cells  
     integer(kind=4):: nc2v,nc2f,nc2c  
     real(kind=8)   :: xc,yc,zc,cv
     integer(kind=4),dimension(:),pointer::c2v
     integer(kind=4),dimension(:),pointer::c2f
     integer(kind=4),dimension(:),pointer::c2c
end type cells  

type(points),dimension(:),pointer::pt
type(points),dimension(:,:),pointer::mesh
type(faces),dimension(:),pointer::fc
type(cells),dimension(:),pointer::cell

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
call connectivity

end program SolverMesh


!=====================================================
subroutine solver_compatible
!=====================================================
use grid
implicit none
integer(kind=4) :: i,j,in,out,count(nc),c1,c2
integer(kind=4) :: n1,n2
real(kind=8)    :: dx,dy,ds,x1,y1,x2,y2,nx,ny
real(kind=8)    :: rx,ry,dotprod,xc,yc,dxc,dyc

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

!  faces forming  a cell
cell(:)%nc2f = 0
do i=1,nf
   in = fc(i)%in
   out = fc(i)%out
   if(in/=0) then
     cell(in)%nc2f = cell(in)%nc2f + 1
   endif
   if(out/=0) then
     cell(out)%nc2f = cell(out)%nc2f + 1
   endif
enddo

do i=1,nc
   allocate(cell(i)%c2f(cell(i)%nc2f))
enddo

cell(:)%nc2f = 0
do i=1,nf
   in = fc(i)%in
   out = fc(i)%out
   if(in/=0) then
     cell(in)%nc2f = cell(in)%nc2f + 1
     cell(in)%c2f(cell(in)%nc2f) = i
   endif
   if(out/=0) then
     cell(out)%nc2f = cell(out)%nc2f + 1
     cell(out)%c2f(cell(out)%nc2f) = i
   endif
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
      else
         fc(i)%in = c2
         fc(i)%out = c1
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
   dx = 0.5d0*(x2+x1) ; dy = y2-y1

   if(in/=0) then
    cell(in)%cv = cell(in)%cv + dx*dy 
   endif
   if(out/=0) then
    cell(out)%cv = cell(out)%cv - dx*dy 
   endif

   dx = x2-x1 ; dy = y2-y1
   ds = dsqrt(dx*dx+dy*dy)
   nx = dy/ds ; ny = -dx/ds
   !fc(i)%sx = nx
   !fc(i)%sy = ny
   fc(i)%sx = dy
   fc(i)%sy = -dx
enddo

ds=0.0
do i=1,nc
   ds=ds+cell(i)%cv
enddo
print*,'Domain volume = ',ds

open(3,file='Vec_plot.dat')
open(4,file='Vec_InCell.dat')
open(5,file='Vec_OutCell.dat')
do i=1,nf
   n1=fc(i)%pt(1)
   n2=fc(i)%pt(2)
   x1 = pt(n1)%x ; y1 = pt(n1)%y
   x2 = pt(n2)%x ; y2 = pt(n2)%y
   write(3,100)x1,y1,x2-x1,y2-y1
   !write(21,*)pt(n2)%x,pt(n2)%y
   !write(21,*)
   x1 = (x2+x1)/2.0d0 ; y1= (y2+y1)/2.0d0
   in=fc(i)%in
   out=fc(i)%out
   if(in>0) write(4,100)x1,y1,cell(in)%xc-x1,cell(in)%yc-y1
   if(out>0) write(5,100)x1,y1,cell(out)%xc-x1,cell(out)%yc-y1
enddo
close(3)
close(4)
close(5)

100 format(1x,4(f15.6,1x))

open(3,file='CC_plot.dat')
do i=1,nc
   write(3,*)cell(i)%xc,cell(i)%yc
   write(3,*)
enddo
close(3)


open(13,file='geometry.inp')
write(13,*)np,nc,nf
!write(13,*)'# node'
do i=1,np
 write(13,*)i,pt(i)%x,pt(i)%y,pt(i)%bc
enddo
!write(13,*)'# cell'
do i=1,nc
 write(13,*)i,cell(i)%xc,cell(i)%yc,cell(i)%cv
enddo
!write(13,*)'# face'
do i=1,nf
 write(13,*)i,fc(i)%pt(1),fc(i)%pt(2),fc(i)%in,fc(i)%out,&
                       fc(i)%sx,fc(i)%sy,fc(i)%bc
enddo
close(13)

end subroutine solver_compatible

!=====================================================
subroutine connectivity
!=====================================================
use grid
implicit none
integer(kind=4) :: i,j,k,in,out,count(nc),c1,c2
integer(kind=4) :: n1,n2,m1,m2,nop,noc,nof
integer(kind=4) :: p1,p2,e1,e2,c,f1,f2
real(kind=8)    :: dx,dy,dxc,dyc,ds,x1,y1,x2,y2,nx,ny
real(kind=8)    :: rx,ry,xc,yc
integer(kind=4) :: indx,bubble
real(kind=8)    :: temp,pi,check

type(points),allocatable,dimension(:)::node
type(faces),allocatable,dimension(:)::face
type(cells),allocatable,dimension(:)::elem


open(13,file='geometry.inp')
read(13,*)nop,noc,nof

allocate(node(nop),face(nof),elem(noc))

print*,nop,noc,nof
do i=1,nop
  read(13,*)j,node(i)%x,node(i)%y,node(i)%bc
enddo
do i=1,noc
 read(13,*)j,elem(i)%xc,elem(i)%yc,elem(i)%cv
enddo
do i=1,nof
   read(13,*)j,face(i)%pt(1),face(i)%pt(2),face(i)%in,face(i)%out,&
                       face(i)%sx,face(i)%sy,face(i)%bc
enddo
close(13)


! cells surrounding a node
print*
print*,"==> Cells surrounding a node"

node(:)%nv2c=0
do i=1,nof
   p1=face(i)%pt(1)
   p2=face(i)%pt(2)
   e1=face(i)%in
   e2=face(i)%out

   if(e1/=0) then
         c=0
         do j=1,node(p1)%nv2c
            if(node(p1)%v2c(j)==e1) c=c+1  
         enddo
         if(c==0) then 
           node(p1)%nv2c=node(p1)%nv2c+1
           call alloc_int_ptr(node(p1)%v2c,node(p1)%nv2c)      
           node(p1)%v2c(node(p1)%nv2c)=e1      
         endif
      
         c=0
         do j=1,node(p2)%nv2c
            if(node(p2)%v2c(j)==e1) c=c+1  
         enddo
         if(c==0) then 
           node(p2)%nv2c=node(p2)%nv2c+1
           call alloc_int_ptr(node(p2)%v2c,node(p2)%nv2c)      
           node(p2)%v2c(node(p2)%nv2c)=e1      
         endif
   endif

   if(e2/=0) then
         c=0
         do j=1,node(p1)%nv2c
            if(node(p1)%v2c(j)==e2) c=c+1  
         enddo
         if(c==0) then 
           node(p1)%nv2c=node(p1)%nv2c+1
           call alloc_int_ptr(node(p1)%v2c,node(p1)%nv2c)      
           node(p1)%v2c(node(p1)%nv2c)=e2      
         endif
      
         c=0
         do j=1,node(p2)%nv2c
            if(node(p2)%v2c(j)==e2) c=c+1  
         enddo
         if(c==0) then 
           node(p2)%nv2c=node(p2)%nv2c+1
           call alloc_int_ptr(node(p2)%v2c,node(p2)%nv2c)      
           node(p2)%v2c(node(p2)%nv2c)=e2      
         endif
   endif

enddo   


!do i=1,nop
!   print*,i,(node(i)%v2c(j), j=1,node(i)%nv2c) 
!enddo
print*, 'Max. nbr=',maxval(node(:)%nv2c),'at node',maxloc(node(:)%nv2c)
print*, 'Min. nbr=',minval(node(:)%nv2c),'at node',minloc(node(:)%nv2c)


!Cells to vertex connectivity 
print*
print*,"==> Faces surrounding a cell "

elem(:)%nc2f=0
do i=1,nof
   e1=face(i)%in
   e2=face(i)%out

   if(e1/=0) then
   elem(e1)%nc2f=elem(e1)%nc2f+1
   call alloc_int_ptr(elem(e1)%c2f,elem(e1)%nc2f)      
   elem(e1)%c2f(elem(e1)%nc2f)=i      
   endif

   if(e2/=0) then
   elem(e2)%nc2f=elem(e2)%nc2f+1
   call alloc_int_ptr(elem(e2)%c2f,elem(e2)%nc2f)      
   elem(e2)%c2f(elem(e2)%nc2f)=i      
   endif
enddo


!Cells to vertex connectivity 
print*
print*,"==> Vertex surrounding a cell "

!elem(:)%nc2v=0
!do i=1,nop
!   do j=1,node(i)%nv2c
!      e1=node(i)%v2c(j)
!      elem(e1)%nc2v=elem(e1)%nc2v+1
!      call alloc_int_ptr(elem(e1)%c2v,elem(e1)%nc2v)      
!      elem(e1)%c2v(elem(e1)%nc2v)=i      
!   enddo    
!enddo

elem(:)%nc2v=0
do i=1,noc
   xc=elem(i)%xc
   yc=elem(i)%yc
      n1=face(elem(i)%c2f(1))%pt(1)
      n2=face(elem(i)%c2f(1))%pt(2)
      check=dotprod(n1,n2,xc,yc)
      !print*,'check',check

      if(check<0.0d0) then
         elem(i)%nc2v=elem(i)%nc2v+1
         call alloc_int_ptr(elem(i)%c2v,elem(i)%nc2v)      
         elem(i)%c2v(elem(i)%nc2v)=n1      

         elem(i)%nc2v=elem(i)%nc2v+1
         call alloc_int_ptr(elem(i)%c2v,elem(i)%nc2v)      
         elem(i)%c2v(elem(i)%nc2v)=n2      
      else
         elem(i)%nc2v=elem(i)%nc2v+1
         call alloc_int_ptr(elem(i)%c2v,elem(i)%nc2v)      
         elem(i)%c2v(elem(i)%nc2v)=n2      

         elem(i)%nc2v=elem(i)%nc2v+1
         call alloc_int_ptr(elem(i)%c2v,elem(i)%nc2v)      
         elem(i)%c2v(elem(i)%nc2v)=n1      
      endif
      
   do k=2,elem(i)%nc2f-1
      n2=elem(i)%c2v(elem(i)%nc2v)
      !print*,(elem(i)%c2v(j), j=1,elem(i)%nc2v) 

      do j=k,elem(i)%nc2f
         m1=face(elem(i)%c2f(j))%pt(1)
         m2=face(elem(i)%c2f(j))%pt(2)
         if(m1==n2.or.m2==n2) then 
            check=dotprod(m1,m2,xc,yc)
            !print*,i,k,j,c,check
            !print*,n1,n2,m1,m2 
            if(check<0.0d0) then
               elem(i)%nc2v=elem(i)%nc2v+1
               call alloc_int_ptr(elem(i)%c2v,elem(i)%nc2v)      
               elem(i)%c2v(elem(i)%nc2v)=m2      
            else
               elem(i)%nc2v=elem(i)%nc2v+1
               call alloc_int_ptr(elem(i)%c2v,elem(i)%nc2v)      
               elem(i)%c2v(elem(i)%nc2v)=m1      
            endif
              
            f1=elem(i)%c2f(j)      
            elem(i)%c2f(j)=elem(i)%c2f(k)      
            elem(i)%c2f(k)=f1      
            !n2=elem(i)%c2v(c)
            cycle 
           ! exit
         endif
      enddo 
   enddo 

!   print*,i,(elem(i)%c2v(j), j=1,elem(i)%nc2v) 
enddo

do i=1,noc,2
   do j=1,elem(i)%nc2v-1
      p1=elem(i)%c2v(j)
      p2=elem(i)%c2v(j+1)
      x1 = node(p1)%x ; y1 = node(p1)%y
      x2 = node(p2)%x ; y2 = node(p2)%y
      write(25,101)x1,y1,x2-x1,y2-y1
      !write(25,*)node(p1)%x,node(p1)%y
      !write(25,*)node(p2)%x,node(p2)%y
      write(25,*)
   enddo
   !print*,i,(elem(i)%c2v(j), j=1,elem(i)%nc2v) 
   !exit
enddo

101 format(1x,4(f15.6,1x))
print*, 'Max. cell vertex=',maxval(elem(:)%nc2v),'at node',maxloc(elem(:)%nc2v)
print*, 'Min. cell vertex=',minval(elem(:)%nc2v),'at node',minloc(elem(:)%nc2v)


! cells surrounding a cell
print*
print*,"==> Cells surrounding a cell "

elem(:)%nc2c=0
do i=1,nof
   e1=face(i)%in
   e2=face(i)%out
 
   if(e1/=0.and.e2/=0) then

     c=0
     do j=1,elem(e1)%nc2c
        if(elem(e1)%c2c(j)==e2) c=c+1  
     enddo

     if(c==0) then 
       elem(e1)%nc2c=elem(e1)%nc2c+1
       call alloc_int_ptr(elem(e1)%c2c,elem(e1)%nc2c)      
       elem(e1)%c2c(elem(e1)%nc2c)=e2      
     endif


     c=0
     do j=1,elem(e2)%nc2c
        if(elem(e2)%c2c(j)==e1) c=c+1  
     enddo

     if(c==0) then 
       elem(e2)%nc2c=elem(e2)%nc2c+1
       call alloc_int_ptr(elem(e2)%c2c,elem(e2)%nc2c)      
       elem(e2)%c2c(elem(e2)%nc2c)=e1      
     endif


   endif
   
enddo


!do i=1,noc
!   print*,i,(elem(i)%c2c(j), j=1,elem(i)%nc2c) 
!enddo

print*, 'Max. cell nbr=',maxval(elem(:)%nc2c),'at node',maxloc(elem(:)%nc2c)
print*, 'Min. cell nbr=',minval(elem(:)%nc2c),'at node',minloc(elem(:)%nc2c)

print*
print*

!===========================================================================
!                      Allocate/extend array
!===========================================================================

contains

subroutine alloc_int_ptr(x,n)
use grid
implicit none
integer(kind=4),intent(in) ::n 
integer(kind=4)::i 
integer(kind=4),dimension(:),pointer::temp
integer(kind=4),dimension(:),pointer::x

if (n <= 0) then
 write(*,*) "my_alloc_int_ptr received non-positive dimension. Stop."
 stop
endif

! If not allocated, allocate and return
if (.not.(associated(x))) then
 allocate(x(n))
 return
endif

! If reallocation, create a pointer with a target of new dimension.
allocate(temp(n))

! (1) Expand the array dimension
if ( n > size(x) ) then
   do i = 1, size(x)
      temp(i) = x(i)
   end do

! (2) Shrink the array dimension: the extra data, x(n+1:size(x)), discarded.
else
   do i = 1, n
      temp(i) = x(i)
   end do
endif

! Destroy the target of x
!  deallocate(x)

! Re-assign the pointer
 x => temp


return

end subroutine alloc_int_ptr

function dotprod(n1,n2,xc,yc)
   implicit none
   integer*4:: i,j,n1,n2
   real*8   :: x1,y1,x2,y2
   real*8   :: dx,dy,ds,xc,yc
   real*8   :: dotprod,nx,ny,rx,ry


   x1 = node(n1)%x ; y1 = node(n1)%y
   x2 = node(n2)%x ; y2 = node(n2)%y
   dx = x2-x1 ; dy = y2-y1
   ds = dsqrt(dx*dx+dy*dy)
   nx = dy/ds ; ny = -dx/ds

   dx = xc-x1 ; dy = yc-y1
   ds = dsqrt(dx*dx+dy*dy)
   rx = dx/ds
   ry = dy/ds
   dotprod = rx*nx+ry*ny

end function dotprod


function crossprod(n1,n2,m1,m2)
   implicit none
   integer*4:: i,j,n1,n2,m1,m2
   real*8   :: x1,y1,x2,y2
   real*8   :: x3,y3,x4,y4
   real*8   :: dx,dy,ds,xc,yc
   real*8   :: crossprod,nx,ny,rx,ry


   x1 = node(n1)%x ; y1 = node(n1)%y
   x2 = node(n2)%x ; y2 = node(n2)%y

   x3 = node(m1)%x ; y3 = node(m1)%y
   x4 = node(m2)%x ; y4 = node(m2)%y

   dx = x2-x1 ; dy = y2-y1
   ds = dsqrt(dx*dx+dy*dy)
   nx = dy/ds ; ny = -dx/ds

   dx = x4-x3 ; dy = y4-y3
   ds = dsqrt(dx*dx+dy*dy)
   rx = dx/ds
   ry = dy/ds

   crossprod = nx*ry-ny*rx

end function crossprod



end subroutine connectivity


!subroutine  renumber

