module grid
implicit none

integer(kind=4):: np,nf,nc 

type points 
     real(kind=8) :: x,y,z
     integer(kind=4):: bc,flag,cn 
     integer(kind=4),dimension(:),pointer::cnbr
end type points 

type faces  
     integer(kind=4):: pt(2) 
     integer(kind=4):: bc,flag,in,out 
     real(kind=8)   :: sx,sy 
end type faces  

type cells  
     integer(kind=4):: nfc,cn,vtx  
     real(kind=8)   :: xc,yc,zc,cv
     integer(kind=4),dimension(:),pointer::fc 
     integer(kind=4),dimension(:),pointer::cnbr
     integer(kind=4),dimension(:),pointer::c2v
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
cell(:)%nfc = 0
do i=1,nf
   in = fc(i)%in
   out = fc(i)%out
   if(in/=0) then
     cell(in)%nfc = cell(in)%nfc + 1
   endif
   if(out/=0) then
     cell(out)%nfc = cell(out)%nfc + 1
   endif
enddo

do i=1,nc
   allocate(cell(i)%fc(cell(i)%nfc))
enddo

cell(:)%nfc = 0
do i=1,nf
   in = fc(i)%in
   out = fc(i)%out
   if(in/=0) then
     cell(in)%nfc = cell(in)%nfc + 1
     cell(in)%fc(cell(in)%nfc) = i
   endif
   if(out/=0) then
     cell(out)%nfc = cell(out)%nfc + 1
     cell(out)%fc(cell(out)%nfc) = i
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
integer(kind=4) :: i,j,in,out,count(nc),c1,c2
integer(kind=4) :: n1,n2,nop,noc,nof
integer(kind=4) :: p1,p2,e1,e2,c
real(kind=8)    :: dx,dy,dxc,dyc,ds,x1,y1,x2,y2,nx,ny
real(kind=8)    :: rx,ry,dotprod,xc,yc
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

!Cells to vertex connectivity 





! cells surrounding a node
print*
print*,"==> Cells surrounding a node"

node(:)%cn=0
do i=1,nof
   p1=face(i)%pt(1)
   p2=face(i)%pt(2)
   e1=face(i)%in
   e2=face(i)%out

   if(e1/=0) then
         c=0
         do j=1,node(p1)%cn
            if(node(p1)%cnbr(j)==e1) c=c+1  
         enddo
         if(c==0) then 
           node(p1)%cn=node(p1)%cn+1
           call alloc_int_ptr(node(p1)%cnbr,node(p1)%cn)      
           node(p1)%cnbr(node(p1)%cn)=e1      
         endif
      
         c=0
         do j=1,node(p2)%cn
            if(node(p2)%cnbr(j)==e1) c=c+1  
         enddo
         if(c==0) then 
           node(p2)%cn=node(p2)%cn+1
           call alloc_int_ptr(node(p2)%cnbr,node(p2)%cn)      
           node(p2)%cnbr(node(p2)%cn)=e1      
         endif
   endif

   if(e2/=0) then
         c=0
         do j=1,node(p1)%cn
            if(node(p1)%cnbr(j)==e2) c=c+1  
         enddo
         if(c==0) then 
           node(p1)%cn=node(p1)%cn+1
           call alloc_int_ptr(node(p1)%cnbr,node(p1)%cn)      
           node(p1)%cnbr(node(p1)%cn)=e2      
         endif
      
         c=0
         do j=1,node(p2)%cn
            if(node(p2)%cnbr(j)==e2) c=c+1  
         enddo
         if(c==0) then 
           node(p2)%cn=node(p2)%cn+1
           call alloc_int_ptr(node(p2)%cnbr,node(p2)%cn)      
           node(p2)%cnbr(node(p2)%cn)=e2      
         endif
   endif

enddo   


!do i=1,nop
!   print*,i,(node(i)%cnbr(j), j=1,node(i)%cn) 
!enddo

print*, 'Max. nbr=',maxval(node(:)%cn),'at node',maxloc(node(:)%cn)
print*, 'Min. nbr=',minval(node(:)%cn),'at node',minloc(node(:)%cn)

! cells surrounding a cell
print*
print*,"==>Cells surrounding a cell "

elem(:)%cn=0
do i=1,nof
   e1=face(i)%in
   e2=face(i)%out
 
   if(e1/=0.and.e2/=0) then

     c=0
     do j=1,elem(e1)%cn
        if(elem(e1)%cnbr(j)==e2) c=c+1  
     enddo

     if(c==0) then 
       elem(e1)%cn=elem(e1)%cn+1
       call alloc_int_ptr(elem(e1)%cnbr,elem(e1)%cn)      
       elem(e1)%cnbr(elem(e1)%cn)=e2      
     endif


     c=0
     do j=1,elem(e2)%cn
        if(elem(e2)%cnbr(j)==e1) c=c+1  
     enddo

     if(c==0) then 
       elem(e2)%cn=elem(e2)%cn+1
       call alloc_int_ptr(elem(e2)%cnbr,elem(e2)%cn)      
       elem(e2)%cnbr(elem(e2)%cn)=e1      
     endif


   endif
   
enddo


!do i=1,noc
!   print*,i,(elem(i)%cnbr(j), j=1,elem(i)%cn) 
!enddo

print*, 'Max. cell nbr=',maxval(elem(:)%cn),'at node',maxloc(elem(:)%cn)
print*, 'Min. cell nbr=',minval(elem(:)%cn),'at node',minloc(elem(:)%cn)


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


end subroutine connectivity


