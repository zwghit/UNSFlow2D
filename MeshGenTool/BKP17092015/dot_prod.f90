program dotproduct
implicit none

integer*4:: i,j,n1,n2,cc,m1,m2,temp
real*8   :: x1,y1,x2,y2,check    
real*8   :: dx,dy,ds,xc,yc
real*8   :: nx,ny,rx,ry 

type points
     real(kind=8) :: x,y
end type points

type faces
     integer(kind=4):: pt(2)
end type faces

type(points),dimension(:),pointer::pt,ptn
type(faces),dimension(:),pointer::fc


allocate(pt(4),fc(4),ptn(4))

pt(1)%x=1.0 ; pt(1)%y=1.0
pt(2)%x=1.5 ; pt(2)%y=1.0
pt(3)%x=2.5 ; pt(3)%y=2.0
pt(4)%x=1.1 ; pt(4)%y=2.3

fc(1)%pt(1)=1 ; fc(1)%pt(2)=2
fc(2)%pt(1)=2 ; fc(2)%pt(2)=3
fc(3)%pt(1)=4 ; fc(3)%pt(2)=3
fc(4)%pt(1)=4 ; fc(4)%pt(2)=1

open(3,file='cell.dat')
do i=1,4
   n1=fc(i)%pt(1)
   n2=fc(i)%pt(2)
   x1 = pt(n1)%x ; y1 = pt(n1)%y
   x2 = pt(n2)%x ; y2 = pt(n2)%y
   write(3,100)x1,y1,x2-x1,y2-y1
   xc=xc+pt(n1)%x
   yc=yc+pt(n1)%y
enddo
close(3)

xc=xc/4.0
yc=yc/4.0


100 format(1x,4(f15.6,1x))

n1=fc(i)%pt(1)
n2=fc(i)%pt(2)
check=dotprod(n1,n2,xc,yc) 
cc=1
if(check > 0.d0) then
  ptn(cc)=pt(fc(1)%pt(2))
  temp=fc(1)%pt(1)
else
  ptn(cc)=pt(fc(1)%pt(1))
  temp=fc(1)%pt(2)
endif

do i=2,4
   cc=cc+1
   n1=fc(i)%pt(1)
   n2=fc(i)%pt(2)
   check=dotprod(n1,n2,xc,yc) 
 
   if(check > 0.d0) then
      print*, i ,check,'CW' 
  
   elseif(check   < 0.d0) then
      print*, i ,check,'ACW' 
   else
      print*,i,'Collinear'
   endif
end do


do i=1,3
   n1=fc(i)%pt(1)
   n2=fc(i)%pt(2)
   m1=fc(i+2)%pt(1)
   m2=fc(i+2)%pt(2)
   check=crossprod(n1,n2,m1,m2) 
 
   if(check > 0.d0) then
      print*, i ,check,'CW',check 
  
   elseif(check   < 0.d0) then
      print*, i ,check,'ACW' ,check
   else
      print*,i,'Collinear',check
   endif
enddo




contains

function dotprod(n1,n2,xc,yc)
   implicit none
   integer*4:: i,j,n1,n2
   real*8   :: x1,y1,x2,y2    
   real*8   :: dx,dy,ds,xc,yc
   real*8   :: dotprod,nx,ny,rx,ry 


   x1 = pt(n1)%x ; y1 = pt(n1)%y
   x2 = pt(n2)%x ; y2 = pt(n2)%y
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


   x1 = pt(n1)%x ; y1 = pt(n1)%y
   x2 = pt(n2)%x ; y2 = pt(n2)%y

   x3 = pt(m1)%x ; y3 = pt(m1)%y
   x4 = pt(m2)%x ; y4 = pt(m2)%y

   dx = x2-x1 ; dy = y2-y1
   ds = dsqrt(dx*dx+dy*dy)
   nx = dy/ds ; ny = -dx/ds

   dx = x4-x3 ; dy = y4-y3
   ds = dsqrt(dx*dx+dy*dy)
   rx = dx/ds
   ry = dy/ds

   crossprod = nx*ry-ny*rx

end function crossprod
end program dotproduct
