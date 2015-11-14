module  data_type
implicit none
integer, parameter :: i1=selected_int_kind(2)
integer, parameter :: i4=selected_int_kind(4)
integer, parameter :: i2=selected_int_kind(9)
integer, parameter :: i8=selected_int_kind(18)
integer, parameter :: sp=selected_real_kind(6,37)
integer, parameter :: dp=selected_real_kind(15,307)
integer, parameter :: qp=selected_real_kind(31,307)

end module data_type

module grid
use data_type
implicit none

integer(kind=i4):: np,nf,nc 

type points 
     real(kind=dp) :: x,y,z
     integer(kind=i4):: bc,flag,nv2c 
     integer(kind=i4),dimension(:),pointer::v2c
end type points 

type faces  
     integer(kind=i4):: pt(2) 
     integer(kind=i4):: in 
     integer(kind=i4):: out 
     integer(kind=i4):: bc
     integer(kind=i4):: flag
     real(kind=dp)   :: sx,sy 
end type faces  

type cells  
     integer(kind=i4):: nc2v,nc2f,nc2c  
     real(kind=dp)   :: xc,yc,zc,cv
     integer(kind=i4),dimension(:),pointer::c2v
     integer(kind=i4),dimension(:),pointer::c2f
     integer(kind=i4),dimension(:),pointer::c2c
end type cells  

end module grid
!========================================================
program Connectivity
use grid
implicit none

integer(kind=i4) :: i,j,k
integer(kind=i4) :: n1,n2,m1,m2,nop,noc,nof
integer(kind=i4) :: p1,p2,e1,e2,c,f1
real(kind=dp)    :: x1,y1,x2,y2
real(kind=dp)    :: xc,yc
real(kind=dp)    :: check

type(points),allocatable,dimension(:)::node
type(faces),allocatable,dimension(:)::face
type(cells),allocatable,dimension(:)::elem

open(3,file='geometry.inp')
read(3,*)nop,noc,nof
allocate(node(nop),face(nof),elem(noc))
do i=1,nop
 read(3,*)j,node(i)%x,node(i)%y,node(i)%bc
enddo
do i=1,noc
 read(3,*)j,elem(i)%xc,elem(i)%yc,elem(i)%cv
enddo
do i=1,nof
 read(3,*)j,face(i)%pt(1),face(i)%pt(2),face(i)%in,face(i)%out,face(i)%sx,face(i)%sy,face(i)%bc
enddo
close(3)


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
       call alloc_int_ptr(elem(e1)%c2c,elem(e1)%nc2c,1_i4)      
       elem(e1)%c2c(elem(e1)%nc2c)=e2      
     endif


     c=0
     do j=1,elem(e2)%nc2c
        if(elem(e2)%c2c(j)==e1) c=c+1  
     enddo

     if(c==0) then 
       elem(e2)%nc2c=elem(e2)%nc2c+1
       call alloc_int_ptr(elem(e2)%c2c,elem(e2)%nc2c,1_i4)      
       elem(e2)%c2c(elem(e2)%nc2c)=e1      
     endif


   endif
   
enddo


print*,'cell'
do i=1,noc
   print*,i,(elem(i)%c2c(j), j=1,elem(i)%nc2c) 
enddo


call renumber

do i=1,noc
   print*,i,(elem(i)%c2c(j), j=1,elem(i)%nc2c) 
enddo

print*, 'Max. cell nbr=',maxval(elem(:)%nc2c),'at node',maxloc(elem(:)%nc2c)
print*, 'Min. cell nbr=',minval(elem(:)%nc2c),'at node',minloc(elem(:)%nc2c)

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
           call alloc_int_ptr(node(p1)%v2c,node(p1)%nv2c,1_i4)      
           node(p1)%v2c(node(p1)%nv2c)=e1      
         endif
      
         c=0
         do j=1,node(p2)%nv2c
            if(node(p2)%v2c(j)==e1) c=c+1  
         enddo
         if(c==0) then 
           node(p2)%nv2c=node(p2)%nv2c+1
           call alloc_int_ptr(node(p2)%v2c,node(p2)%nv2c,1_i4)      
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
           call alloc_int_ptr(node(p1)%v2c,node(p1)%nv2c,1_i4)      
           node(p1)%v2c(node(p1)%nv2c)=e2      
         endif
      
         c=0
         do j=1,node(p2)%nv2c
            if(node(p2)%v2c(j)==e2) c=c+1  
         enddo
         if(c==0) then 
           node(p2)%nv2c=node(p2)%nv2c+1
           call alloc_int_ptr(node(p2)%v2c,node(p2)%nv2c,1_i4)      
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
   call alloc_int_ptr(elem(e1)%c2f,elem(e1)%nc2f,1_i4)      
   elem(e1)%c2f(elem(e1)%nc2f)=i      
   endif

   if(e2/=0) then
   elem(e2)%nc2f=elem(e2)%nc2f+1
   call alloc_int_ptr(elem(e2)%c2f,elem(e2)%nc2f,1_i4)      
   elem(e2)%c2f(elem(e2)%nc2f)=i      
   endif
enddo

print*, 'Max. cell faces =',maxval(elem(:)%nc2f),'at node',maxloc(elem(:)%nc2f)
print*, 'Min. cell faces =',minval(elem(:)%nc2f),'at node',minloc(elem(:)%nc2f)

!Cells to vertex connectivity 
print*
print*,"==> Vertex surrounding a cell "

!elem(:)%nc2v=0
!do i=1,nop
!   do j=1,node(i)%nv2c
!      e1=node(i)%v2c(j)
!      elem(e1)%nc2v=elem(e1)%nc2v+1
!      call alloc_int_ptr(elem(e1)%c2v,elem(e1)%nc2v,1_i4)      
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
         call alloc_int_ptr(elem(i)%c2v,elem(i)%nc2v,1_i4)      
         elem(i)%c2v(elem(i)%nc2v)=n1      

         elem(i)%nc2v=elem(i)%nc2v+1
         call alloc_int_ptr(elem(i)%c2v,elem(i)%nc2v,1_i4)      
         elem(i)%c2v(elem(i)%nc2v)=n2      
      else
         elem(i)%nc2v=elem(i)%nc2v+1
         call alloc_int_ptr(elem(i)%c2v,elem(i)%nc2v,1_i4)      
         elem(i)%c2v(elem(i)%nc2v)=n2      

         elem(i)%nc2v=elem(i)%nc2v+1
         call alloc_int_ptr(elem(i)%c2v,elem(i)%nc2v,1_i4)      
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
               call alloc_int_ptr(elem(i)%c2v,elem(i)%nc2v,1_i4)      
               elem(i)%c2v(elem(i)%nc2v)=m2      
            else
               elem(i)%nc2v=elem(i)%nc2v+1
               call alloc_int_ptr(elem(i)%c2v,elem(i)%nc2v,1_i4)      
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
   print*,i
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


print*
print*

!===========================================================================
!                      Allocate/extend array
!===========================================================================

contains

subroutine alloc_int_ptr(x,n,key)
implicit none
integer(kind=i4),intent(in) ::n 
integer(kind=i4)::i,key 
integer(kind=i4),dimension(:),pointer::temp
integer(kind=i4),dimension(:),pointer::x

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

!deallocate(temp)

return

end subroutine alloc_int_ptr

function dotprod(n1,n2,xc,yc)
   implicit none
   integer(kind=i4):: n1,n2
   real(kind=dp)   :: x1,y1,x2,y2
   real(kind=dp)   :: dx,dy,ds,xc,yc
   real(kind=dp)   :: dotprod,nx,ny,rx,ry


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
   integer(kind=i4):: n1,n2,m1,m2
   real(kind=dp)   :: x1,y1,x2,y2
   real(kind=dp)   :: x3,y3,x4,y4
   real(kind=dp)   :: dx,dy,ds
   real(kind=dp)   :: crossprod,nx,ny,rx,ry


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





subroutine  renumber
implicit none
integer(kind=i4) :: i,j,c,it,t1
integer(kind=i4),allocatable:: oldnum(:), newnum(:)

type(cells),allocatable,dimension(:)::elmn


allocate(oldnum(noc),newnum(noc))
allocate(elmn(noc))


print*,'noc=',noc
do i=1,noc
   oldnum(i) = 0
   newnum(i) = 0
enddo

oldnum(1)=1
newnum(1)=1
c=1


do i=1,noc
   it=oldnum(i)
   if(it .eq. 0)then
      print*,'renumber: Fatal error. it is zero for i=',i
      stop
   endif

   do j=1,elem(it)%nc2c
      t1=elem(it)%c2c(j)
      if(t1.gt.0)  then 
      if(newnum(t1).eq.0) then
        c=c+1
        oldnum(c)=t1   
        newnum(t1)=c   
      endif 
      endif 
   enddo 

enddo
if(noc .ne. c)then
  print*,'renumber: count does not match nt.'
  print*,'          Possible bug'
  stop
endif


do i=1,noc
   elmn(i)%nc2c=0
   elmn(i)%xc=elem(i)%xc
   elmn(i)%yc=elem(i)%yc
   elmn(i)%cv=elem(i)%cv
   do j=1,elem(i)%nc2c
      elmn(i)%nc2c=elmn(i)%nc2c+1
      call alloc_int_ptr(elmn(i)%c2c,elmn(i)%nc2c,1_i4)      
      elmn(i)%c2c(elmn(i)%nc2c)=elem(i)%c2c(j)      
   enddo 
enddo

!do i=1,noc
!   print*,'old',i,(elmn(i)%c2c(j),j=1,elmn(i)%nc2c)
!   print*,'new',i,(newnum(elmn(i)%c2c(j)),j=1,elmn(i)%nc2c)
!enddo
do i=1,nof
   e1=face(i)%in
   e2=face(i)%out

   if(e1/=0) then
     face(i)%in=newnum(e1)
   endif

   if(e2/=0) then
     face(i)%out=newnum(e2)
   endif
enddo    
 



do i=1,noc
   it=oldnum(i)
   elem(i)%nc2c=0
   elem(i)%xc=elmn(it)%xc
   elem(i)%yc=elmn(it)%yc
   elem(i)%cv=elmn(it)%cv
   do j=1,elmn(it)%nc2c
      t1=elmn(it)%c2c(j)
      if(t1.gt.0) then  
         elem(i)%nc2c=elem(i)%nc2c+1
         call alloc_int_ptr(elem(i)%c2c,elem(i)%nc2c,1_i4)      
         elem(i)%c2c(elem(i)%nc2c)=newnum(t1)
      else
         elem(i)%nc2c=elem(i)%nc2c+1
         call alloc_int_ptr(elem(i)%c2c,elem(i)%nc2c,1_i4)      
         elem(i)%c2c(elem(i)%nc2c)=t1
      endif
   enddo 
enddo

end subroutine  renumber

end program Connectivity
