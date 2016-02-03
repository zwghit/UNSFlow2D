module grid
implicit none

integer(kind=4):: np,nf,nc 

type points 
     real(kind=8) :: x,y,z
     integer(kind=4):: bc,flag,nv2c 
     integer(kind=4),pointer,dimension(:)::v2c
end type points 

type faces  
     integer(kind=4):: pt(2) 
     integer(kind=4):: in 
     integer(kind=4):: out 
     real(kind=8)   :: sx,sy 
     integer(kind=4):: bc
end type faces  

type cells  
     integer(kind=4):: nc2v,nc2f,nc2c  
     real(kind=8)   :: xc,yc,zc,cv
     integer(kind=4),pointer,dimension(:)::c2v
     integer(kind=4),pointer,dimension(:)::c2f
     integer(kind=4),pointer,dimension(:)::c2c
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
integer(kind=4) :: i,j
integer(kind=4) :: nop,noc,nof
integer(kind=4) :: e1,e2,c

type(points),allocatable,dimension(:)::node
type(faces),allocatable,dimension(:)::face
type(cells),allocatable,dimension(:)::elem,elm
type(cells),allocatable,dimension(:)::elms


nop=32
noc=24
nof=56

open(3,file='geometry.inp')
read(3,*)nop,noc,nof
allocate(node(nop),face(nof),elem(noc),elm(noc))
!allocate(elms(noc))

print*,nop,noc,nof
do i=1,nop
  read(3,*)j,node(i)%x,node(i)%y,node(i)%bc
enddo
do i=1,noc
 read(3,*)j,elem(i)%xc,elem(i)%yc,elem(i)%cv
enddo
do i=1,nof
 read(3,*)j,face(i)%pt(1),face(i)%pt(2),face(i)%in,face(i)%out,face(i)%sx,face(i)%sy,face(i)%bc
! print*,j,face(i)%pt(1),face(i)%pt(2),face(i)%in,face(i)%out,face(i)%sx,face(i)%sy,face(i)%bc
enddo
close(3)




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



elm(:)%nc2c=0
do i=1,nof
   e1=face(i)%in
   e2=face(i)%out
 
   if(e1/=0.and.e2/=0) then

     c=0
     do j=1,elm(e1)%nc2c
        if(elm(e1)%c2c(j)==e2) c=c+1  
     enddo

     if(c==0) then 
       elm(e1)%nc2c=elm(e1)%nc2c+1
       call alloc_int_ptr(elm(e1)%c2c,elm(e1)%nc2c)      
       elm(e1)%c2c(elm(e1)%nc2c)=e2      
     endif


     c=0
     do j=1,elm(e2)%nc2c
        if(elm(e2)%c2c(j)==e1) c=c+1  
     enddo

     if(c==0) then 
       elm(e2)%nc2c=elm(e2)%nc2c+1
       call alloc_int_ptr(elm(e2)%c2c,elm(e2)%nc2c)      
       elm(e2)%c2c(elm(e2)%nc2c)=e1      
     endif


   endif
   
enddo


allocate(elms(noc))


elms(:)%nc2c=0
do i=1,nof
   e1=face(i)%in
   e2=face(i)%out
 
   if(e1/=0.and.e2/=0) then

     c=0
     do j=1,elms(e1)%nc2c
        if(elms(e1)%c2c(j)==e2) c=c+1  
     enddo

     if(c==0) then 
       elms(e1)%nc2c=elms(e1)%nc2c+1
       call alloc_int_ptr(elms(e1)%c2c,elms(e1)%nc2c)      
       elms(e1)%c2c(elms(e1)%nc2c)=e2      
     endif


     c=0
     do j=1,elms(e2)%nc2c
        if(elms(e2)%c2c(j)==e1) c=c+1  
     enddo

     if(c==0) then 
       elms(e2)%nc2c=elms(e2)%nc2c+1
       call alloc_int_ptr(elms(e2)%c2c,elms(e2)%nc2c)      
       elms(e2)%c2c(elms(e2)%nc2c)=e1      
     endif


   endif
   
enddo


contains


!===========================================================================
!                      Allocate/extend array
!===========================================================================

subroutine alloc_int_ptr(x,n)
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

!deallocate(temp)

return

end subroutine alloc_int_ptr

end program SolverMesh
