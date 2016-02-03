program mesh
implicit none 
integer(kind=4):: i,j,k,m,n,nz
real(kind=8)   :: z,dummy
real(kind=8),allocatable,dimension(:,:)::x,y

open(3,file='naca0012.mesh-00257.2d.x',form='unformatted')
!open(3,file='naca0012.mesh-00009.2d.x')
READ(3)m,n
print*,m,n
allocate(x(m,n),y(m,n))
READ(3)((X(I,J), I=1,m), J=1,N),&
       ((Y(I,J), I=1,m), J=1,N)

nz=1
dummy=0.d0
open(3,file='grid.grd')
write(3,*)1
write(3,*) m , n , nz
write(3,*) (((x(i,j), i=1,m), j=1,n), k=1,nz), &
          (((y(i,j), i=1,m), j=1,n), k=1,nz), &
          (((dummy , i=1,m), j=1,n), k=1,nz)
close(3)

end program mesh
