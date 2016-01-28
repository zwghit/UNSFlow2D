subroutine solid_flux(ie,cc)
use commons
use pri
use grid
implicit none
integer(kind=i4):: ie,cc
real(kind=dp) :: qc(nvar)

real(kind=dp) :: nx, ny,con(nvar)

!con(:)=cell(cc)%qc(:)
!call con2prim(con)

!nx=fc(ie)%ldx
!ny=fc(ie)%ldy
!if(fc(ie)%ldx==0.d0.and.fc(ie)%ldy==0.d0) then
!nx=fc(ie)%rdx
!ny=fc(ie)%rdy
!endif

con(:)=cell(cc)%qp(:)!+(cell(cc)%qx(:)*nx+cell(cc)%qy(:)*ny)

p=con(4)

nx=fc(ie)%sx
ny=fc(ie)%sy
!ds=dsqrt(nx*nx+ny*ny)

cell(cc)%res(3) =cell(cc)%res(3)+p*nx
cell(cc)%res(4) =cell(cc)%res(4)+p*ny
!print*,'wall',ie,cc
end
