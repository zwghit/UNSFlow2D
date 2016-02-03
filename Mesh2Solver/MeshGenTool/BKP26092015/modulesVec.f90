module class_Vector
implicit none
! filename: class_Vector.inc
! public, everything by default, but can specify any
type Vector
private
      integer :: size ! vector length
      real, pointer, dimension(:) :: data ! component values
end type Vector

! Overload common operators

! add others later
interface operator (+) 
    module procedure add_Vector, add_Real_to_Vector 
end interface

! add unary versions later
interface operator (-) 
    module procedure subtract_Vector, subtract_Real 
end interface

! overload *
interface operator (*) 
    module procedure dot_Vector, real_mult_Vector, Vector_mult_real
end interface

! overload =
interface assignment (=) 
    module procedure equal_Real 
end interface

! overload ==
interface operator (==) 
   module procedure is_equal_to  
end interface

! functions & operators
contains 

! overload +
function add_Real_to_Vector (v, r) result (new) 
implicit none
type (Vector), intent(in) :: v
real, intent(in) :: r
type (Vector) :: new ! new = v + r
     if ( v%size < 1 ) stop "No sizes in add_Real_to_Vector"
          allocate ( new%data(v%size) )  
                     new%size = v%size
     ! new%data = v%data + r ! as array operation, or use implied loop
       new%data(1:v%size) = v%data(1:v%size) + r  
end function

! vector + vector
function add_Vector (a, b) result (new) 
implicit none
type (Vector), intent(in) :: a, b
type (Vector) :: new ! new = a + b
     if ( a%size /= b%size ) stop "Sizes differ in add_Vector"
     allocate ( new%data(a%size) )  
     new%size = a%size
     new%data = a%data + b%data 
end function add_Vector

! array to vector constructor
function assign (values) result (name) 
implicit none
real, intent(in) :: values(:) ! given rank 1 array
integer :: length ! array size
type (Vector) :: name ! Vector to create
     length = size(values)
     allocate ( name%data(length) )
     name % size = length
     name % data = values
end function assign


function copy_Vector (name) result (new)
implicit none
type (Vector), intent(in) :: name
type (Vector) :: new
      allocate ( new%data(name%size) ) 
      new%size = name%size
      new%data = name%data 
end function copy_Vector

! deallocate allocated items
subroutine delete_Vector (name) 
implicit none
type (Vector), intent(inout) :: name
integer :: ok ! check deallocate status
      deallocate (name%data, stat = ok )
      if ( ok /= 0 ) stop "Vector not allocated in delete_Vector"
      name%size = 0 
end subroutine delete_Vector

! overload *
function dot_Vector (a, b) result (c) 
implicit none
type (Vector), intent(in) :: a, b
real :: c
     if ( a%size /= b%size ) stop "Sizes differ in dot_Vector"
     c = dot_product (a%data, b%data) 
end function dot_Vector

! overload =, real to vector
subroutine equal_Real (new, R) 
implicit none
type (Vector), intent(inout) :: new
real, intent(in) :: R
     if ( associated (new%data) ) deallocate (new%data)
     allocate ( new%data(1) ) 
     new%size = 1
     new%data = R 
end subroutine equal_Real

! overload ==
logical function is_equal_to (a, b) result (t_f) 
implicit none
type (Vector), intent(in) :: a, b ! left & right of ==
      t_f = .false. ! initialize
      if ( a%size /= b%size ) return ! same size ?
      t_f = all ( a%data == b%data ) ! and all values match
end function is_equal_to

! accessor member
function length (name) result (n) 
implicit none
type (Vector), intent(in) :: name
integer :: n
      n = name % size 
end function length

! accessor member
subroutine list (name) 
implicit none
type (Vector), intent(in) :: name
     print *,"[", name % data(1:name%size), "]"
end subroutine list

! Optional Constructor
function make_Vector (len, values) result(v) 
implicit none
integer, optional, intent(in) :: len ! number of values
real, optional, intent(in) :: values(:) ! given values
type (Vector) :: v
      if ( present (len) ) then ! create vector data
          v%size = len 
          allocate ( v%data(len) )
          if ( present (values)) then 
              v%data = values ! vector
          else 
              v%data = 0.d0 ! null vector
          end if ! values present
      else ! scalar constant
          v%size = 1 
          allocate ( v%data(1) ) ! default
          if ( present (values)) then 
             v%data(1) = values(1) ! scalar
          else 
             v%data(1) = 0.d0 ! null
          end if ! value present
      end if ! len present
end function make_Vector


function normalize_Vector (name) result (new)
implicit none
type (Vector), intent(in) :: name
type (Vector) :: new
real :: total, nil = epsilon(nil) ! tolerance
     allocate ( new%data(name%size) ) 
     new%size = name%size
     total = dsqrt ( sum ( name%data**2 ) ) ! intrinsic functions
     if ( total < nil ) then
         new%data = 0.d0 ! avoid division by 0
     else 
         new%data = name%data/total
     end if 
end function normalize_Vector

! read array, assign
subroutine read_Vector (name) 
implicit none
type (Vector), intent(inout) :: name
integer, parameter :: max = 999
integer :: length
       read (*,'(i1)', advance = 'no') length
       if ( length <= 0 ) stop "Invalid length in read_Vector"
       if ( length >= max ) stop "Maximum length in read_Vector"
       allocate ( name % data(length) ) 
       name % size = length
       read *, name % data(1:length) 
end subroutine read_Vector

 ! overload *
function real_mult_Vector (r, v) result (new)
implicit none
real, intent(in) :: r
type (Vector), intent(in) :: v
type (Vector) :: new ! new = r * v
     if ( v%size < 1 ) stop "Zero size in real_mult_Vector"
     allocate ( new%data(v%size) ) 
     new%size = v%size
     new%data = r * v%data 
end function real_mult_Vector

! accessor member
function size_Vector (name) result (n) 
implicit none
type (Vector), intent(in) :: name
integer :: n
        n = name % size 
end function size_Vector

! vector-real, overload -
function subtract_Real (v, r) result (new) 
implicit none
type (Vector), intent(in) :: v
real, intent(in) :: r
type (Vector) :: new ! new = v + r
     if ( v%size < 1 ) stop "Zero length in subtract_Real"
     allocate ( new%data(v%size) ) 
     new%size = v%size
     new%data = v%data - r 
end function subtract_Real

! overload -
function subtract_Vector (a, b) result (new) 
implicit none
type (Vector), intent(in) :: a, b
type (Vector) :: new
     if ( a%size /= b%size ) stop "Sizes differ in subtract_Vector"
     allocate ( new%data(a%size) ) 
     new%size = a%size
     new%data = a%data - b%data 
end function subtract_Vector

! accessor member
function values (name) result (array)
implicit none
type (Vector), intent(in) :: name
real :: array(name%size)
      array = name % data 
end function values

! Public constructor
function Vector_ (length, values) result(name) 
implicit none
integer, intent(in) :: length ! array size
real, target, intent(in) :: values(length) ! given array
real, pointer :: pt_to_val(:) ! pointer to array
type (Vector) :: name ! Vector to create
integer :: get_m ! allocate flag
     allocate ( pt_to_val (length), stat = get_m ) ! allocate
     if ( get_m /= 0 ) stop 'allocate error' ! check
     pt_to_val = values ! dereference values
     name = Vector(length, pt_to_val) ! intrinsic constructor
end function Vector_


function Vector_max_value (a) result (v) ! accessor member
implicit none
type (Vector), intent(in) :: a
real :: v
    v = maxval ( a%data(1:a%size) )
end function Vector_max_value


function Vector_min_value (a) result (v) ! accessor member
implicit none
type (Vector), intent(in) :: a
real :: v
     v = minval ( a%data(1:a%size) ) 
end function Vector_min_value


function Vector_mult_real (v, r) result (new) ! vector*real, overload *
implicit none
type (Vector), intent(in) :: v
real, intent(in) :: r
type (Vector) :: new ! new = v * r
     if ( v%size < 1 ) stop "Zero size in Vector_mult_real"
     new = Real_mult_Vector (r, v) 
end function Vector_mult_real

end module class_Vector


! Testing Vector Class Constructors & Operators
!include 'class_Vector.f90' ! see previous figure
program check_vector_class
use class_Vector
implicit none
type (Vector) :: x, y, z

! test optional constructors: assign, and copy

x = make_Vector () ! single scalar zero
write (*,'("made scalar x = ")', advance='no')
call list (x)
call delete_Vector (x) 

y = make_Vector (4) ! 4 zero values
write (*,'("made null y = ")', advance='no')
call list (y)

z = make_Vector (4, (/11., 12., 13., 14./) ) ! 4 non-zero values
write (*,'("made full z = ")', advance='no')
call list (z)
write (*,'("assign [ 31., 32., 33., 34. ] to x")')

x = assign( (/31., 32., 33., 34./) ) ! (4) non-zeros
write (*,'("assigned x = ")', advance='no')
call list (x)
x = Vector_(4, (/31., 32., 33., 34./) ) ! 4 non-zero values
write (*,'("public x = ")', advance='no')
call list (x)
write (*,'("copy x to y =")', advance='no')
y = copy_Vector (x)  
call list (y) ! copy

! test overloaded operators
write (*,'("z * x gives ")', advance='no')
print *, z*x ! dot
write (*,'("z + x gives ")', advance='no')
call list (z+x) ! add
y = 25.6 ! real to vector
write (*,'("y = 25.6 gives ")', advance='no')
call list (y)
y = z ! equality
write (*,'("y = z gives y as ")', advance='no')
call list (y)
write (*,'("logic y == x gives ")', advance='no')
print *, y==x
write (*,'("logic y == z gives ")', advance='no')
print *, y==z

! test destructor, accessors
call delete_Vector (y) ! destructor
write (*,'("deleting y gives y = ")', advance='no'); call list (y)
print *, "size of x is ", length (x) ! accessor
print *, "data in x are [", values (x), "]" ! accessor
write (*,'("2. times x is ")', advance='no'); call list (2.0*x)
write (*,'("x times 2. is ")', advance='no'); call list (x*2.0)
call delete_Vector (x); call delete_Vector (z) ! clean up
end program check_vector_class


! Running gives the output: 
! made scalar x = [0.]
! made null y = [0., 0., 0., 0.] 
!made full z = [11., 12., 13., 14.]
! assign [31., 32., 33., 34.] to x 
!assigned x = [31., 32., 33., 34.]
! public x = [31., 32., 33., 34.] 
! copy x to y = [31., 32., 33., 34.]
! z * x gives 1630. 
! z + x gives [42., 44., 46., 48.]
! y = 25.6 gives [25.6000004] 
! y = z, y = [11., 12., 13., 14.]
! logic y == x gives F 
! logic y == z gives T
! deleting y gives y = [] 
! size of x is 4
! data in x : [31., 32., 33., 34.] 
! 2. times x is [62., 64., 66., 68.]
! x times 2. is [62., 64., 66., 68.]
