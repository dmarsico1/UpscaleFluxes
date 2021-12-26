module grid

implicit none

real*8, allocatable, dimension(:,:,:,:) :: var,write_array,f !f is the array of forcings for each equation
real*8, allocatable, dimension(:,:,:) :: g
real*8, allocatable, dimension(:,:,:) :: S !S is the heating
real*8, allocatable, dimension(:) :: xt,zt,yt

integer :: N_vert=50
integer :: N_x=500
integer :: N_y=500
integer :: NGx=1
integer :: NGy=1
integer :: NGz=1
integer :: x_dim,y_dim

character(len=305) :: file_name='grid1'

real*8 :: hz !vertical grid spacing
real*8 :: hx !zonal grid spacing
real*8 :: hy !meridional grid spacing

real*8 :: dtheta=0
real*8 :: dalpha=5.0
real*8 :: dbeta=0.1
real*8 :: dL=10.0
real*8 :: dLz=1.0
real*8 :: drad=0.25
real*8 :: tol=0.5*10.0**(-2)

contains

subroutine define_name_vars

implicit none

namelist /model/ &
	N_x,N_y,N_vert,file_name,dtheta,dalpha,dbeta,drad
open(unit=1,file='NAMELIST')
read(1,nml=model)
 close(1)

end subroutine define_name_vars


end module
