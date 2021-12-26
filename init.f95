module init

use grid
use mpi_interface

implicit none

contains

subroutine init_grids

implicit none

integer :: i,j

do i=1,N_vert
	zt(i)=hz*(dble(i)-0.5)
enddo

if(outcoords(2) > 0)then
	do i=1,x_dim
		xt(i)=-dL+(N_x-(outdims(2)-1)*(N_x/outdims(2)))*hx+x_dim*(outcoords(2)-1)*hx+hx*(dble(i)-0.5)
	enddo
else
	do i=1,x_dim
		xt(i)=-dL+hx*(dble(i)-0.5)
	enddo
endif

if(outcoords(1)>0)then
	do i=1,y_dim
		yt(i)=-dL+(N_y-(outdims(1)-1)*(N_y/outdims(1)))*hy+y_dim*(outcoords(1)-1)*hy+hy*(dble(i)-0.5)
	enddo
else
	do i=1,y_dim
		yt(i)=-dL+hy*(dble(i)-0.5)
	enddo
endif

end subroutine init_grids

subroutine init_heating(S,M,N,L)

implicit none

real*8, dimension(M+2,N+2,L+2), intent(out) :: S
integer, intent(in) :: M,N,L
real*8 :: r
integer :: i,j,k

do i=2,M+1
	do j=2,N+1
		do k=2,L+1
			r=sqrt(xt(j-1)**2 + yt(k-1)**2)
			if(r<drad)then
				S(i,j,k)=sin(4.0*atan(1.0)*zt(i-1))
			!	S(i,j,k)=sin(4.0*atan(1.0)*zt(i-1))*sin(4.0*atan(1.0)*xt(j-1))*sin(4.0*atan(1.0)*yt(k-1))
			else
				S(i,j,k)=0.0
			endif
		enddo
	enddo
enddo


end subroutine init_heating

end module


