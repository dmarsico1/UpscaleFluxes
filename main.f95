program main

use grid
use mpi_interface
use netcdf_write
use init
use matrix

implicit none

integer :: i

call define_name_vars

hx=2.0*dL/N_x
hy=2.0*dL/N_y
hz=dLz/N_vert

call mpi_begin

if(myid==0)then
	call init_netcdf_vars
endif

call init_grids

call neighbors

call get_ranks

if(myid==0)then
	allocate(write_array(N_vert,N_x,N_y,7))
endif

S=0.0

!initialize heating
call init_heating(S,N_vert,x_dim,y_dim)

!set forcing for g equation
f(:,:,:,1)=S(2:N_vert+1,2:x_dim+1,2:y_dim+1)

!solve for g
call ConjGrad(g,f(:,:,:,1),tol,x_dim,N_vert,y_dim,hx,hy,hz)

!exchange boundary data
call v_exchange(g,1)

!set upper and lower boundary conditions in z
g(1,2:x_dim+1,2:y_dim+1)=-g(2,2:x_dim+1,2:y_dim+1)

g(N_vert+2,2:x_dim+1,2:y_dim+1)=-g(N_vert+1,2:x_dim+1,2:y_dim+1)


!solve for w,b,p,Phi,psi by applying the appropriate forcing operator to g
do i=1,5
	
	call MatMultForcing(var(2:N_vert+1,2:x_dim+1,2:y_dim+1,i),g,N_vert,x_dim,y_dim,hx,hy,hz,i)

enddo

!exchange boundary data for phi and psi
do i=4,5 !only need to exchange phi and psi boundaries
	call v_exchange(var(:,:,:,i),1)
enddo

!solve for horizontal velocities
call CalcVelocity(var(2:N_vert+1,:,:,4),var(2:N_vert+1,:,:,5),var(2:N_Vert+1,2:x_dim+1,2:y_dim+1,6),&
				  var(2:N_vert+1,2:x_dim+1,2:y_dim+1,7),N_vert,x_dim,y_dim)

!gather data to single write array
call get_full_array(var(2:N_vert+1,2:x_dim+1,2:y_dim+1,:))

!write data to file
call write_netcdf_vars(write_array(:,:,:,1),write_array(:,:,:,2),write_array(:,:,:,3),&
									write_array(:,:,:,4),write_array(:,:,:,5),write_array(:,:,:,6),write_array(:,:,:,7))

if(myid==0)then
	deallocate(write_array)
endif

if(myid==0)then
	call end_netcdf
endif

call mpi_end

end program main
