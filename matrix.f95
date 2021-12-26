module matrix

use grid
use mpi_interface

contains

subroutine CalcVelocity(phi,psi,u,v,M,N,L)

implicit none

real*8, intent(in), dimension(M,N+2,L+2) :: phi,psi
real*8, intent(out), dimension(M,N,L) :: u,v
integer, intent(in) :: M,N,L
real*8, dimension(M,N,L) :: phix,phiy,psix,psiy

call Grad(phi,phix,phiy,M,N,L)

call Grad(psi,psix,psiy,M,N,L)

u=-phix-psiy

v=-phiy+psix

end subroutine CalcVelocity

subroutine Grad(p,px,py,M,N,L)

real*8, dimension(M,N+2,L+2), intent(in) :: p
real*8, dimension(M,N,L), intent(out) :: px,py
integer, intent(in) :: M,N,L
integer :: i

do i=2,N+1
	px(:,i-1,:)=(0.5/hx)*(p(:,i+1,2:L+1)-p(:,i-1,2:L+1))
enddo

do i=2,L+1
	py(:,:,i-1)=(0.5/hy)*(p(:,2:N+1,i+1)-p(:,2:N+1,i-1))
enddo

end subroutine Grad

subroutine MatMultForcing(p,q,M,N,L,hx,hy,hz,opt)

implicit none
real*8, dimension(M+2,N+2,L+2), intent(in) :: q !this is what gets INPUT.  Layer of ghost cells around q.
real*8, dimension(M,N,L), intent(out) :: p !this is what gets OUTPUT.  No ghost cells around p.
real*8, intent(in) :: hx,hz,hy
integer, intent(in) :: M,N,L,opt
integer :: i,j,k
real*8 :: pxx,pyy,pzz,pyz,theta,alpha,beta,px,py,pz,p1y,p2y,p3y

alpha=dbeta/dalpha!d2/d1
beta=dalpha*dbeta


if(opt==1)then !w equation
	do j=2,N+1
		do i=2,M+1
			do k=2,L+1
				pxx=(1.0/hx**2)*(q(i,j+1,k)-2.0*q(i,j,k)+q(i,j-1,k))
				pyy=(1.0/hy**2)*(q(i,j,k+1)-2.0*q(i,j,k)+q(i,j,k-1))
				pzz=(1.0/hz**2)*(q(i+1,j,k)-2.0*q(i,j,k)+q(i-1,j,k))
				p(i-1,j-1,k-1)=(pxx+pyy)/(1.0+beta)
			enddo
		enddo
	enddo

elseif(opt==2)then !b equation
	do j=2,N+1
		do i=2,M+1
			do k=2,L+1
				pxx=(1.0/hx**2)*(q(i,j+1,k)-2.0*q(i,j,k)+q(i,j-1,k))
				pyy=(1.0/hy**2)*(q(i,j,k+1)-2.0*q(i,j,k)+q(i,j,k-1))
				if(i==2)then
					pzz=(1.0/hz**2)*(2.0*q(i,j,k)-5.0*q(i+1,j,k)+4.0*q(i+2,j,k)-q(i+3,j,k))
					p1y=(0.5/hy)*(q(i,j,k+1)-q(i,j,k-1))
					p2y=(0.5/hy)*(q(i+1,j,k+1)-q(i+1,j,k-1))
					p3y=(0.5/hy)*(q(i+2,j,k+1)-q(i+2,j,k-1))
					pyz=(1.0/hz)*(-1.5*p1y+2.0*p2y-0.5*p3y)
				elseif(i==M+1)then
					pzz=(1.0/hz**2)*(-q(i-3,j,k)+4.0*q(i-2,j,k)-5.0*q(i-1,j,k)+2.0*q(i,j,k))
					p1y=(0.5/hy)*(q(i,j,k+1)-q(i,j,k-1))
					p2y=(0.5/hy)*(q(i-1,j,k+1)-q(i-1,j,k-1))
					p3y=(0.5/hy)*(q(i-2,j,k+1)-q(i-2,j,k-1))
					pyz=(1.0/hz)*(0.5*p3y-2.0*p2y+1.5*p1y)
				else
					pzz=(1.0/hz**2)*(q(i+1,j,k)-2.0*q(i,j,k)+q(i-1,j,k))
					pyz=(1.0/(4.0*hy*hz))*(q(i+1,j,k+1)-q(i+1,j,k-1)-q(i-1,j,k+1)+q(i-1,j,k-1))
				endif
				theta=(4.0*atan(1.0)/180.0)*dtheta
				p(i-1,j-1,k-1)=(dalpha/(1.0+beta))*(pxx+pyy+pzz) &
								+(1.0/(dalpha*(1.0+beta)))*((cos(theta)**2)*pyy+2.0*cos(theta)*sin(theta)*pyz+(sin(theta)**2)*pzz)
			enddo
		enddo
	enddo

elseif(opt==3)then !p equation
	do j=2,N+1
		do i=2,M+1
			do k=2,L+1
				px=(q(i,j+1,k)-q(i,j-1,k))/(2.0*hx)
				py=(q(i,j,k+1)-q(i,j,k-1))/(2.0*hy)
				pz=(q(i+1,j,k)-q(i-1,j,k))/(2.0*hz)
				theta=(4.0*atan(1.0)/180.0)*dtheta
				p(i-1,j-1,k-1)=(1.0/(1.0+beta))*((((sin(theta)**2)/(dalpha))+dalpha)*pz + (1.0/dalpha)*sin(theta)*cos(theta)*py - cos(theta)*px)
								
			enddo
		enddo
	enddo

elseif(opt==4)then !phi equation
	do j=2,N+1
		do i=2,M+1
			do k=2,L+1			
				px=(q(i,j+1,k)-q(i,j-1,k))/(2.0*hx)
				py=(q(i,j,k+1)-q(i,j,k-1))/(2.0*hy)
				pz=(q(i+1,j,k)-q(i-1,j,k))/(2.0*hz)
				theta=(4.0*atan(1.0)/180.0)*dtheta
				p(i-1,j-1,k-1)=(1.0/(1.0+beta))*pz		
			enddo
		enddo
	enddo

elseif(opt==5)then !psi equation
	do j=2,N+1
		do i=2,M+1
			do k=2,L+1
				px=(q(i,j+1,k)-q(i,j-1,k))/(2.0*hx)
				py=(q(i,j,k+1)-q(i,j,k-1))/(2.0*hy)
				pz=(q(i+1,j,k)-q(i-1,j,k))/(2.0*hz)
				theta=(4.0*atan(1.0)/180.0)*dtheta
				p(i-1,j-1,k-1)=(1.0/(dalpha*(1.0+beta)))*(sin(theta)*pz+cos(theta)*py)		
			enddo
		enddo
	enddo
	
	
elseif(opt==6)then !u equation

	do j=2,N+1
		do i=2,M+1
			do k=2,L+1
				px=(q(i,j+1,k)-q(i,j-1,k))/(2.0*hx)
				py=(q(i,j,k+1)-q(i,j,k-1))/(2.0*hy)
				pz=(q(i+1,j,k)-q(i-1,j,k))/(2.0*hz)
				theta=(4.0*atan(1.0)/180.0)*dtheta
				p(i-1,j-1,k-1)=(1.0/(dalpha*(1.0+beta)))*(sin(theta)*pz+cos(theta)*py)		
			enddo
		enddo
	enddo


endif

end subroutine MatMultForcing



subroutine MatMult(p,q,M,N,L,hx,hy,hz)

implicit none
real*8, dimension(M+2,N+2,L+2), intent(in) :: q !this is what gets INPUT.  Layer of ghost cells around q.
real*8, dimension(M,N,L), intent(out) :: p !this is what gets OUTPUT.  No ghost cells around p.
real*8, intent(in) :: hx,hz,hy
integer, intent(in) :: M,N,L
integer :: i,j,k
real*8 :: pxx,pyy,pzz,pyz,theta,alpha,beta

alpha=dbeta/dalpha
beta=dalpha*dbeta


do j=2,N+1
	do i=2,M+1
		do k=2,L+1
			pxx=(1.0/hx**2)*(q(i,j+1,k)-2.0*q(i,j,k)+q(i,j-1,k))
			pyy=(1.0/hy**2)*(q(i,j,k+1)-2.0*q(i,j,k)+q(i,j,k-1))
			pzz=(1.0/hz**2)*(q(i+1,j,k)-2.0*q(i,j,k)+q(i-1,j,k))
			pyz=(1.0/(4.0*hy*hz))*(q(i+1,j,k+1)-q(i+1,j,k-1)-q(i-1,j,k+1)+q(i-1,j,k-1))
			theta=(4.0*atan(1.0)/180.0)*dtheta
			p(i-1,j-1,k-1)=pxx+pyy &
							+(dalpha/(1.0+beta))*((cos(theta)**2)*pyy+2.0*cos(theta)*sin(theta)*pyz+(sin(theta)**2)*pzz) &
							+(beta/(1.0+beta))*pzz
		enddo
	enddo
enddo

end subroutine MatMult

subroutine ConjGrad(x,f,tol,N,M,L,hx,hy,hz)

implicit none
real*8, dimension(M+2,N+2,L+2), intent(out) :: x !This is what gets output
real*8, dimension(M,N,L), intent(in) :: f !This is the right hand side
real*8, intent(in) :: tol
integer :: i,j,k
real*8, dimension(M,N,L) :: r,q,Ax
integer, intent(in) :: N,M,L
real*8 :: hx,hz,hy

!x(2:M+1,2:N+1,2:L+1) = f

r=f

rhs_sum=sum(r)

call mpi_allreduce(rhs_sum,rhs_sum_tot,1,mpi_real8,mpi_sum,comm2d,ierr)

!r=r-(rhs_sum_tot/(N_x*N_vert*N_y))

p(2:M+1,2:N+1,2:L+1)=r

rold=sum(r*r)

call mpi_allreduce(rold,rold_tot,1,mpi_real8,mpi_sum,comm2d,ierr)


do k=1,N_x**2
	call conj_exchange(p,k)
	p(1,:,:)=-p(2,:,:)
	p(M+2,:,:)=-p(M+1,:,:)
	!write(*,*) size(p,1),size(p,2),size(p,3)
	call matmult(q,p,M,N,L,hx,hy,hz)
	gam=sum(q*p(2:M+1,2:N+1,2:L+1))
	call mpi_allreduce(gam,gam_tot,1,mpi_real8,mpi_sum,comm2d,ierr)
	alpha_tot=rold_tot/gam_tot
	x(2:M+1,2:N+1,2:L+1) = x(2:M+1,2:N+1,2:L+1)+alpha_tot*p(2:M+1,2:N+1,2:L+1)
	r=r-alpha_tot*q
	rnew=sum(r*r)
	call mpi_allreduce(rnew,rnew_tot,1,mpi_real8,mpi_sum,comm2d,ierr)
	if(sqrt(rnew_tot)<tol)then
		exit
	endif
	p(2:M+1,2:N+1,2:L+1) = r + (rnew_tot/rold_tot)*p(2:M+1,2:N+1,2:L+1)
	rold_tot=rnew_tot
	write(*,*) k,rnew_tot
enddo

write(*,*) k, rnew
	
end subroutine ConjGrad

end module matrix
