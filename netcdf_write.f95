module netcdf_write

use grid
use mpi_interface
use netcdf

implicit none

integer :: ncid
integer :: u_varid,v_varid,w_varid,phi_varid,pai_varid,b_varid,p_varid

contains

subroutine init_netcdf_vars

implicit none

integer,parameter :: ndims=3
integer :: u_dimids(ndims),v_dimids(ndims),w_dimids(ndims),b_dimids(ndims),xm_dimid,zm_dimid,xt_dimid,zt_dimid
integer :: phi_dimids(ndims),psi_dimids(ndims),p_dimids(ndims)
integer :: yt_dimid,ym_dimid,nc_rec,x_dim

nc_rec = nf90_create(file_name,NF90_CLOBBER,ncid)
!nc_rec = nf90_def_dim(ncid,"xm",N_x,xm_dimid)
!nc_rec = nf90_def_dim(ncid,"ym",N_y,ym_dimid)
!nc_rec = nf90_def_dim(ncid,"zm",N_vert,zm_dimid)
nc_rec = nf90_def_dim(ncid,"xt",N_x,xt_dimid)
nc_rec = nf90_def_dim(ncid,"yt",N_y,yt_dimid)
nc_rec = nf90_def_dim(ncid,"zt",N_vert,zt_dimid)
!nc_rec = nf90_def_dim(ncid,"time",nf90_unlimited,rec_dimid)

w_dimids=(/zt_dimid,xt_dimid,yt_dimid/) 
b_dimids=(/zt_dimid,xt_dimid,yt_dimid/)
p_dimids=(/zt_dimid,xt_dimid,yt_dimid/)
psi_dimids=(/zt_dimid,xt_dimid,yt_dimid/)
phi_dimids=(/zt_dimid,xt_dimid,yt_dimid/)
u_dimids=(/zt_dimid,xt_dimid,yt_dimid/)
v_dimids=(/zt_dimid,xt_dimid,yt_dimid/)


nc_rec = nf90_def_var(ncid,"w",nf90_double,w_dimids,w_varid)
nc_rec = nf90_def_var(ncid,"b",nf90_double,b_dimids,b_varid)
nc_rec = nf90_def_var(ncid,"p",nf90_double,p_dimids,p_varid)
nc_rec = nf90_def_var(ncid,"phi",nf90_double,phi_dimids,phi_varid)
nc_rec = nf90_def_var(ncid,"psi",nf90_double,psi_dimids,pai_varid)
nc_rec = nf90_def_var(ncid,"u",nf90_double,u_dimids,u_varid)
nc_rec = nf90_def_var(ncid,"v",nf90_double,v_dimids,v_varid)

nc_rec = nf90_enddef(ncid)

end subroutine init_netcdf_vars


subroutine write_netcdf_vars(W,B,P,Phi,Psi,U,V)

implicit none

integer :: nc_rec,varid
real*8 :: time
real*8, dimension(:,:,:) :: U,V,W,B,Phi,Psi,P

nc_rec = nf90_inq_varid(ncid,"w",varid)
nc_rec= nf90_put_var(ncid,varid,W)
nc_rec = nf90_inq_varid(ncid,"b",varid)
nc_rec = nf90_put_var(ncid,varid,B)
nc_rec = nf90_inq_varid(ncid,"p",varid)
nc_rec = nf90_put_var(ncid,varid,P)
nc_rec = nf90_inq_varid(ncid,"phi",varid)
nc_rec = nf90_put_var(ncid,varid,Phi)
nc_rec = nf90_inq_varid(ncid,"psi",varid)
nc_rec = nf90_put_var(ncid,varid,Psi)
nc_rec = nf90_inq_varid(ncid,"u",varid)
nc_rec= nf90_put_var(ncid,varid,U)
nc_rec = nf90_inq_varid(ncid,"v",varid)
nc_rec= nf90_put_var(ncid,varid,V)
!nc_rec = nf90_close(ncid)


end subroutine write_netcdf_vars

subroutine end_netcdf

implicit none

integer :: nc_rec

nc_rec = nf90_close(ncid)

end subroutine end_netcdf

end module
