function plot_variables_yz(filename)

ncid = netcdf.open(filename,'NC_NOWRITE');

[dimname0,dimlen0] = netcdf.inqDim(ncid,0);
[dimname1,dimlen1] = netcdf.inqDim(ncid,1);
[dimname1,dimlen2] = netcdf.inqDim(ncid,2);

vecStrVariableList=["w","b","u","v"];

matVar=zeros(4,dimlen2,dimlen0,dimlen1);

dL=10.0;
a=1.0;
Nx=dimlen0;
Ny=dimlen1;
hx=2*dL/Nx;
hy=2*dL/Ny;
hz=1.0/dimlen2;

y_slice=dimlen1/2;

ub=(Nx/2)+2*a/hx
lb=(Nx/2)-2*a/hy

xt=-dL+hx/2:hx:dL-hx/2;
yt=-dL+hy/2:hy:dL-hy/2;
zt=hz/2:hz:1-hz/2;

matVar(1,:,:,:)=ncread(filename,'w');
matVar(2,:,:,:)=ncread(filename,'b');
matVar(3,:,:,:)=ncread(filename,'u');
matVar(4,:,:,:)=ncread(filename,'v');

%get plot bounds for u,v,w

bound_vec=zeros(1,4);

for i=1:4
    maxVar=max(max(squeeze(matVar(i,:,y_slice,:))))
    minVar=min(min(squeeze(matVar(i,:,y_slice,:))))
    bound_vec(i)=min([abs(maxVar),abs(minVar)])
    
end

bound_vec

h_vec=zeros(1,4);

for i=1:4
   
    h_vec(i)=2*bound_vec(i)/10;
    
end

figure(1)
for i=1:4
    
    subplot(2,2,i)
    %if(i==2)
     %   contourf(yt(lb:ub),zt,(squeeze(matVar(i,:,dimlen0/2,lb:ub))),10)

    %else
        contourf(yt(lb:ub),zt,(squeeze(matVar(i,:,dimlen0/2,lb:ub))),-bound_vec(i):h_vec(i):bound_vec(i))
    %end
    colormap(redblueTecplot)
    colorbar
    title(vecStrVariableList(i),'interpreter','latex')
    xlabel('y')
    ylabel('z')
    hold on
    
end

hold off

end
