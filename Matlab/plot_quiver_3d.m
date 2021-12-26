function plot_quiver_3d(theta_vec)

ncid = netcdf.open(strcat('grid_d0.25_d0.25_rad1.0_t',num2str(theta_vec(1))),'NC_NOWRITE');

[dimname0,dimlen0] = netcdf.inqDim(ncid,0);
[dimname1,dimlen1] = netcdf.inqDim(ncid,1);
[dimname1,dimlen2] = netcdf.inqDim(ncid,2);


vecStrVariableList=["w","b","u","v"];

dL=10.0; %length of the domain (ni each horizontal direction)
a=1.0; %radius of heating region
Nx=dimlen0;
Ny=dimlen1;
hx=2*dL/Nx;
hy=2*dL/Ny;
hz=1.0/dimlen2;
eps=0.01;
len_y=13;
len_x=13;
len_z=6;

ub=(Nx/2)+2*a/hx;
lb=(Nx/2)-2*a/hy;

xt=-dL+hx/2:hx:dL-hx/2;
yt=-dL+hy/2:hy:dL-hy/2;
zt=hz/2:hz:1-hz/2;

xt_plot=xt(lb:len_x:ub);
yt_plot=yt(lb:len_y:ub);
zt_plot=zt(1:len_z:end);

matVar=zeros(3,length(theta_vec),dimlen2,dimlen0,dimlen1);

for j=1:length(theta_vec)

    matVar(1,j,:,:,:)=ncread(strcat('grid_d0.25_d0.25_rad1.0_t',num2str(theta_vec(j))),...
                            convertStringsToChars(vecStrVariableList(1)));
    matVar(2,j,:,:,:)=ncread(strcat('grid_d0.25_d0.25_rad1.0_t',num2str(theta_vec(j))),...
                            convertStringsToChars(vecStrVariableList(3)));
    matVar(3,j,:,:,:)=ncread(strcat('grid_d0.25_d0.25_rad1.0_t',num2str(theta_vec(j))),...
                            convertStringsToChars(vecStrVariableList(4)));
        
end

[X,Z,Y]=meshgrid(xt_plot,zt_plot,yt_plot);
size(X)

size(squeeze(matVar(2,1,1:len_z:end,lb:len_x:ub,lb:len_y:ub)))

figure(1)

subplot(2,2,1)

quiver3(X,Y,Z,squeeze(matVar(3,1,1:len_z:end,lb:len_x:ub,lb:len_y:ub)),...
        squeeze(matVar(2,1,1:len_z:end,lb:len_x:ub,lb:len_y:ub)),...
        squeeze(matVar(1,1,1:len_z:end,lb:len_x:ub,lb:len_y:ub)),'LineWidth',1,'MaxHeadSize',1.85)
xlim([-2-4*eps 2+2*eps])
ylim([-2-4*eps 2+2*eps])

xlabel('x')
ylabel('z')
title(strcat('$\theta=$ \,',num2str(theta_vec(1))),'interpreter','latex')

subplot(2,2,2) 

quiver3(X,Y,Z,squeeze(matVar(3,2,1:len_z:end,lb:len_x:ub,lb:len_y:ub)),...
        squeeze(matVar(2,2,1:len_z:end,lb:len_x:ub,lb:len_y:ub)),...
        squeeze(matVar(1,2,1:len_z:end,lb:len_x:ub,lb:len_y:ub)),'LineWidth',1,'MaxHeadSize',1.85)
xlim([-2-4*eps 2+2*eps])
ylim([-2-4*eps 2+2*eps])

xlabel('x')
ylabel('z')
title(strcat('$\theta=$ \,',num2str(theta_vec(2))),'interpreter','latex')

subplot(2,2,3.5)

quiver3(X,Y,Z,squeeze(matVar(3,3,1:len_z:end,lb:len_x:ub,lb:len_y:ub)),...
        squeeze(matVar(2,3,1:len_z:end,lb:len_x:ub,lb:len_y:ub)),...
        squeeze(matVar(1,3,1:len_z:end,lb:len_x:ub,lb:len_y:ub)),'LineWidth',1,'MaxHeadSize',1.85)
xlim([-2-4*eps 2+2*eps])
ylim([-2-4*eps 2+2*eps])

xlabel('x')
ylabel('z')
title(strcat('$\theta=$ \,',num2str(theta_vec(3))),'interpreter','latex')
