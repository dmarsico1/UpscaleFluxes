%plot the upscale fluxes as functions of latitude

function all_fluxes_theta(theta_vec)

ncid = netcdf.open(strcat('grid_d0.25_d0.25_rad1.0_t22.5'),'NC_NOWRITE');

[dimname0,dimlen0] = netcdf.inqDim(ncid,0);
[dimname1,dimlen1] = netcdf.inqDim(ncid,1);
[dimname1,dimlen2] = netcdf.inqDim(ncid,2);

vecStrVariableList=['w','b','u','v'];

N=length(theta_vec);

matVar=zeros(N,4,dimlen2,dimlen0,dimlen1);

dL=10.0;
Nx=dimlen0;
Ny=dimlen1;
hx=2*dL/Nx;
hy=2*dL/Ny;
hz=1.0/dimlen2;

zt=hz/2:hz:1-hz/2;

for i=1:4
    
    for j=1:N
    
        matVar(j,i,:,:,:)=ncread(strcat('grid_d0.25_d0.25_rad1.0_t',num2str(theta_vec(j))),vecStrVariableList(i));
        
    end
    
end

matFlux=zeros(N,4,dimlen2);

%compute products

for i=1:4
   
    for j=1:dimlen2
       
        for k=1:N
            
           matFlux(k,i,j)=-sum(sum(matVar(k,1,j,:,:).*matVar(k,i,j,:,:)))*hx*hy;
        
        end
        
    end
    
end

matDiffFlux=zeros(N,4,dimlen2-1);

for i=1:4
    
    for k=1:N
        
        matDiffFlux(k,i,:)=diff(matFlux(k,i,:))/hz;
    
    end
    
end

%compute fluxes
for i=1:4
    subplot(2,2,i)
    for k=1:N
        plot(squeeze(matDiffFlux(k,i,:)),zt(1:end-1),'LineWidth',1.5)
        title(strcat('$-(',vecStrVariableList(i),vecStrVariableList(1),')_z$'),'interpreter','latex')
        zlabel('z')
        hold on
    end
    legend('22.5','45','67.5')
    hold off
end

figure(2)

for i=1:4
    for k=1:N
        plot(squeeze(matDiffFlux(k,i,:)),zt(1:end-1),'LineWidth',1.5)
        title('Upscale Flux of Vertical Momentum','interpreter','latex','FontSize',13)
        ylabel('z','FontSize',12)
        hold on
    end
    legend({'22.5','45','67.5'},'FontSize',13)
    hold off
end

figure(3)

i=1;
subplot(1,2,1)
for k=1:N
    plot(squeeze(matDiffFlux(k,i,:)),zt(1:end-1),'LineWidth',1.5)
    title("Upscale Flux of"+newline+" Vertical Momentum",'interpreter','latex','FontSize',13)
    zlabel('z')
    hold on
end
legend('22.5','45','67.5')
hold off

i=3; 
subplot(1,2,2)
for k=1:N
    plot(squeeze(matDiffFlux(k,i,:)),zt(1:end-1),'LineWidth',1.5)
    title("Upscale Flux of"+newline+" Zonal Momentum",'interpreter','latex','FontSize',13)
    zlabel('z')
    hold on
end
legend('22.5','45','67.5')
hold off

hold off


end
