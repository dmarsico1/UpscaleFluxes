function fluxes(filename)

ncid = netcdf.open(filename,'NC_NOWRITE');

[dimname0,dimlen0] = netcdf.inqDim(ncid,0);
[dimname1,dimlen1] = netcdf.inqDim(ncid,1);
[dimname1,dimlen2] = netcdf.inqDim(ncid,2);

vecStrVariableList=['w','b','u','v'];

matVar=zeros(4,dimlen2,dimlen0,dimlen1);

dL=10.0;
Nx=dimlen0;
Ny=dimlen1;
hx=2*dL/Nx;
hy=2*dL/Ny;
hz=1.0/dimlen2;

zt=hz/2:hz:1-hz/2;

for i=1:4
    
    matVar(i,:,:,:)=ncread(filename,vecStrVariableList(i));
    
end

matFlux=zeros(4,dimlen2);

%compute products

for i=1:4
   
    for j=1:dimlen2
       
        matFlux(i,j)=-sum(sum(matVar(1,j,:,:).*matVar(i,j,:,:)))*hx*hy;
        
    end
    
end

matDiffFlux=zeros(4,dimlen2-1);

for i=1:4
    
   matDiffFlux(i,:)=diff(matFlux(i,:))/hz;
    
end

%compute fluxes
for i=1:4
    
    subplot(2,2,i)
    plot(matDiffFlux(i,:),zt(1:end-1))
    title(strcat('$(',vecStrVariableList(i),vecStrVariableList(1),')_z$'),'interpreter','latex')
    hold on
    
end

hold off


end
