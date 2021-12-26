% plot the upscale fluxes of all the variables as functions of the 
% damping coefficients d1 and d2

function all_fluxes_d1d2(d1_vec,d2_vec)

ncid = netcdf.open(strcat('grid_d',num2str(d1_vec(1)),'_d',num2str(d2_vec(1)),'_rad1.0_t0'),'NC_NOWRITE');

[dimname0,dimlen0] = netcdf.inqDim(ncid,0);
[dimname1,dimlen1] = netcdf.inqDim(ncid,1);
[dimname1,dimlen2] = netcdf.inqDim(ncid,2);

vecStrVariableList=['w','b','u','v'];

N=length(d1_vec);
N1=length(d2_vec);

matVar=zeros(N,N1,4,dimlen2,dimlen0,dimlen1);

dL=10.0;
Nx=dimlen0;
Ny=dimlen1;
hx=2*dL/Nx;
hy=2*dL/Ny;
hz=1.0/dimlen2;

zt=hz/2:hz:1-hz/2;

for i=1:3
    
    for j=1:N
        
        for k=1:N1
    
            matVar(j,k,i,:,:,:)=ncread(strcat('grid_d',num2str(d1_vec(j)),'_d',num2str(d2_vec(k)),'_rad1.0_t0'),vecStrVariableList(i));
            
        end
        
    end
    
end

matFlux=zeros(N,N1,4,dimlen2);

%compute products

for i=1:3
   
    for j=1:dimlen2
       
        for k=1:N
            
            for l=1:N1
            
                matFlux(k,l,i,j)=-sum(sum(matVar(k,l,1,j,:,:).*matVar(k,l,i,j,:,:)))*hx*hy;
                
            end
        
        end
        
    end
    
end

matDiffFlux=zeros(N,N1,4,dimlen2-1);

for i=1:3
    
    for k=1:N
        
        for l=1:N1
            
            matDiffFlux(k,l,i,:)=diff(matFlux(k,l,i,:))/hz;
            
        end
    
    end
    
end

%plot fluxes
for i=1:3
    figure(i)
    title(vecStrVariableList(i))
    for k=1:N
        subplot(2,2,k)
        for j=1:N1
            plot(squeeze(matDiffFlux(k,j,i,:)),zt(1:end-1),'LineWidth',1.5)
            %title(strcat('$(',vecStrVariableList(i),vecStrVariableList(1),')_z$'),'interpreter','latex')
            hold on
        end
        title(strcat('$d_{1}=$',num2str(d1_vec(k))),'interpreter','latex')
        legend('0.0625','0.25','1','4')
        hold off
    end
    %legend('0','10','45','90')
    %hold off
end

hold off

end
