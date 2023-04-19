
function test_profile_psi

close all;
%% Find x-point
% % Temperature
% n_rho1 = 5e-1;
% n_rho_max = 1e-1;
% n_rho_0 = 35e0;

% Density
n_rho1 = 5e18;
n_rho_max = 1e18;
n_rho_0 = 5e19;

% Scaling Parameters
a = 0.5;
b = 0.8;
rho_max = 1.3;
sol_rolloff = 0; % some number between [-50,0] say. 

    function n = myfun(rho)
        
        mask1 = rho<1;
        mask2 = rho>=1;

      
     
        N1 = numel(rho(:,1));
        N2= numel(rho(1,:));
        for i=1:N1
            for j=1:N2
        n1(i,j) = real(fun1(rho(i,j)));
        n2(i,j) = fun2(rho(i,j));
            end
        end
        
        n2 = n2.*mask2; 
        n2 = n2 .* ((1-sol_rolloff.*rho_max)+sol_rolloff .* rho); % make sinh steeper near the pedestal
        n2 = (n2-min(n2)) ./ (max(n2)-min(n2)) .* (n_rho1-n_rho_max) + n_rho1;
        
        n = n1.*mask1  + n2.*mask2;
        
% %          plot(rho,((1-sol_rolloff*rho_max)+sol_rolloff .* rho))
%         figure
%         plot(r,n1.*mask1)
%         figure
%         plot(r,n2)
% %         
    end

fun1 = @(rho) n_rho1+(n_rho_0-n_rho1).*((1-rho).^a).^b;
fun2 = @(rho) (-sinh((rho-1)/(rho_max-1)*2*pi-pi));
fun3 = @myfun;


% Read an g-eqdsk file
filename= 'g1050426022.01300';
g = readg_g3d(filename);

nR_n = 100;
nZ_n = 200;
offset = 3;
r=linspace(g.r(1+offset),g.r(end-offset),nR_n);
z=linspace(g.z(1+offset),g.z(end-offset),nZ_n);
[r2D,z2D] = meshgrid(r,z);
[psiN,psi] = calc_psiN(g,r2D(:),z2D(:),[]);

psi=reshape(psi,[nZ_n,nR_n]);
psiN=reshape(psiN,[nZ_n,nR_n]);
rho=sqrt(psiN);
xpt_info = find_xpt_jl(g,1,1,1e-8,1);
xz1 = xpt_info.zx;




ne = real(fun3(rho));
ne(r2D<2.4 & z2D<-0.65)=0.5;
ni=ne;
% Plot 2D density
figure; imagesc(r, z, ne);

%Plot rho w/ R
% figure; plot(r,rho(100,:))

% Plot density w/ R at midplane
% figure
% pcolor(r,z,real(fun1(rho)))
figure;
semilogy(r,real(fun3(rho(50,:))))
% xlim([0.8 1.2])
figure;
plot(r,real(fun3(rho(100,:))))
% xlim([0.8 1.2])

writematrix(ne,'ne.csv');
a=[r,z];
writematrix(a, 'axes.csv');

ncid = netcdf.create(('./density_west.nc'),'NC_WRITE');
dimR = netcdf.defDim(ncid,'nX_n',nR_n);

dimZ = netcdf.defDim(ncid,'nZ_n',nZ_n);

gridRnc = netcdf.defVar(ncid,'gridx_n','float',dimR);

gridZnc = netcdf.defVar(ncid,'gridz_n','float',dimZ);

nenc = netcdf.defVar(ncid,'ne','float',[dimR dimZ]);
ninc = netcdf.defVar(ncid,'ni','float',[dimR dimZ]);


netcdf.endDef(ncid);
% 
netcdf.putVar(ncid,gridRnc,r);
netcdf.putVar(ncid,gridZnc,z);


netcdf.putVar(ncid,nenc,ne);
netcdf.putVar(ncid,ninc,ni);


netcdf.close(ncid);






end
% 
%  x=linspace(2.8,3.3,100);
% y1=-sinh(x-2.98);
% figure,plot(x,y1)
% y2=(-sinh(x-2.98)/(max(x)-2.98));
% figure,plot(x,y2)
% 
% y3=(-sinh(x-1)/(3.3-3)+3.5*pi);
% figure,plot(x,y3)
% y4=(y3-min(y3))./(max(y3)-min(y3));
% figure,plot(x,y4)
% n_rho_max=1e18;
% n_rho1=5e18;
% y5=y4.*(n_rho1-n_rho_max)+n_rho1;
% figure,plot(x,y5)

