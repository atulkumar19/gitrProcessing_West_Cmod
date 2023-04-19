function test_profile

close all;
clear all;
clc;

n_rho1 = 5e18;
n_rho_max = 1e18;
n_rho_0 = 5e19;
a = 0.5;
b = 0.8;
rho_max = 3.2;
sol_rolloff = 0; % some number between [-50,0] say. 

    function n = myfun(rho)
        
        mask1 = rho<1;
        mask2 = rho>=1;
        
        N = numel(rho);
        for i=1:N
            n1(i) = fun1(rho(i));
            n2(i) = fun2(rho(i));
        end
        
        n2 = n2 .* mask2;
        n2 = n2 .* ((1-sol_rolloff.*rho_max)+sol_rolloff .* rho); % make sinh steeper near the pedestal
        n2 = (n2-min(n2)) ./ (max(n2)-min(n2)) .* (n_rho1-n_rho_max) + n_rho1;
        
        n = n1 .* mask1 + n2 .* mask2;
        
         plot(rho,((1-sol_rolloff*rho_max)+sol_rolloff .* rho))
        figure
        plot(rho,n1.*mask1)
        figure
        plot(rho,n2.*mask2)
        
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
size(rho)

xpt_info = find_xpt_jl(g,1,1,1e-8,1);
xz1 = xpt_info.zx;

% rho=linspace(0.8,1.3,100);
figure
semilogy(rho,fun3(rho),'Marker','s');
figure
plot(rho,fun3(rho),'Marker','s');

end


