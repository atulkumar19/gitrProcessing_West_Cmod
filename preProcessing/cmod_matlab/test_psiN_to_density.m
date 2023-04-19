clearvars;

gfile_name = 'g1050426022.01300';
g = readg_g3d(gfile_name);

a=1;
b=1;

%% Find x-point
xpt_info = find_xpt_jl(g,1,1,1e-8,1);
xz1 = xpt_info.zx;

%% Set up rectangular evaluation mesh
nr = 201;
nz = 201;
r = linspace(0.4,1.1,nr);
z = linspace(-0.6,0.6,nz);

%% Build some cells and cell centers
icount = 0;
for i = 1:length(r)-1
    for j = 1:length(z)-1
        icount = icount + 1;
        rCell(:,icount) = [r(i),r(i+1),r(i+1),r(i)];
        zCell(:,icount) = [z(j),z(j),z(j+1),z(j+1)];
        rCenter(j,i) = (r(i)+r(i+1))/2;
        zCenter(j,i) = (z(j)+z(j+1))/2;
    end
end

psiNCenter = reshape(calc_psiN(g,rCenter(:),zCenter(:)),nz-1,nr-1);


%% Put a good n = f(psiN) here!
density = 5e18 + (5e19 - 5e18).*((max(psiNCenter(:))-psiNCenter).^a).^b;

%% Set a constant PFR density (could add a function here)
density(zCenter < xz1 & psiNCenter < 1) = 1e18;

%% Some mask on R
density(rCenter > 1 )  = 10e19;

figure; hold on; box on; grid on; set(gcf,'color','w'); set(gca,'fontsize',14)
h = patch(rCell,zCell,density(:),'edgecolor','none');
contour(g.r,g.z,[(g.psirz-g.ssimag)/(g.ssibry-g.ssimag)].',[1,1],'k-','linewidth',2);
plot(g.lim(1,:),g.lim(2,:),'b-','linewidth',2)
colorbar;
axis equal;

figure; 
plot(r, [density(:, 100);0])

    