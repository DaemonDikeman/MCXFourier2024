tic

clear;

cfg.nphoton=1e7;
cfg.maxdetphoton = 1e7; % save up to this many photons

%Dimension of domain
dim=100;
%Generate voxel positions
[xi,yi,zi]=meshgrid(1:dim,1:dim,1:dim);

% % cfg.vol=ones(size(xi));
cfg.vol=ones(size(xi));
cfg.vol(:,:,1:70) = 2;
%cfg.vol(:,:,60) = 0; % pad a layer of 0s to get diffuse reflectance
cfg.vol = uint8(cfg.vol);

cfg.gpuid=1;
cfg.autopilot=1;

mua = 0.01; % mm^-1
mus = 10; % mm^-1
g = 0.9; % anisotropy
n = 1.4; % refractive index

muavol = [0; mua; 2*mua];
musvol = [0; mus; mus];
gvol = [1; g; g];
nvol = [1; n; n];

cfg.prop=[muavol musvol gvol nvol];

cfg.seed=99999;
cfg.issaveexit=1; %Save the position and direction of detected photons (can't hurt)'
cfg.issaveref=1; %Also save diffuse reflectance

cfg.unitinmm = 1/10; %Each grid unit in mm

cfg.srctype='fourier';
cfg.srcpos=[0 0 dim+10];
cfg.srcdir=[0 0 -1];

cfg.tstart=0;
cfg.tend=1e-9;
cfg.tstep=1e-9;

cfg.isreflect=0; %Don''t reflect from outside the boundary (assume photons can escape)
cfg.isrefint=1; %Do consider internal reflections
cfg.outputtype = 'flux';


cfg1 = cfg;
cfg2 = cfg;
cfg3 = cfg;
%%%%%%%%%%%%%%%%%%%%%
fx1 = 0; % mm^-1 %%
fx2 = 0.01; % mm^-1 %%
fx3 = 0.3; % mm^-1 %%
%%%%%%%%%%%%%%%%%%%%%
cfg1.srcparam1=[dim 0 0 fx1*dim];
cfg1.srcparam2=[0 dim 0 0];

cfg2.srcparam1=[dim 0 0 fx2*dim];
cfg2.srcparam2=[0 dim 0 0];

cfg3.srcparam1=[dim 0 0 fx3*dim];
cfg3.srcparam2=[0 dim 0 0];

[flux1,detpt1,vol1,seeds1,traj1]=mcxlab(cfg1);
[flux2,detpt2,vol2,seeds2,traj2]=mcxlab(cfg2);
[flux3,detpt3,vol3,seeds3,traj3]=mcxlab(cfg3);

norm1 = flux1.stat.normalizer;
norm2 = flux2.stat.normalizer;
norm3 = flux3.stat.normalizer;

% fluence (1/mm^2 or J/mm^2) = integral(flux*dt) = sum(flux.data,4), where
% 4th dimension is time
fcw1 = flux1.data*cfg.tstep/norm1;
fcw2 = flux2.data*cfg.tstep/norm2;
fcw3 = flux3.data*cfg.tstep/norm3;

%%Plot entire volumes
figure;
% % hs=slice(abs(double(fcw1)),1,1,dim);
% % caxis([0,5e-6])
hs=slice(log10(abs(double(fcw1))),1,1,dim);
caxis([-7 -5.5])
set(hs,'linestyle','none');
axis equal; colorbar;box on;
title(['Fluence SFDI f_x = ' num2str(fx1) ' mm^{-1}']);

figure;
% % hs=slice(abs(double(fcw2)),1,1,dim);
% % caxis([0,5e-6])
hs=slice(log10(abs(double(fcw2))),1,1,dim);
caxis([-7 -5.5])
set(hs,'linestyle','none');
axis equal; colorbar;box on;
title(['Fluence SFDI f_x = ' num2str(fx2) ' mm^{-1}']);

figure;
% % hs=slice(abs(double(fcw2)),1,1,dim);
% % caxis([0,5e-6])
hs=slice(log10(abs(double(fcw3))),1,1,dim);
caxis([-7 -5.5])
set(hs,'linestyle','none');
axis equal; colorbar;box on;
title(['Fluence SFDI f_x = ' num2str(fx3) ' mm^{-1}']);

%%Plot mean fluence(depth)
f1 = zeros(1,dim);
f2 = f1;
f3 = f1;
for ii = 1:dim
    f1(ii) = mean2(fcw1(:,:,dim - (ii-1)));
    f2(ii) = mean2(fcw2(:,:,dim - (ii-1)));
    f3(ii) = mean2(fcw3(:,:,dim - (ii-1)));
end
figure;
depth = 0:cfg.unitinmm:(dim*cfg.unitinmm-cfg.unitinmm);
plot(depth,f1);
% % plot(log10(f1));
hold on
plot(depth,f2);
% % plot(log10(f2));
plot(depth,f3);
title('normalized fluence')
xlabel('depth (mm)')
legend(num2str(fx1),num2str(fx2),num2str(fx3))

%% Energy deposited manually calculated (method 1)

% % % energy deposited (unitless or J) = fluence*mua*voxel_volume
% % Edep1 = flux1.data*cfg.tstep*mua*(cfg.unitinmm^3)/norm1;
% % Edep2 = flux2.data*cfg.tstep*mua*(cfg.unitinmm^3)/norm2;
% % 
% % %%Plot entire volumes
% % figure;
% % hs=slice(log10(abs(double(Edep1))),1,1,dim);
% % caxis([-12 -10.5])
% % set(hs,'linestyle','none');
% % axis equal; colorbar;box on;
% % title(['Energy deposited SFDI f_x = ' num2str(fx1) ' mm^{-1}']);
% % 
% % figure;
% % hs=slice(log10(abs(double(Edep2))),1,1,dim);
% % caxis([-12 -10.5])
% % set(hs,'linestyle','none');
% % axis equal; colorbar;box on;
% % title(['Energy deposited SFDI f_x = ' num2str(fx2) ' mm^{-1}']);
% % 
% % %%Plot mean energy(depth)
% % E1 = zeros(1,dim);
% % E2 = E1;
% % for i = 1:dim
% %     E1(i) = mean2(Edep1(:,:,dim - (i-1)));
% %     E2(i) = mean2(Edep2(:,:,dim - (i-1)));
% % end
% % figure;
% % depth = 0:cfg.unitinmm:(dim*cfg.unitinmm-cfg.unitinmm);
% % plot(depth,E1);
% % % % plot(log10(f1));
% % hold on
% % plot(depth,E2);
% % % % plot(log10(f2));
% % title('energy deposited manually calculated')
% % xlabel('depth (mm)')
% % legend(num2str(fx1),num2str(fx2))


%% Energy deposited auto calculated by mcx (method 2)

% energy deposited (unitless or J) = fluence*mua*voxel_volume
cfg1.outputtype = 'energy';
cfg2.outputtype = 'energy';
cfg3.outputtype = 'energy';
energy1 = mcxlab(cfg1);
energy2 = mcxlab(cfg2);
energy3 = mcxlab(cfg3);
norm1 = energy1.stat.normalizer;
norm2 = energy2.stat.normalizer;
norm3 = energy3.stat.normalizer;
Edep1 = energy1.data/norm1;
Edep2 = energy2.data/norm2;
Edep3 = energy3.data/norm3;
% % Edep1 = energy1.data*cfg.tstep*mua/norm1;
% % Edep2 = energy2.data*cfg.tstep*mua/norm2;

%%Plot entire volumes
figure;
% % hs=slice(log10(abs(double(Edep1))),1,1,dim);
% % caxis([-12 -10.5])
hs=slice(abs(double(Edep1)),1,1,dim);
caxis([0 2.5])
set(hs,'linestyle','none');
axis equal; colorbar;box on;
title(['Energy deposited SFDI f_x = ' num2str(fx1) ' mm^{-1}']);

figure;
% % hs=slice(log10(abs(double(Edep2))),1,1,dim);
% % caxis([-12 -10.5])
hs=slice(abs(double(Edep2)),1,1,dim);
caxis([0 2.5])
set(hs,'linestyle','none');
axis equal; colorbar;box on;
title(['Energy deposited SFDI f_x = ' num2str(fx2) ' mm^{-1}']);

figure;
% % hs=slice(log10(abs(double(Edep2))),1,1,dim);
% % caxis([-12 -10.5])
hs=slice(abs(double(Edep3)),1,1,dim);
caxis([0 2.5])
set(hs,'linestyle','none');
axis equal; colorbar;box on;
title(['Energy deposited SFDI f_x = ' num2str(fx3) ' mm^{-1}']);

%%Plot mean energy(depth)
E1 = zeros(1,dim);
E2 = E1;
E3 = E1;
for ii = 1:dim
    E1(ii) = mean2(Edep1(:,:,dim - (ii-1)));
    E2(ii) = mean2(Edep2(:,:,dim - (ii-1)));
    E3(ii) = mean2(Edep3(:,:,dim - (ii-1)));
end
figure;
depth = 0:cfg.unitinmm:(dim*cfg.unitinmm-cfg.unitinmm);
plot(depth,E1);
% % plot(log10(f1));
hold on
plot(depth,E2);
% % plot(log10(f2));
plot(depth,E3);
title('energy deposited auto calculated by mcx')
xlabel('depth (mm)')
legend(num2str(fx1),num2str(fx2),num2str(fx3))


toc
