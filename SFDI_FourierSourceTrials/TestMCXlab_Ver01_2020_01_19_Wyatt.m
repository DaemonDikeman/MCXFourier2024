tic

clear all;

cfg.nphoton=1e7;
cfg.maxdetphoton = 1e7; % save up to this many photons

%Dimension of domain
dim=60;
%Generate voxel positions
[xi,yi,zi]=meshgrid(1:dim,1:dim,1:dim);

cfg.vol=ones(size(xi));
%cfg.vol(:,:,60) = 0; % pad a layer of 0s to get diffuse reflectance
cfg.vol = uint8(cfg.vol);

cfg.gpuid=1;
cfg.autopilot=1;

cfg.prop=[0 0 1 1;0.005 8 0.9 1.4];
cfg.seed=99999;

cfg.issaveexit=1; %Save the position and direction of detected photons (can't hurt)'
cfg.issaveref=1; %Also save diffuse reflectance

cfg.unitinmm = 1; %Each grid unit is 1 mm

cfg.srctype='fourier';
cfg.srcpos=[0 0 70];
cfg.srcdir=[0 0 -1];

cfg.tstart=0;
cfg.tend=1e-9;
cfg.tstep=1e-9;

cfg.isreflect=0; %Don''t reflect from outside the boundary (assume photons can escape)
cfg.isrefint=1; %Do consider internal reflections
cfg.outputtype = 'flux';


cfg1 = cfg;
cfg2 = cfg;
%%%%%%%%%%%%%%%%%%%%%
fx1 = 0.05; % mm^-1 %%
fx2 = 0.3; % mm^-1 %%
%%%%%%%%%%%%%%%%%%%%%
cfg1.srcparam1=[60 0 0 fx1*60];
cfg1.srcparam2=[0 60 0 0];

cfg2.srcparam1=[60 0 0 fx2*60];
cfg2.srcparam2=[0 60 0 0];

[flux1,detpt1,vol1,seeds1,traj1]=mcxlab(cfg1);
[flux2,detpt2,vol2,seeds2,traj2]=mcxlab(cfg2);

norm1 = flux1.stat.normalizer;
norm2 = flux2.stat.normalizer;

fcw1 = flux1.data*cfg.tstep/norm1;
fcw2 = flux2.data*cfg.tstep/norm2;
%%Plot entire volumes
figure;
hs=slice(abs(double(fcw1)),1,1,60);
caxis([0,1e-5])
set(hs,'linestyle','none');
axis equal; colorbar;box on;
title(['a spatial frequency domain source ' num2str(fx1) ' mm^{-1}']);

figure;
hs=slice(abs(double(fcw2)),1,1,60);
caxis([0,1e-5])
set(hs,'linestyle','none');
axis equal; colorbar;box on;
title(['a spatial frequency domain source ' num2str(fx2) ' mm^{-1}']);

%%Plot mean fluence(depth)
f1 = zeros(1,60);
f2 = f1;
for i = 1:60
    f1(i) = mean2(fcw1(:,:,60 - (i-1)));
    f2(i) = mean2(fcw2(:,:,60 - (i-1)));
end
figure
plot(f1);
hold on
plot(f2);
title('normalized fluence')
xlabel('depth (mm)')
legend(num2str(fx1),num2str(fx2))