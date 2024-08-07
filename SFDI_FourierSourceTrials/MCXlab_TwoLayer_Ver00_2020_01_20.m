tic

clear

cfg.nphoton=1e7; % number of photons
cfg.maxdetphoton = 1e7; % make sure this value is equal to cfg.nphoton

% domain settings
cfg.unitinmm = 1; % side length of each voxel (mm)
cfg.vol= 2*ones(60,60,60); % size of entire domain (voxels)
cfg.vol(:,:,1) = 1; % first layer (1 mm thick)
% % cfg.bc = 'aaaaar'; % boundary conditions (every side is absorbing except -z)
cfg.prop=[0 0 1 1; 0.01 10 0.9 1.4; 0.02 10 0.9 1.4]; % [mua mus g n]

% source settings
cfg.srctype='fourier';
cfg.srcdir=[0 0 -1];
cfg.srcpos=[0 0 70];
fx = 0.5; % cycles/mm
cfg.srcparam1=[60 0 0 fx*60];
cfg.srcparam2=[0 60 0 0];
% % cfg.issrcfrom0= 1; % 1 - first voxel is [0 0 0], 0 - first voxel is [1 1 1]

% detector settings
cfg.detpos = [30 30 0 120]; % ([x y z radius] of detector sphere)
cfg.savedetflag = 'px'; % save partial pathlengths and exit position

% time gate settings (seconds, ensure enough time for most of the photons to make it out of domain)
cfg.tstart=0;
cfg.tend=1e-6;
cfg.tstep=1e-6;

% processor info
cfg.gpuid= 1; % select which gpu driver to use
cfg.autopilot= 1; % 1-automatically set threads and blocks, 0-use nthread/nblocksize

% mcx simulation settings
cfg.seed=99999;
% % cfg.isnormalized = 1; % [1]-normalize the output fluence to unitary source, 0-no reflection

% run simulation
[flux,detp]=mcxlab(cfg);
fcw=flux.data*cfg.tstep;



% visualize results
figure;
hs=slice(log10(abs(double(fcw))),1,1,60);
set(hs,'linestyle','none');
axis equal; colorbar;box on;
% % caxis([-30 0])
title('2nd layer starts 1 mm, 0.5 mm^{-1}');

toc


