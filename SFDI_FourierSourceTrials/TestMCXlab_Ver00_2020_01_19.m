tic

clear;

cfg.nphoton=1e7;
cfg.maxdetphoton = 1e7; % save up to this many photons

%Dimension of domain
dim=60;
%Generate voxel positions
[xi,yi,zi]=meshgrid(1:dim,1:dim,1:dim);

cfg.vol=ones(size(xi));
cfg.vol(:,:,60) = 0; % pad a layer of 0s to get diffuse reflectance
cfg.vol = uint8(cfg.vol);

cfg.gpuid=1;
cfg.autopilot=1;

cfg.prop=[0 0 1 1;0.01 10 0.9 1.4];
cfg.seed=99999;

cfg.issaveexit=1; %Save the position and direction of detected photons (can't hurt)'
cfg.issaveref=1; %Also save diffuse reflectance

cfg.unitinmm = 1; %Each grid unit is 1 mm

cfg.srctype='fourier';
cfg.srcpos=[0 0 70];
cfg.srcdir=[0 0 -1];
fx = 0.05; % mm^-1
cfg.srcparam1=[60 0 0 fx*60];
cfg.srcparam2=[0 60 0 0];

cfg.tstart=0;
cfg.tend=1e-9;
cfg.tstep=1e-9;

cfg.isreflect=0; %Don''t reflect from outside the boundary (assume photons can escape)
cfg.isrefint=1; %Do consider internal reflections

[flux,detpt,vol,seeds,traj]=mcxlab(cfg);

fcw=flux.data*cfg.tstep;

figure;
hs=slice(log10(abs(double(fcw))),1,1,59);
set(hs,'linestyle','none');
axis equal; colorbar;box on;
title('a spatial frequency domain source 0.05 mm^{-1}');

figure;
imagesc(squeeze(log10(abs(double(fcw(:,30,:))))));
title('a spatial frequency domain source 0.05 mm^{-1}');

% % figure(101);
% % plot(squeeze(log10(abs(double(fcw(30,30,:))))));
% % title('a spatial frequency domain source');
% % fx = 0.3; % mm^-1
% % cfg.srcparam1=[60 0 0 fx*60];
% % flux=mcxlab(cfg);
% % fcw=flux.data*cfg.tstep;
% % hold all
% % plot(squeeze(log10(abs(double(fcw(30,30,:))))));
% % legend('0.05 mm^{-1}','0.3 mm^{-1}')

toc