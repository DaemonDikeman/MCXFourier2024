% a regular pencil beam at the center of the volume
clear cfg;
cfg.nphoton=1e7;
cfg.vol=uint8(ones(60,60,60));
cfg.srcpos=[30 30 -10];
cfg.srcdir=[0 0 1];
cfg.gpuid=1;
cfg.autopilot=1;
cfg.prop=[0 0 1 1;0.005 1 0.8 1.37];
cfg.tstart=0;
cfg.tend=1e-9;
cfg.tstep=1e-9;
%cfg.printnum=10;
cfg.seed=99999;
cfg.srctype='isotropic';
flux=mcxlab(cfg);
fcw=flux.data*cfg.tstep;

figure;
hs=slice(log10(abs(double(fcw))),[],30,1);
set(hs,'linestyle','none');
axis equal; colorbar
title('isotrpoic source at [30 30 -10]');