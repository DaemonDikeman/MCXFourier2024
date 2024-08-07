%==========================================================================
%  Large detector  pencil beam in semi infinite medium.  This should be the
%  simplest MC case that I can compare with an analytical solution in order
%  to test diffuse reflectance calculations
%==========================================================================
clear cfg

%number of photons to launch
cfg.nphoton=1e7; 
%Dimension of domain
dim=60;
%Generate voxel positions
[xi,yi,zi]=meshgrid(1:dim,1:dim,1:dim);
%dist=sqrt((xi-30).^2+(zi-30).^2);

%Set the volume 
cfg.vol=ones(size(xi));
cfg.vol(:,:,1)=0;
cfg.vol(:,:,5:60) = 2;
cfg.vol=uint8(cfg.vol);

%Tell it what gpu to use
cfg.gpuid='1'; %GPU ID
cfg.autopilot=1; %Automatically choose how many threads, etc.

cfg.srctype='pencil';  %pencil beam
cfg.srcpos=[30 30 1]; %Location of pencil beam
cfg.srcdir=[0 0 1]; %Headed in the +z direction

cfg.maxdetphoton=1e7; %save up to this many photons
cfg.issaveexit=1; %Save the position and direction of detected photons (can't hurt)'
cfg.issaveref=1; %Also save diffuse reflectance

cfg.unitinmm=0.5; %Each grid unit is 1 mm

%%%Set up timing (mostly for visualizing flux)
cfg.tstart=0;   %Time the simulation starts
cfg.tstep=5e-12; %Steps to take
cfg.tend = 5e-9; %When to end.  The output will have [tstart:tstep:tend] slices at each of the different time points

%Photons are only detected at interfaces!
cfg.detpos = [30,40,1,2;
              30,52, 1, 2]; %detector at center, 2 mm layer
          
cfg.isreflect=0; %Don''t reflect from outside the boundary (assume photons can escape)
cfg.isrefint=1; %Do consider internal reflections


cfg.prop=[0 0 1 1;  %Properties of the materials.  anywhere cfg.vol == 0 will have properties in the first row (medium type 0)
          0.01,10,.9,1.4;
          1,20,.9,1.4]; %mu_a mu_s (!!!not mus prime!!!), g, n
     

 [flux,det]=mcxlab(cfg);
 save(fullfile('/project/botlab/Matt','intro2MCX_seminfinite'),'flux','det','cfg','-v7.3');


%Let's plot some fluxes
%NB: Flux measurements depend on mua ~= 0
h=figure;
set(h,'Position',[36,672,1208,425])
subplot(141)
imagesc(squeeze(flux.data(30,:,:,5))')
title('t=25 ps')
subplot(142)
imagesc(squeeze(flux.data(30,:,:,50))')
title('t=250 ps')
subplot(143)
imagesc(squeeze(flux.data(30,:,:,100))')
title('t=500 ps')
subplot(144)
imagesc(squeeze(flux.data(30,:,:,150))')
title('t=750 ps')

%You can also plot the positions of detected photons
det1Idx = find(det.detid == 1);
det2Idx = find(det.detid == 2);
h=figure;
subplot(131)
histogram(det.p(det2Idx,1))
title('det photons (X)')
subplot(132)
histogram(det.p(det2Idx,2))
title('det photons (Y)')
subplot(133)
histogram(det.p(det2Idx,3))
title('det photons (Z)')

%Using the partial pathlengths in each medium you can get the temporal
%point spread function
n=1.4; %index of refraction
c_mmps = 3e11/n; %Speed of light in medium
secPerSamp = 1e-12; %How many bins to split the tpsf into
fs = 1/secPerSamp; %effective sampling rate

edges = 0:secPerSamp:10e-9; %Edges of the histogram
%time steps
tSteps = linspace(secPerSamp/2,max(edges)-secPerSamp/2, length(edges)-1);
%Distance steps
pathSteps = tSteps * c_mmps;
%Frequency bins
freqBins = ([0:length(tSteps)-1] * fs / length(tSteps))/1e6;
%Transit times
times = double(det.ppath(det1Idx))/c_mmps;
%Calculate tpsf
tpsf = histcounts(times,edges);

%Plot TPSF
figure
plot(tSteps/1e-9,tpsf)
ylabel('Detected Photons')
xlabel('Time (ns)')
title('Temporal point spread function')

%Translate tpsf to frequency response
paddedTPSF = zeros(2^15,1); %Pad TPSF for smoothness
fResp = zeros(length(tSteps),1);
paddedTPSF(1:length(tSteps)) = tpsf;

paddedFResp = fft(paddedTPSF);
%Padded frequency bins
paddedBins = ([0:(length(paddedFResp)-1)] * fs / length(paddedFResp))/1e6;

figure
subplot(121)
plot(paddedBins,abs(paddedFResp))
xlim([50,2000])
subplot(122)
plot(paddedBins,-angle(paddedFResp))
xlim([50,2000])









%%%%Same thing but with two layers
clear cfg
clear cfgVec
%number of photons to launch
cfg.nphoton=1e7; 
%Dimension of domain
dim=60;
%Generate voxel positions
[xi,yi,zi]=meshgrid(1:dim,1:dim,1:dim);
%dist=sqrt((xi-30).^2+(zi-30).^2);

%Set the volume 
cfg.vol=ones(size(xi));
cfg.vol(:,:,1)=0;
cfg.vol(:,50:end)= 2;
cfg.vol=uint8(cfg.vol);

figure
imagesc(squeeze(cfg.vol(30,:,:)))
%Tell it what gpu to use
cfg.gpuid='1'; %GPU ID
cfg.autopilot=1; %Automatically choose how many threads, etc.

cfg.srctype='pencil';  %pencil beam
cfg.srcpos=[30 30 1]; %Location of pencil beam
cfg.srcdir=[0 0 1]; %Headed in the +z direction

cfg.maxdetphoton=1e7; %save up to this many photons
cfg.issaveexit=1; %Save the position and direction of detected photons (can't hurt)'
cfg.issaveref=1; %Also save diffuse reflectance

cfg.unitinmm=1; %Each grid unit is 1 mm

%%%Set up timing (mostly for visualizing flux)
cfg.tstart=0;   %Time the simulation starts
cfg.tstep=5e-12; %Steps to take
cfg.tend = 5e-9; %When to end.  The output will have [tstart:tstep:tend] slices at each of the different time points

%Photons are only detected at interfaces!
cfg.detpos = [30,40,1,2;
              30,52, 1, 2]; %detector at center, 2 mm layer
          
cfg.isreflect=0; %Don''t reflect from outside the boundary (assume photons can escape)
cfg.isrefint=1; %Do consider internal reflections
cfg.prop=[0 0 1 1;  %Properties of the materials.  anywhere cfg.vol == 0 will have properties in the first row (medium type 0)
          0.01,10,.9,1.4]; %mu_a mu_s (!!!not mus prime!!!), g, n
     

 [flux,det]=mcxlab(cfg);
 save(fullfile('/project/botlab/Matt','intro2MCX_seminfinite'),'flux','det','cfg','-v7.3');


%Let's plot some fluxes
%NB: Flux measurements depend on mua ~= 0
h=figure;
set(h,'Position',[36,672,1208,425])
subplot(141)
imagesc(squeeze(flux.data(30,:,:,5))')
title('t=5 ps')
subplot(142)
imagesc(squeeze(flux.data(30,:,:,50))')
title('t=50 ps')
subplot(143)
imagesc(squeeze(flux.data(30,:,:,150))')
title('t=150 ps')
subplot(144)
imagesc(squeeze(flux.data(30,:,:,250))')
title('t=250 ps')

%You can also plot the positions of detected photons
det1Idx = find(det.detid == 1);
det2Idx = find(det.detid == 2);
h=figure;
subplot(131)
histogram(det.p(det1Idx,1))
title('det photons (X)')
subplot(132)
histogram(det.p(det1Idx,2))
title('det photons (Y)')
subplot(133)
histogram(det.p(det1Idx,3))
title('det photons (Z)')
