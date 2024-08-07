%% Two_Layer_LUT_Generation_MCX
addpath('mcxlab','Generated_LUTs');
LUT_name = 'my_LUTb';

%Scattering anisotropy
g = 0.9;
%Index of refraction
n = 1.4;

static_scalar = 1;

%Absorption for static layer
mua_static = 0.096 * static_scalar; % in mm^-1
%Scattering for static layer
musp_static = 0.78; % in mm^-1
mus_static = static_scalar * (musp_static/(1-g));

%Absorption range for dynamic layer
mua_dynamic = linspace(0,0.2,400);
%Scattering range for dynamic layer
musp_dynamic = linspace(0.01,3,200);

%Spatial frequencies (mm^-1)
fx1 = 0;
fx2 = 0.1;

% Number of Photons    %Make sure this value is equal to cfg.nphoton
cfg.nphoton= 1e7;      cfg.maxdetphoton = 1e7;

% Domain Settings
xlen = 100;
ylen = 100;
zlen = 100;
cfg.unitinmm = 1; %Side length of each voxel (mm)
cfg.vol = 2*ones(xlen,ylen,50); %Size of entire domain (voxels)
cfg.vol(:,:,1) = 1; %static first layer
cfg.bc = 'aaraaa'; %Boundary Conditions (every side is absorbing except -z)

% Source Settings
cfg.srctype='pencil';
cfg.srcpos=[xlen/2 ylen/2 0]; %(voxels)
direction = [0 0 1]; % vector
cfg.issrcfrom0=1;
cfg.srcdir=direction/norm(direction);

% Detector Settings
cfg.detpos = [xlen/2 ylen/2 0 (max([xlen ylen])*1.5)]; %([x y z radius] of detector sphere)
cfg.savedetflag = 'px'; % save partial pathlengths and exit position

% Time Gate Settings (ensure enough time for most of the photons to make it out of domain)
cfg.tstart=0;
cfg.tend=1e-6;
cfg.tstep=1e-6;

% Processor Info
cfg.gpuid=1; %Select which GPU driver to use
cfg.autopilot=1;

%LUT allocation
Mua = zeros(length(musp_dynamic),length(mua_dynamic));
Musp = zeros(length(musp_dynamic),length(mua_dynamic));
M1 = zeros(length(musp_dynamic),length(mua_dynamic));
M2 = zeros(length(musp_dynamic),length(mua_dynamic));

%Initial Random Seed #
s = 0;
for i = 1:length(musp_dynamic)
    Musp(i,:) = musp_dynamic(i);
    s = s + 1;
    cfg.seed = s;
    mus = musp_dynamic(i)/(1-g);
    cfg.prop = [0 0 1 1; 0 mus_static g 1.4; 0 mus g 1.4]; %[mua mus g n]
    tic
    [~,dp]=mcxlab(cfg);
    toc
    
    good_photons = logical(dp.p(:,3)<0); % ???
    gopho_x = gpuArray(dp.p(good_photons,1));
    gopho_y = gpuArray(dp.p(good_photons,1));
    gppath = gpuArray(dp.ppath(good_photons));
    rho_x = gopho_x*cfg.unitinmm;
    rho_y = gopho_y*cfg.unitinmm;
    ppath = cfg.unitinmm*gppath; %Partial pathlength of first layer
    photon_displacement_x = rho_x - mean(rho_x); %mean(rho_x) is the location of the source (radial symmetry)
    photon_displacement_y = rho_y - mean(rho_y);
    
    good_photons = logical(dp.p(:,3)<0);
    rho_x = cfg.unitinmm*dp.p(good_photons,1)';
    ppath1 = cfg.unitinmm*dp.ppath(good_photons,1)'; %Partial pathlength of first layer
    ppath2 = cfg.unitinmm*dp.ppath(good_photons,2)'; %Partial pathlength of second (dynamic) layer
    photon_displacement = rho_x - mean(rho_x); %mean(rho_x) is the location of the source (radial symmetry)
    for j = 1:length(mua_dynamic)
        Mua(:,j) = mua_dynamic(j);
        M1(i,j) = mean(exp(-mua_static*double(ppath1) - mua_dynamic(j)*double(ppath2)).*cos(2*pi*fx1*double(photon_displacement)));
        M2(i,j) = mean(exp(-mua_static*double(ppath1) - mua_dynamic(j)*double(ppath2)).*cos(2*pi*fx2*double(photon_displacement)));
        %progressbar(j/length(mua_dynamic));
    end
    %progressbar(i/length(musp_dynamic));
end

%Export LUT
LUT.M1 = M1;
LUT.M2 = M2;
LUT.Mua = Mua;
LUT.Musp = Musp;
save(['Generated_LUTs\' LUT_name],'LUT');