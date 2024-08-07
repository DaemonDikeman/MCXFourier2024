%% Single_Layer_Homogeneous_LUT_Generation
%addpath('mcxlab','mcxlab\fangq-iso2mesh-14d64c7');
LUT_name = 'my_LUT';

%Scattering anisotropy
g = 0.71;
%Index of refraction
n = 1.33;

%Absorption range for dynamic layer
mua_dynamic = linspace(0,0.2,400); %(0,0.2, 400)
%Scattering range for dynamic layer
musp_dynamic = linspace(0.01,3,200); %(0.01,3, 200)

%Spatial frequencies (mm^-1)
fx1 = 0;
fx2 = 0.1;

% Number of Photons    %Make sure this value is equal to cfg.nphoton
cfg.nphoton= 1e7;      cfg.maxdetphoton = 1e7;

% Domain Settings
cfg.unitinmm = 10; %Side length of each voxel (mm)
cfg.vol = ones(20,20,20); %Size of entire domain (voxels)
cfg.bc = 'aaraaa'; %Boundary Conditions (every side is absorbing except -z)

% Source Settings
cfg.srctype='pencil';
cfg.srcpos=[10 10 0]; %(voxels)
direction = [0 0 1]; % vector
cfg.issrcfrom0=1;
cfg.srcdir=direction/norm(direction);

% Detector Settings
cfg.detpos = [10 10 0 30]; %([x y z radius] of detector sphere)
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
for i = 100:100%length(musp_dynamic)
    Musp(i,:) = musp_dynamic(i);
    s = s + 1;
    cfg.seed = s;
    mus = musp_dynamic(i)/(1-g);
    cfg.prop = [0 0 1 1; 0 mus g 1.4]; %[mua mus g n]
    tic
    [~,dp]=mcxlab(cfg);
    toc
    good_photons = logical(dp.p(:,3)<0);
    rho_x = cfg.unitinmm*dp.p(good_photons,1)';
    rho_y = cfg.unitinmm*dp.p(good_photons,2)';
    
    ppath = cfg.unitinmm*dp.ppath(good_photons)'; %Partial pathlength of first layer
    photon_displacement_x = rho_x - mean(rho_x); %mean(rho_x) is the location of the source (radial symmetry)
    photon_displacement_y = rho_y - mean(rho_y);
    photon_displacement = sqrt(photon_displacement_x.^2 + photon_displacement_y.^2);
    for j = 1:length(mua_dynamic)
        Mua(:,j) = mua_dynamic(j);
        M1(i,j) = mean(exp(-mua_dynamic(j)*double(ppath)).*cos(2*pi*fx1*double(photon_displacement)));
        M2(i,j) = mean(exp(-mua_dynamic(j)*double(ppath)).*cos(2*pi*fx2*double(photon_displacement)));
        %progressbar(j/length(mua_dynamic));
    end
    %progressbar(i/length(musp_dynamic));
end
mcxpreview(cfg); title('domain preview')
% %Export LUT
% LUT.M1 = M1;
% LUT.M2 = M2;
% LUT.Mua = Mua;
% LUT.Musp = Musp;
% if ~exist('./Generated_LUTs', 'dir')
%    mkdir('./Generated_LUTs')
% end
% save(['./Generated_LUTs/' LUT_name],'LUT');
