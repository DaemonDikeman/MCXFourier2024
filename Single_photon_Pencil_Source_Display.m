%% Single_Layer_Homogenous_Diffuse_Reflectance_Comparison
LUT_name = 'my_LUT';
close all
%Scattering anisotropy
g = 0.9;
%Index of refraction
n = 1.4;

L_star_dynamic = [1];

%Absorption range for dynamic layer
mua_dynamic = 1./(L_star_dynamic.*101);
%Scattering range for dynamic layer
musp_dynamic = mua_dynamic.*100;

% Number of Photons    %Make sure this value is equal to cfg.nphoton
cfg.nphoton= 1;      cfg.maxdetphoton = 1e3;

% Domain Settings
cfg.unitinmm = 1; %Side length of each voxel (mm)
cfg.vol = ones(40,40,20); %Size of entire domain (voxels)
cfg.bc = 'aaraaa'; %Boundary Conditions (every side is absorbing except -z)

% Source Settings
cfg.srctype='pencil';
cfg.srcpos=[20 20 0]; %(voxels)
direction = [0 0 1]; % vector
cfg.issrcfrom0=1;
cfg.srcdir=direction/norm(direction);

% Detector Settings
cfg.detpos = [20 20 0 30]; %([x y z radius] of detector sphere)
cfg.savedetflag = 'dpx'; % save partial pathlengths and exit position

% Time Gate Settings (ensure enough time for most of the photons to make it out of domain)
cfg.tstart=0;
cfg.tend=1e-9;
cfg.tstep=1e-9;

% Processor Info
cfg.gpuid=1; %Select which GPU driver to use
cfg.autopilot=1;

%LUT allocation
% Mua = zeros(length(musp_dynamic),length(mua_dynamic));
% Musp = zeros(length(musp_dynamic),length(mua_dynamic));
% M1 = zeros(length(musp_dynamic),length(mua_dynamic));
% M2 = zeros(length(musp_dynamic),length(mua_dynamic));

%% Run Loop
%Initial Random Seed #
s = 0;
for i = 1:length(musp_dynamic) %length(musp_dynamic)
        
    Musp(i,:) = musp_dynamic(i);
    s = s + 1;
    cfg.seed = s;
    mua = mua_dynamic(i);
    mus = musp_dynamic(i)/(1-g);
    cfg.prop = [0 0 1 1; mua mus g n]; %[mua mus g n]
        
    j = 1;
    
    tic
    [flux, dp(i,j), vol, seeds]=mcxlab(cfg);
    toc

    good_photons = logical(dp(i,j).p(:,3)<0); % ???
    rho_x = cfg.unitinmm*dp(i,j).p(good_photons,1)';
    ppath = cfg.unitinmm*dp(i,j).ppath(good_photons)'; %Partial pathlength of first layer
    photon_displacement = rho_x - mean(rho_x); %mean(rho_x) is the location of the source (radial symmetry)
    %drefmc(i,j) = mcxcwdref(dp(i,j), cfg);

    Mua(:,i) = mua_dynamic(i);
    %progressbar(j/length(mua_dynamic));

    %progressbar(i/length(musp_dynamic));
end

%%

newcfg=cfg;
newcfg.seed=seeds.data;
newcfg.outputtype='jacobian';
newcfg.detphotons=dp(i,j);
[flux2, detp2, vol2, seeds2, trajectory]=mcxlab(newcfg);
jac=sum(flux2.data,4);
imagesc(log10(abs(squeeze(jac(:,30,:)))))

newcfg.outputtype='rf';
newcfg.omega=2*pi*100e6; % 100 MHz RF modulation
newcfg.detphotons=dp(i,j);
rfjac=mcxlab(newcfg);
jac=sum(rfjac.data,4);
figure;imagesc(log10(abs(squeeze(jac(:,30,:)))))

ax = figure(3);
p = plot3(trajectory.pos(:,1), trajectory.pos(:,2), trajectory.pos(:,3));
p.LineStyle = "-";
grid on
%p.Marker = ".";
xlim([0 40])
ylim([0 40])
zlim([0 20])
camOrbit('fps', 5, 'time', 15);

figure(4)
p2 = plot3(trajectory.pos(:,1), trajectory.pos(:,2), trajectory.pos(:,3));
p2.LineStyle = "-";
p2.Marker = ".";
grid on


%camOrbit(p2, 'fps', 10, 'time', 10);

% FINNSTAR7 (2022). camOrbit (https://www.mathworks.com/matlabcentral/fileexchange/74975-camorbit), MATLAB Central File Exchange. Retrieved March 24, 2022.