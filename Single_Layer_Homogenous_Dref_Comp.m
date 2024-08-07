%% Single_Layer_Homogenous_Diffuse_Reflectance_Comparison
clear
close all

LUT_name = 'my_LUT';

%Scattering anisotropy
g = 0.71; % alt base: 0.71 OR 0.9
%Index of refraction
n = 1.33; % alt base: 1.33 OR 1.4

L_star_dynamic = [0.5,1,2,4]; %[0.5,1,2,4]

freqfoc = 0; % 1 -> frequencies from 0 to 0.10; 0 -> frequencies from 0 to 0.50

if freqfoc == 1
    freq_range = [0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.10];
else
    freq_range = [0, 0.05/3, 0.10/3, 0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5];
end

axnum = 2; % 1 -> single axis; 2 -> dual axis; 3 -> omnidirectional

%Absorption range for dynamic layer
mua_dynamic = 1./(L_star_dynamic.*101);
%Scattering range for dynamic layer
musp_dynamic = mua_dynamic.*100;

% Number of Photons    %Make sure this value is equal to cfg.nphoton
cfg.nphoton= 1e7;      cfg.maxdetphoton = 1e7;

% Domain Settings
cfg.unitinmm = 1; %Side length of each voxel (mm)
cfg.vol = ones(100,100,40); %Size of entire domain (voxels)
cfg.vol(:,:,1) = 0;
cfg.bc = 'aaraaa'; %Boundary Conditions (every side is absorbing except -z)

% Source Settings
cfg.srctype='pencil';
cfg.srcpos=[50 50 1]; %(voxels)
direction = [0 0 1]; % vector
cfg.issrcfrom0=1;
cfg.srcdir=direction/norm(direction);

% Detector Settings
cfg.detpos = [50 50 0 80]; %([x y z radius] of detector sphere)
cfg.savedetflag = 'dpx'; % save partial pathlengths and exit position

% Time Gate Settings (ensure enough time for most of the photons to make it out of domain)
cfg.tstart=0;
cfg.tend=1e-9;
cfg.tstep=1e-9;

% Processor Info
cfg.gpuid=1; %Select which GPU driver to use
cfg.autopilot=1;

%LUT allocation
Mua = zeros(length(musp_dynamic),length(mua_dynamic));
Musp = zeros(length(musp_dynamic),length(mua_dynamic));
M1 = zeros(length(musp_dynamic),length(mua_dynamic));
% M2 = zeros(length(musp_dynamic),length(mua_dynamic));

tTotal = zeros(length(musp_dynamic),1);
tAvg = zeros(length(musp_dynamic),1);

%% Run Loop
%Initial Random Seed #
s = 0;
for i = 1:length(musp_dynamic) %length(musp_dynamic)
    
    Musp(i,:) = musp_dynamic(i);
    s = s + 1;
    cfg.seed = s;
    mua = mua_dynamic(i);
    mus = musp_dynamic(i)/(1-g);
    %cfg.prop = [0 0 1 1; 0.067, 0.64, 0.82, 1.4; mua mus g n]; %[mua mus g n]
    cfg.prop = [0 0 1 1; mua mus g n]; %[mua mus g n]
    T = zeros(length(freq_range),1);
    tStart = tic;
    for j = 1:length(freq_range)

        tRun = tic;
        %[~,dp(i,j)]=mcxlab(cfg);
        [~,dp]=mcxlab(cfg);
        
        good_photons = logical(dp.p(:,3)<1); % ???
        Mua(:,i) = mua_dynamic(i);
        
        if axnum == 1
        rho_x = cfg.unitinmm*dp.p(good_photons,1);
        ppath = cfg.unitinmm*dp.ppath(good_photons); %Partial pathlength of first layer
        photon_displacement_x = rho_x - mean(rho_x); %mean(rho_x) is the location of the source (radial symmetry)
        M1(i,j) = sum(exp(-mua_dynamic(i)*double(ppath)).*(cos(2*pi*freq_range(j)*double(photon_displacement_x))))/cfg.nphoton;
        
        elseif axnum == 2
        rho_x = cfg.unitinmm*dp.p(good_photons,1);
        rho_y = cfg.unitinmm*dp.p(good_photons,2);
        ppath = cfg.unitinmm*dp.ppath(good_photons); %Partial pathlength of first layer
        photon_displacement_x = rho_x - mean(rho_x); %mean(rho_x) is the location of the source (radial symmetry)
        photon_displacement_y = rho_y - mean(rho_y);

        M1(i,j) = sum(exp(-mua_dynamic(i)*double(ppath)).*(cos(2*pi*freq_range(j)*double(photon_displacement_x)).*cos(2*pi*freq_range(j)*double(photon_displacement_y))))/cfg.nphoton;
        
        elseif axnum == 3 % omnidirectional
        rho_x = cfg.unitinmm*dp.p(good_photons,1);
        rho_y = cfg.unitinmm*dp.p(good_photons,2);
        ppath = cfg.unitinmm*dp.ppath(good_photons); %Partial pathlength of first layer
        photon_displacement_x = rho_x - mean(rho_x); %mean(rho_x) is the location of the source (radial symmetry)
        photon_displacement_y = rho_y - mean(rho_y);
        photon_displacement = sqrt(photon_displacement_x.^2 + photon_displacement_y.^2);
        % drefmc(i,j) = mcxcwdref(dp(i,j), cfg);

        M1(i,j) = sum(exp(-mua_dynamic(i)*double(ppath)).*(cos(2*pi*freq_range(j)*double(photon_displacement))))/cfg.nphoton;
        
        end
        
        progressbar(j/length(freq_range));
        T(j) = toc(tRun);
        Tmessage = ['Run ', num2str(j), ' Time elapsed: ', num2str(T(j)), ' seconds'];
        disp(Tmessage);

    end
    tTotal(i) = toc(tStart);
    tAvg(i) = sum(T)/length(freq_range);
    Tmessage = ['Average Run Time: ',num2str(tAvg(i)),' seconds'];
    disp(Tmessage);
    progressbar(i/length(musp_dynamic));
end

%% Export LUT
if axnum == 1
    axstring = "_single_axis";
elseif axnum == 2
    axstring = "_dual_axis";
elseif axnum == 3
    axstring = "_omni_axis";
else
    axstring= "_unknown_axis";
end

if freqfoc == 1
    axstring = axstring + "_focused_freq";
end

filename = "SFDI_drefcalc\\" + "SFDI_drefcalc_g_" + num2str(g,'%.2f') + axstring + "_nphoton_" + num2str(cfg.nphoton, '%.1e') + "_vx" + num2str(cfg.unitinmm) + "mm.csv";
% writematrix(M1,filename);

%%

figure()
hold on
for i = 1:4
plot(freq_range, log10(M1(i,:)))
end
hold off
% % ylim([-2 0])
% % xlim([0 0.5])
% % ylabel("Log10 Rd")
% % xlabel("Spatial freq. (mm^{-1})")
% % yticks([-2 -1.5 -1 -0.5 0])
% % yticklabels(["10^{-2}," "10^{-1.5}", "10^{-1}", "10^{-0.5}", "10^{0}"])
% % grid on
title(['Diffuse Reflectance vs. Spatial freq. (mm^{-1}) at diff. L* values'])

mcxpreview(cfg); title('domain preview')

