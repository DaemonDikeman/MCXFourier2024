close all
clear
clear cfg

%% Single_Layer_Homogenous_Diffuse_Reflectance_Comparison
LUT_name = 'my_LUT';

%Scattering anisotropy
g = 0.71;
%Index of refraction
n = 1.33;

% L_star_dynamic = [0.5];
L_star_dynamic = [0.5,1,2,4];

freq_range = [0,0.1];
%[0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.10];
%[0, 0.05/3, 0.10/3, 0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5]

%Absorption range for dynamic layer
mua_dynamic = 1./(L_star_dynamic.*101);
%Scattering range for dynamic layer
musp_dynamic = mua_dynamic.*100;

% Number of Photons    %Make sure this value is equal to cfg.nphoton
cfg.nphoton= 1e7;      cfg.maxdetphoton = 1e7;

% Domain Settings
cfg.unitinmm = 1; %Side length of each voxel (mm)
cfg.vol = ones(100,100,40); %Size of entire domain (voxels)
cfg.vol(:,:,:)=2;
cfg.vol(:,:,1)=1;
% cfg.vol(:,:,1)=0;
% cfg.vol(:,:,2)=1;
% cfg.vol(:,:,3:40)=2;
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
gdev = gpuDevice(1);

%LUT allocation
% Mua = zeros(length(musp_dynamic),length(mua_dynamic));
% Musp = zeros(length(musp_dynamic),length(mua_dynamic));
% M1 = zeros(length(musp_dynamic),length(mua_dynamic));
% M2 = zeros(length(musp_dynamic),length(mua_dynamic));

%% Run Loop
%Initial Random Seed #
s = 0;
for i = 1:length(musp_dynamic) %length(musp_dynamic)
        
    
    % 0.067 0.64 0.82 1.4
    Musp(i,:) = musp_dynamic(i);
    mua_static = 0.067; 
    s = s + 1;
    cfg.seed = s;
    mua = mua_dynamic(i);
    mus = musp_dynamic(i)/(1-g);
    cfg.prop = [0 0 1 1; mua_static 0.64 0.82 1.4; mua mus g n]; %[mua mus g n]
    mcxpreview(cfg); title('domain preview')  
    for j = 1:length(freq_range)

        tic
        [~,dp]=mcxlab(cfg);
        toc
        % reset(gdev);
        
         good_photons = logical(dp.p(:,3)<0); % ???
%         gopho_x = gpuArray(dp.p(good_photons,1));
%         gopho_y = gpuArray(dp.p(good_photons,2));
%         gppath1 = gpuArray(dp.ppath(good_photons,1));
%         gppath2 = gpuArray(dp.ppath(good_photons,2));
%         rho_x = gopho_x*cfg.unitinmm;
%         rho_y = gopho_y*cfg.unitinmm;
%         ppath1 = cfg.unitinmm*gppath1; %Partial pathlength of first layer
%         ppath2 = cfg.unitinmm*gppath2;
        rho_x = cfg.unitinmm*dp.p(good_photons,1);
        rho_y = cfg.unitinmm*dp.p(good_photons,2);
        ppath1 = cfg.unitinmm*dp.ppath(good_photons,1); %Partial pathlength of first layer
        ppath2 = cfg.unitinmm*dp.ppath(good_photons,2); %Partial pathlength of first layer
        photon_displacement_x = rho_x - mean(rho_x); %mean(rho_x) is the location of the source (radial symmetry)
        photon_displacement_y = rho_y - mean(rho_y);
        photon_displacement = sqrt(photon_displacement_y.^2 + photon_displacement_x.^2);
        %drefmc(i,j) = mcxcwdref(dp(i,j), cfg);

        Mua(:,i) = mua;
%         M1(i,j) = mean(exp(-mua_static*double(ppath1) - mua*double(ppath2)).*cos(2*pi*freq_range(j)*double(photon_displacement)));
        M1a = (cos(double(photon_displacement_x)*2*pi*freq_range(j)).*cos(double(photon_displacement_y)*2*pi*freq_range(j)));
        M1b = exp(-1*mua_static*double(ppath1) - mua*double(ppath2));
        M1(i,j) = gather(sum(M1a.*M1b)/cfg.nphoton);
        %progressbar(j/length(mua_dynamic));

    end
    %progressbar(i/length(musp_dynamic));
end

%Export LUT
% LUT.M1 = M1;
% LUT.M2 = M2;
% LUT.Mua = Mua;
% LUT.Musp = Musp;
% if ~exist('./Generated_LUTs', 'dir')
%    mkdir('./Generated_LUTs')
% end
% save(['./Generated_LUTs/' LUT_name],'LUT');
filename = "SFDI_drefcalc_g_adjust_multilayer2_dual_axis_nphoton_" + num2str(cfg.nphoton, '%.1e') + "_voxel_" + num2str(cfg.unitinmm) + "mm.csv";
%writematrix(M1,filename);

%%
figure()
hold on
for i = 1:4
plot(freq_range, log10(M1(i,:)))
end
hold off
% ylim([-2 0])
% xlim([0 0.5])
% ylabel("Log10 Rd")
% xlabel("Spatial freq. (mm^{-1})")
% yticks([-2 -1.5 -1 -0.5 0])
% yticklabels(["10^{-2}," "10^{-1.5}", "10^{-1}", "10^{-0.5}", "10^{0}"])
% grid on
title(['Diffuse Reflectance vs. Spatial freq. (mm^{-1}) at diff. L* values'])

