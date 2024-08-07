%% SFDI 2-layer white MLMC 
% Karissa Tilbury & Christian Crane 07_27_2021
% Daemon Dikeman 11_17_2021
close all % close all windows opened by previous scripts
clear
clear cfg
rng('default')

%% CONTROL PANEL 
%info=mcxlab('gpuinfo');
%Total Volume Creation (mm)
x                = 140;     % Total x dimension of box
y                = 100;     % Total y dimension of box
z                = 20;      % Total z dimension of box
lstar            = 1;     % 1/(mua+musp) ; constant ratio of musp/mua = 100 used
cfg.margin       = 20;      % Margin between edge of box and edge of pattern source
cfg.unitinmm     = 1;       % length of a voxel in mm can't be larger than layer thickness estimation (ex. 5.1mm thickness would require voxel size 0.1 or smaller)
cfg.nphoton      = 1e7;     % number of photons
cfg.srctype = 'fourier';    % type of pattern/source to project; fourier pattern is our interest
freq = [0, 0.10]; % set of spatial frequencies to run in mm^-1
% [0, 0.05/3, 0.10/3, 0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5]
layer2mult = 2;
layer3mult = 0.5;
%cfg.minenergy = 0.0001;
zerolayerflag = 1;          % determines whether we use a layer of zero-value voxels to save diffuse reflectance; 1 is yes
trajflag = 0;               % determines whether to calculate and display trajectories or not.
displayflag = 0;            % determines whether or not to display results appropriate to the run type
multilayerflag = 1;
SACflag = 0;

if multilayerflag == 1
    filename = "SFDI_multilayer_drefval_2layermodelBCthin_nphoton_" + num2str(cfg.nphoton, '%.1e') + "_voxel_" + num2str(cfg.unitinmm) + "mm.csv";
end 
if multilayerflag == 0
    filename = "SFDI_LUT_g_0.8_homogenous_Lstar_" + num2str(lstar) + "_nphoton_" + num2str(cfg.nphoton, '%.1e') + "_voxel_" + num2str(cfg.unitinmm) + "mm.csv";
end

%% Layer Creation
numLayer = 3; % number of tissue layers (not including top air layer (layer 1))
layerthickness(2)=1; %thickness of epidermis in mm (layer 1 is air)
layerthickness(3)=2; %thickness of dermis in mm

mua_dynamic = linspace(0,0.2,100); % what absorbance values do we want to iterate through?

musp_dynamic = linspace(0.01,5,100); % what scattering values do we want to iterate through?

n = 1.4; % index of refraction
g = 0.82; % anisotropy
% constant ratio of musp/mua = 100
% Lstar = 1/(musp+mua)
mua=1/(lstar*101); % from above eq.: Lstar = 1/(101*mua)
musp=100*mua; % from above eq.
mus=musp/(1-g); % mu_s = (mu_s_prime) / (1-anisotropy)

%% Properties from Skin Anatomy Calculations Table
% [mua, mus, g, n]
epiderm_mus     = 6.413;
epiderm_mualit  = 0.66766;
epiderm_muadrk  = 3.26464;
epiderm_g       = 0.82;
epiderm_n       = 1.4;
epiderm_thick   = 0.10;

epimus_work     = epiderm_mus*(epiderm_thick/ceil(epiderm_thick));
epimua_worklit  = epiderm_mualit*(epiderm_thick/ceil(epiderm_thick));
epimua_workdrk  = epiderm_muadrk*(epiderm_thick/ceil(epiderm_thick));
epiderm_laylit  = [epimua_worklit, epimus_work, epiderm_g, epiderm_n];
epiderm_laydrk  = [epimua_workdrk, epimus_work, epiderm_g, epiderm_n];

dermis_mus      = 2.627;
dermis_mua      = 0.04292;
dermis_g        = 0.82;
dermis_n        = 1.4;
dermis_thick    = 1.8;

dermmus_work    = dermis_mus*(dermis_thick/ceil(dermis_thick));
dermmua_work    = dermis_mua*(dermis_thick/ceil(dermis_thick));
dermis_layer    = [dermmua_work, dermmus_work, dermis_g, dermis_n];

subcut_mus      = 0.500;
subcut_mua      = 0.03273;
subcut_g        = 0.75;
subcut_n        = 1.44;

subcut_layer    = [subcut_mua, subcut_mus, subcut_g, subcut_n];

%% CONSTANTS DO NOT CHANGE=================================================
cfg.gpuid           = 1;
cfg.autopilot       = 1;
cfg.isreflect       = 0;       % Don't reflect from outside the boundary
cfg.isrefint        = 1;       % Do consider internal reflections
cfg.maxdetphoton    = 1e7;     % save up to this many photons (effects memory significantly) !DO NOT INCREASE ABOVE 1E^8!
cfg.issaveexit      = 0;       % Save the position and direction of detected photons (can't hurt)
cfg.issaveref       = zerolayerflag;       % Also save diffuse reflectance
cfg.isnormalized    = 0;
cfg.issrcfrom0      = 0;

% Timing information
cfg.tstart          = 0;       % Time the simulation starts
cfg.tstep           = 5e-9;    % Step size
cfg.tend            = 5e-9;    % When to end.  The output will have [tstart:tstep:tend] slices at each of the different time points
%Scaling to Voxel
scale               = cfg.unitinmm;
xscaled             = x/scale;
yscaled             = y/scale;
zscaled             = z/scale;
edge                = cfg.margin/scale;
% simulation parameters & settings
cfg.srcpos = [1+edge, 1+edge, 0];           % source position
cfg.srcdir = [0,  0, 1];            % source direction
%cfg.detpos = [xscaled/2, yscaled/2, 0, yscaled/2];    % detector position and radius (x,y,z,r)*****************
cfg.bc='______001000';  % capture photons exiting from z=0 face
cfg.savedetflag='dpxvwm';
illumY           = yscaled - 2*edge;      % size of illumination y
illumX           = xscaled - 2*edge;
% trajectory check simulation calculations
if trajflag == 1
    bz=Dz/scale;
    bxc=Dpos(1)/scale;
    byc=Dpos(2)/scale;
    bx=Dx/scale;
    by=Dy/scale;
    bx0=bxc - bx/2;
    bx1=bxc + bx/2;
    by0=byc - by/2;
    by1=byc + by/2;
end
%==========================================================================  

%% Model Running
icx = 0;
idx = 0;
ifx = 0;

if zerolayerflag == 0
    expected_runs = length(musp_dynamic)*length(freq)*3;
    trix(expected_runs,length(mua_dynamic)) = struct();
    Tg                                      = zeros(expected_runs,length(mua_dynamic));
else
    expected_runs = length(musp_dynamic)*length(freq)*3*length(mua_dynamic);
    trix(expected_runs) = struct();
end

LUTtrack(length(freq))                  = struct();
%fluence(expected_runs)                  = struct('data',[],'stat',[], 'dref', []);
%detphoton(expected_runs)                = struct('detid', [], 'ppath', [], 'mom', [], 'p', [], 'v', [], 'w0', [], 'prop', [], 'unitinmm', [], 'data', []);
pix(length(freq))                       = struct();
T                                       = zeros(1,expected_runs);

% vol(expected_runs)                      = struct();
% seed(expected_runs)                     = struct();
% trajectory(expected_runs)               = struct();
% drefmc                                  = zeros(1,expected_runs/3);
% dreflog                                 = zeros(1,expected_runs/3);

%LUT allocation
MuaMatrix = zeros(length(musp_dynamic),length(mua_dynamic));
MuspMatrix = zeros(length(musp_dynamic),length(mua_dynamic));
M1 = zeros(length(musp_dynamic),length(mua_dynamic));
M2 = zeros(length(musp_dynamic),length(mua_dynamic));


gdev = gpuDevice(1);

tStart = tic;
for fy = freq % spatial frequencies , 0.0-0.1
    ifx = ifx + 1;
    
    for i = 1:length(musp_dynamic) % Running simulations for each scattering value
        air= [0,0,1,1]; % air (formatted as [mua, mus, g, n]) Layer #1
        zf = 1 + zerolayerflag;
        MuspMatrix(i,:) = musp_dynamic(i);
        mus_cur = musp_dynamic(i)/(1-g);
        % Using a layer of zeros to catch photons and determine dref
        if zerolayerflag == 1
            for j = 1:length(mua_dynamic)
                MuaMatrix(:,j) = mua_dynamic(j);
                idx = idx + 1;
                % Layer Properties [mua, mus, g, n]
                
                if multilayerflag == 1 && SACflag == 0
%                     mua     = 0.048; % Set to values used in Darren et al (0.096)
%                     musp    = 0.78; % set to values used in Darren et al
%                     skin_thick = 0.3125; % based on values used by Darren et al
%                     mua     = 8.68;
%                     musp    = 6.413 * (1-g);
                    mua     = 0.002793;
                    musp    = 1.168783; % DOSI measurements at 659nm for Karissa Phantom
                    skin_thick = 0.1;
                    mua     = mua * skin_thick;
                    musp    = musp * skin_thick;

                    mus     = musp/(1-g); % mu_s = (mu_s_prime) / (1-anisotropy)

                    layer2  = [mua, mus, g, n]; % outermost layer (e.g. epidermis OR epidermis+dermis)
                else
                    layer2  = [mua_dynamic(j), mus_cur, g, n]; % outermost layer (e.g. epidermis OR epidermis+dermis)
                end
                
                if multilayerflag == 1 && SACflag == 1
                    layer2  = epiderm_laylit;
                end
                
                if multilayerflag == 1 && SACflag == 2
                    layer2  = epiderm_laydrk;
                end

                layer3      = [mua_dynamic(j), mus_cur, g, n]; % middle layer (e.g. dermis)
                layer4      = [mua_dynamic(j), mus_cur, g, n]; % semi-infinite layer (e.g. adipose, tumor, etc.)
                
                layer2(1) 	= layer2(1)*scale;
                layer2(2) 	= layer2(2)*scale;
                layer3(1)	= layer3(1)*scale;
                layer3(2) 	= layer3(2)*scale;
                layer4(1)   = layer4(1)*scale;
                layer4(2)   = layer4(2)*scale;
                cfg.prop    = [     air;      
                                    layer2;      
                                    layer3;     
                                    layer4];     
                layerscaled(1)      = 1; %One voxel for air on top layer of volume
                layerscaled(2)      = layerthickness(2)/scale;
                layerscaled(3)      = layerthickness(3)/scale;
                cfg.vol = uint8(ones(xscaled,yscaled,zscaled))*numLayer; %Fills remaining space to requested box dimensions with layer 4
                cfg.vol(:,:,layerscaled(1)) = 0; %pad a layer of 0x to get diffuse reflectance
                cfg.vol(:,:,zf:zerolayerflag+(layerscaled(2))) = 1; %layer 2
                cfg.vol(:,:,(zf+layerscaled(2)):(zerolayerflag+layerscaled(2)+layerscaled(3))) = 2; % layer 3
                for phase = [0, 1/3, 2/3] % phase 
                    phs = phase*3 + 1;
                    % Approach to SFDI Pattern Projection
                    rng('shuffle');
                    kx = 0;
                    ky = fy * illumY * scale; % number of cycles
                     % Approach to SFDI pattern projection
                    cfg.srcparam1 = [illumX, 0, 0, kx+phase]; % [30, 0, 0] is offset of the 2nd corner to the 1st
                    cfg.srcparam2 = [0, illumY, 0, ky]; %  S=0.5*[1+M*cos(2*pi*(fx*x+fy*y)+phi0)], (0<=x,y,M<=1)
                    curdx = idx*3 - 2 + (phase*3);

                    % Run MCX
                    % Also track how long mcx takes to run
                    tRun = tic;
                    [fluence(phs),~]=mcxlab(cfg, 'cuda');
                    T(curdx) = toc(tRun);
                    Tmessage = ['Run ', num2str(curdx), ' Time elapsed: ', num2str(T(curdx)), ' seconds'];
                    disp(Tmessage);
                end
                [LUTtrack(ifx).data(i,j), ~, ~] = mcxdemodRd(fluence, cfg, 1);
            end
        end

        % Using detected photons instead
        if zerolayerflag == 0
            % Layer Properties
            layer2=[0, mus_cur, g, 1.4];
            layer3=[0, mus_cur, g, 1.4];
            layer4=[0, mus_cur, g, 1.4];
            layer2(1)=layer2(1)*scale;
            layer2(2)=layer2(2)*scale;
            layer3(1)=layer3(1)*scale;
            layer3(2)=layer3(2)*scale;
            layer4(1)=layer4(1)*scale;
            layer4(2)=layer4(2)*scale;
            cfg.prop = [    air;      
                            layer2;      
                            layer3;     
                            layer4];     
            layerscaled(2)=layerthickness(2)/scale;
            layerscaled(3)=layerthickness(3)/scale;
            cfg.vol = uint8(ones(xscaled,yscaled,zscaled))*numLayer; %Fills remaining space to requested box dimensions with layer 4
            cfg.vol(:,:,zf:zerolayerflag+(layerscaled(2))) = 1; %layer 2
            cfg.vol(:,:,(zf+layerscaled(2)):(zerolayerflag+layerscaled(2)+layerscaled(3))) = 2; % layer 3
            
            for phase = [0, 1/3, 2/3] % phase 
                % Approach to SFDI Pattern Projection
                phs = phase*3 + 1;
                rng('shuffle');
                kx = 0;
                ky = fy * illumY * scale; % number of cycles
                 % Approach to SFDI pattern projection
                cfg.srcparam1 = [illumX, 0, 0, kx+phase]; % [30, 0, 0] is offset of the 2nd corner to the 1st
                cfg.srcparam2 = [0, illumY, 0, ky]; %  S=0.5*[1+M*cos(2*pi*(fx*x+fy*y)+phi0)], (0<=x,y,M<=1)
                curdx = idx*3 - 2 + (phase*3);

                tRun = tic;
                [fluence(curdx),detphoton(phs)]=mcxlab(cfg, 'cuda');
                T(curdx) = toc(tRun);
                Tmessage = ['Run ', num2str(curdx), ' Time elapsed: ', num2str(T(curdx)), ' seconds'];
                disp(Tmessage);
                reset(gdev);

                % [fluence(idx),detphoton(idx),vol(idx),seed(idx),trajectory(idx)]=mcxlab(cfg, 'cuda');
                % detphoton(idx).detid(detphoton(idx).detid < 0) = 1;
                % for loop through different mua values for our LUT
                for j = 1:length(mua_dynamic)
                    % set up matrix to serve as our bins based on x & y exit
                    % positions - this recreates the diffuse reflectance setting,
                    % letting us continue to use demodulation
                    trix(curdx,j).dref = zeros(xscaled, yscaled);
                    trix(curdx,j).stat = fluence(curdx).stat;

                    % gpuArrays with detected photons to look at partial path 
                    % lengths, attenuated by mua of interest
                    toplayer_mua = 0.067;
                    midlayer_mua = mua_dynamic(j);
                    botlayer_mua = mua_dynamic(j);

                    tRun = tic;
                    gpDetP = gpuArray(detphoton(phs).w0);

                    gpTopPath = gpuArray(detphoton(phs).ppath(:,1));
                    gpMidPath = gpuArray(detphoton(phs).ppath(:,2));
                    gpBotPath = gpuArray(detphoton(phs).ppath(:,3));

                    % Running this code, but on a GPU
                    % detphoton(curdx.w0).*exp(-1*(                 ...
                    % detphoton(curdx).ppath(:,1).*toplayer_mua +   ...
                    % detphoton(curdx).ppath(:,2).*midlayer_mua +   ... 
                    % detphoton(curdx).ppath(:,3).*botlayer_mua));

                    gpTopPath = gpTopPath.*toplayer_mua;
                    gpMidPath = gpMidPath.*midlayer_mua;
                    gpBotPath = gpBotPath.*botlayer_mua;

                    gpPath = gpTopPath + gpMidPath + gpBotPath;
                    gpPath = gpPath.*-1;
                    gpPath = exp(gpPath);

                    gpWeights = gpPath.*gpDetP;

                    % sum weights based on exit position, add to matrix

                    gpXpos = gpuArray(detphoton(phs).p(:,1));
                    gpYpos = gpuArray(detphoton(phs).p(:,2));

                    gpXpos = ceil(gpXpos);
                    gpYpos = ceil(gpYpos);

                    gpWeightSums = zeros(140,100, "gpuArray");
                    for px = 1:xscaled
                        for py = 1:yscaled
                            gpWeightSums(px,py) = sum(gpWeights(gpXpos == px & gpYpos == py));
                        end
                    end
                    trix(curdx,j).dref = gather(gpWeightSums);

                    Tg(curdx,j) = toc(tRun);
                    Tmessage = ['Run ', num2str(curdx), ' GPU Cycle ', num2str(j), ' Time elapsed: ', num2str(Tg(curdx,j)), ' seconds'];
                    disp(Tmessage);
    %               mask = ceil(detphoton(curdx).p(:,1)) == px & ceil(detphoton(curdx).p(:,2)) == py;
    %               trix(curdx,j).dref(px,py) = sum(temp_weights(mask));
                end
            end
            %cwdref(idx)=sum(fluence(idx).dref,4);
            % once all photons have been processed accordingly, take that
            % matrix and demodulate it similarly to the boundary zero layer
            for j = 1:length(mua_dynamic)
                %pix(ifx).data=zeros(xscaled, yscaled);
                % store result in LUT
                [LUTtrack(ifx).data(i,j), ~ , ~] = mcxdemodRd(transpose(trix(:,j)), cfg, idx);
            end
        end
    end
end

tTotal = toc(tStart);
tAvg = sum(T)/expected_runs;
Tmessage = ['Average Run Time: ',num2str(tAvg),' seconds'];
disp(Tmessage);

%% Save Output
if multilayerflag == 1
    if SACflag == 1
        layer_info = "_SAC_multilayer_light";
    else
        if SACflag == 2
            layer_info = "_SAC_multilayer_tan";
        else
            layer_info = "_multilayer_toplayval_[" + num2str(mua) + ";" + num2str(musp) + "]";
        end
    end
else
    layer_info = "_homogenous";
end

n_sample = "_" + num2str(length(musp_dynamic)) + "x" + num2str(length(mua_dynamic));

LUT_name = 'my_LUT_Fourier' + layer_info + n_sample;

LUT.M1 = double(LUTtrack(1).data(:,:)); % Frequency of 0
LUT.M2 = double(LUTtrack(2).data(:,:)); % frequency of 0.1
LUT.Mua = MuaMatrix;
LUT.Musp = MuspMatrix;
save(('Generated_LUTs\' + LUT_name),'-struct','LUT');


%% Fluence+ Display
pf=1;
if rem(displayflag, 2) == 1
    if zerolayerflag == 0
        for fy = freq
            figure(pf)
            fcw=fluence((pf*3)-2).data*cfg.tstep;
            imagesc(log10(abs(squeeze(fcw(:,:,1)))))
            axis equal; colorbar
            %caxis([-5,-4])
            title(['flux ', num2str(freq(pf)), ' freq']);
            pf=pf+1;
        end
    else
        for fy = freq
            figure(pf)
            fcw=fluence((pf*3)-2).dref*cfg.tstep;
            imagesc(log10(abs(squeeze(fcw(:,:,1)))))
            axis equal; colorbar
            %caxis([-8,-5])
            title(['dref ', num2str(freq(pf)), ' freq']);
            pf=pf+1;
        end
    end
    
end
pft = pf;
    
    %% Secondary Displays
pf = pft;
if or(displayflag == 2, displayflag == 3)
    n=2; % which spatial frequency to display info about?
    if zerolayerflag == 0
        for i = 1:length(freq)
            figure(pf)
            heatmap(log10(LUT(i).data))
            
        end
    else
        figure(pf)
        plot(freq, drefmc)
        ylim([-2 0])
        xlim([0 0.5])
        ylabel("Log10 Rd")
        xlabel("Spatial freq. (mm^{-1})")
        yticks([-2 -1.5 -1 -0.5 0])
        yticklabels(["10^{-2}," "10^{-1.5}", "10^{-1}", "10^{-0.5}", "10^{0}"])
        grid on
        title(['L* = ', num2str(lstar), ' mm^{-1}'])
        pf=pf+1;
        
        figure(pf)
        fcw=pix(n).data*cfg.tstep;
        imagesc(log10(abs(squeeze(fcw(:,:,1)))))
        axis equal; colorbar
        %caxis([-12,-5.9])
        title(['demod ', num2str(freq(n)), ' freq']);
        pf=pf+1;
        
        illstart = 1+edge;
        illyend = yscaled - edge;
        illxend = xscaled - edge;
        
        figure(pf)
        mcv0=pix(n).data(xscaled/2,illstart:illyend);
        mcv1=mean(pix(n).data(illstart:illxend,illstart:illyend));
        mcv2=mean(fluence((n*3)-2).dref(illstart:illxend,illstart:illyend,1));
        mcv3=mean(fluence((n*3)-1).dref(illstart:illxend,illstart:illyend,1));
        mcv4=mean(fluence(n*3).dref(illstart:illxend,illstart:illyend,1));
        plot(1:illumY, log10(mcv0.*cfg.tstep),...
            1:illumY, log10(mcv1.*cfg.tstep),...
            1:illumY, log10(mcv2.*cfg.tstep),...
            1:illumY, log10(mcv3.*cfg.tstep),...
            1:illumY, log10(mcv4.*cfg.tstep))
        title([num2str(freq(n)), ' freq log10 cross-section comparisons'])
        legend("Demodulation", "Demod. Mean", "Phase Shift 0", "Phase Shift 1", "Phase Shift 2")
        pf=pf+1;
    end
    figure(pf)
    subplot(231);
    mcxpreview(cfg);title('domain preview');
    if zerolayerflag == 0
        subplot(232);
        histogram(detphoton((n*3)-2).ppath(:,1),20); title('partial path tissue#1');
    end
    subplot(233);
    plot(squeeze(fluence((n*3)-2).data(30/scale,30/scale,10/scale,:)),'-o');title('TPSF at [30,30,10]');
    subplot(234);
    imagesc(squeeze(log(fluence((n*3)-2).data(xscaled/2,:,:))));title('fluence at x=',x/2,'mm');
    subplot(235);
    imagesc(squeeze(log(fluence((n*3)-2).data(:,yscaled/2,:))));title('fluence at y=',y/2,'mm');
    subplot(236);
    imagesc(squeeze(log(fluence((n*3)-2).data(:,:,zscaled/2))));title('fluence at z=',z/2,'mm');
end