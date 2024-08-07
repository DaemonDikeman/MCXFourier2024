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
x                = 140;     % Total x dimension of box in mm
y                = 100;     % Total y dimension of box in mm
z                = 20;      % Total z dimension of box in mm
lstar            = 2;       % 1/(mua+musp) ; constant ratio of musp/mua = 100 used
g                = 0.71;    % anisotropy; base alt: 0.71/0.9
n                = 1.4;    % base alt: 1.330/1.4
cfg.margin       = 20;      % Margin between edge of box and edge of pattern source in mm
cfg.unitinmm     = 1;       % length of a voxel in mm can't be larger than layer thickness estimation (ex. 5.1mm thickness would require voxel size 0.1 or smaller)
cfg.nphoton      = 1e7;     % number of photons
cfg.srctype      = 'fourier';    % type of pattern/source to project; fourier pattern is our interest
freq             =  0:0.01:0.5; % set of spatial frequencies to run in mm^-1
%       [0, 0.05/3, 0.10/3, 0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5]
%       0:0.01:0.5
layer2mult = 1;
layer3mult = 1;
%cfg.minenergy = 0.0001;
zerolayerflag = 1;          % determines whether we use a layer of zero-value voxels to save diffuse reflectance; 1 is yes
trajflag = 0;               % determines whether to calculate and display trajectories or not.
displayflag = 0;            % determines whether or not to display results appropriate to the run type
multilayerflag = 0;

% foldername = "SFDI_dref_nphoton_e7_test\\";
foldername = "SFDI_dref_nphoton_e7_test\\";

if multilayerflag == 1
    filename = "SFDI_multilayer_test_drefval_2layermodelBCthin_Lstar_" + num2str(lstar) + "_nphoton_" + num2str(cfg.nphoton, '%.1e') + "_voxel_" + num2str(cfg.unitinmm) + "mm.csv";
end 
if multilayerflag == 0
    filename = "SFDI_drefval_test_g_" + num2str(g,'%.2f') + "_Lstar_" + num2str(lstar) + "_nphoton_" + num2str(cfg.nphoton, '%.1e') + "_vx_" + num2str(cfg.unitinmm) + "mm.csv";
    %filename = "testing.csv";
end

%% Layer Creation
numLayer = 3; % number of tissue layers (not including top air layer (layer 1))
layerthickness(2)=1; %thickness of epidermis in mm (layer 1 is air)
layerthickness(3)=2; %thickness of dermis in mm

% constant ratio of musp/mua = 100
% Lstar = 1/(musp+mua)
mua=1/(lstar*101); % from above eq.: Lstar = 1/(101*mua)
musp=100*mua; % from above eq.
mus=musp/(1-g); % mu_s = (mu_s_prime) / (1-anisotropy)

% [0.067, 0.64, 0.82, 1.4] A - Dermis
% [0.043, 2.63, 0.82, 1.4] B - Epidermis
% [0.033, 0.5, 0.75, 1.44] C - Adipose

layer2=[mua, mus, g, n];% Layer #2 Epidermis optical properties [mua,mus,g,n]
if multilayerflag == 1
    layer2=[0.067, 0.64, 0.82, 1.4];
    mua2=1/(lstar*layer2mult*101);
    musp2=100*mua2;
    mus2=musp2/(1-g); %musp = mus*(1-g)
    mua3=1/(lstar*layer3mult*101);
    musp3=100*mua3;
    mus3=musp3/(1-g); %musp = mus*(1-g)
    layer3=[mua2, mus2, g, 1.33]; % Layer #3 Dermis optical properties [mua, mus,g,n]
    layer4=[mua3, mus3, g, 1.33];% Layer #4 Adipose optical properties (remaining volume)
else
    layer3=[mua, mus, g, n]; % Layer #3 Dermis optical properties [mua, mus,g,n]
    layer4=[mua, mus, g, n];% Layer #4 Adipose optical properties (remaining volume)
end

%% Photon diffusion approximation
% if multilayerflag == 0
%     t = 0;
%     for f = freq
%         t = t + 1;
%         approxRd(t) = diffApproxSFD(mua, musp, n, f);
%     end
% end
% 
% appfilename = "SFDI_approximation_Lstar_" + num2str(lstar) + ".csv";
% 
% writematrix(approxRd,foldername+appfilename);

%% Trajectory Detector Creation
if trajflag == 1
    Dx=20;  %Detector size in mm along x axis
    Dy=20; %Detector size in mm along y axis
    Dz=600; %Detector distance from box in mm
    Dpos=[x/2 y/2];
end

%% CONSTANTS DO NOT CHANGE=================================================
cfg.gpuid        = 1;
cfg.autopilot    = 1;
cfg.isreflect    = 0;       % Don't reflect from outside the boundary
cfg.isrefint     = 1;       % Do consider internal reflections
cfg.maxdetphoton = 1e7;     % save up to this many photons (effects memory significantly) !DO NOT INCREASE ABOVE 1E^8!
cfg.issaveexit   = 0;       % Save the position and direction of detected photons (can't hurt)
cfg.issaveref    = zerolayerflag;       % Also save diffuse reflectance
cfg.isnormalized = 0;
layerscaled(1)=1; %One voxel for air on top layer of volume
% Timing information
cfg.tstart = 0;       % Time the simulation starts
cfg.tstep  = 5e-9;    % Step size
cfg.tend   = 5e-9;    % When to end.  The output will have [tstart:tstep:tend] slices at each of the different time points
%Scaling to Voxel
scale  = cfg.unitinmm;
xscaled=x/scale;
yscaled=y/scale;
zscaled=z/scale;
edge=cfg.margin/scale;
% layer2(1)=layer2(1)*scale;
% layer2(2)=layer2(2)*scale;
% layer3(1)=layer3(1)*scale;
% layer3(2)=layer3(2)*scale;
% layer4(1)=layer4(1)*scale;
% layer4(2)=layer4(2)*scale;

% Layer Properties
air= [0,0,1,1];% air (formatted as [mua, mus, g, n]) Layer #1
cfg.prop = [    air;      
                layer2;      
                layer3;     
                layer4];     
layerscaled(2)=layerthickness(2)/scale;
layerscaled(3)=layerthickness(3)/scale;
cfg.vol = uint8(ones(xscaled,yscaled,zscaled))*numLayer; %Fills remaining space to requested box dimensions with layer 4
if zerolayerflag == 1
    cfg.vol(:,:,layerscaled(1)) = 0; %pad a layer of 0x to get diffuse reflectance
end
zf = 1 + zerolayerflag;
cfg.vol(:,:,zf:zerolayerflag+(layerscaled(2))) = 1; %layer 2
cfg.vol(:,:,(zf+layerscaled(2)):(zerolayerflag+layerscaled(2)+layerscaled(3))) = 2; % layer 3   
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
idx = 0;
icx = 0;

% fluence = struct(length(freq)*3);
% detphoton = zeros(length(freq)*3);

% if zerolayerflag == 0
%     drefmc = zeros(length(freq)*3);
% end
% if zerolayerflag == 1
%     drefmc = zeros(length(freq));
%     pix = zeros(length(freq));
%     dreflog = zeros (length(freq));
% end
tStart = tic;
for fy = freq % spatial frequencies , 0.1-0.3
    for phase = [0, 1/3, 2/3] % phase 
        rng('shuffle');
        kx = 0;
        ky = fy * illumY * scale; % number of cycles
         % Approach to SFDI pattern projection
        cfg.srcparam1 = [illumX, 0, 0, kx+phase]; % [30, 0, 0] is offset of the 2nd corner to the 1st
        cfg.srcparam2 = [0, illumY, 0, ky]; %  S=0.5*[1+M*cos(2*pi*(fx*x+fy*y)+phi0)], (0<=x,y,M<=1)                   

        idx = idx + 1;
        
        tRun = tic;
        [fluence(idx),detphoton(idx)]=mcxlab(cfg, 'cuda');
        T(idx) = toc(tRun);
        Tmessage = ['Run ', num2str(idx), ' Time elapsed: ', num2str(T(idx)), ' seconds'];
        disp(Tmessage);
        
        %[fluence(idx),detphoton(idx),vol(idx),seed(idx),trajectory(idx)]=mcxlab(cfg, 'cuda');
        
        %detphoton(idx).detid(detphoton(idx).detid < 0) = 1;
        if zerolayerflag == 0
            drefmc(idx) = mcxcwbbdref(detphoton(idx), cfg, fluence(idx));
        end
        %cwdref(idx)=sum(fluence(idx).dref,4);
    end
    if zerolayerflag == 1
        icx = icx + 1;
        pix(icx).data=zeros(xscaled, yscaled);
        [drefmc(icx), dreflog(icx), pix(icx).data] = mcxdemodRd(fluence, cfg, icx);
    end
end
tTotal = toc(tStart);
tAvg = sum(T)/idx;
Tmessage = ['Average Run Time: ',num2str(tAvg),' seconds'];
disp(Tmessage);

%% Trajectory Prep
if trajflag == 1
    fn = fieldnames(detphoton);
    propf=find(strcmp('prop',fn));
    if(~isempty(propf))
        fn(propf,:)=[];
    end
    unitf=find(strcmp('unitinmm',fn));
    if(~isempty(unitf))
        fn(unitf,:)=[];
    end
    dataf=find(strcmp('data',fn));
    if(~isempty(dataf))
        fn(dataf,:)=[];
    end
end

%% Trajectory calculator
% walon = zeros(length(detphoton));
% detset = zeros(length(detphoton));
if trajflag == 1
    for i = 1:length(detphoton)
        walon(i).t=bz./detphoton(i).v(:,3);
        walon(i).tx=detphoton(i).p(:,1)+detphoton(i).v(:,1).*walon(i).t;
        walon(i).ty=detphoton(i).p(:,2)+detphoton(i).v(:,2).*walon(i).t;
        walon(i).b = (walon(i).tx>bx0) & (walon(i).tx<bx1) & (walon(i).ty>by0) & (walon(i).ty<by1);
        for k=1:numel(fn)
            detset(i).(fn{k})=detphoton(i).(fn{k})(walon(i).b,:);
        end
        detset(i).fp=[walon(i).tx(walon(i).b) walon(i).ty(walon(i).b) repmat(bz,sum(walon(i).b),1)];

    end
end
%% Trajectory Display
if trajflag == 1 && displayflag == 1
    detptraj=length(detset(1).p);
    %displaying trajectories
    %Xf = detset(1).p(1:detptraj,1)+(10*detset(1).v(1:detptraj,1)); % Finding end points after 10 additional time steps for X
    detX = [detset(1).p(1:detptraj,1) detset(1).fp(1:detptraj,1)]; %
    %Yf = detset(1).p(1:detptraj,2)+(10*detset(1).v(1:detptraj,2));
    detY = [detset(1).p(1:detptraj,2) detset(1).fp(1:detptraj,2)];
    %Zf = detset(1).p(1:detptraj,3)+(10*detset(1).v(1:detptraj,3));
    detZ = [detset(1).p(1:detptraj,3) detset(1).fp(1:detptraj,3)];
    figure(3)
    plot3(detX',detY',detZ')
    hold on
    plot3(detX',detY',detZ','.')
    hold off
end

%% Save Output
%writematrix(drefmc,foldername+filename);


%% Fluence+ Display
 pf=1;
if or(displayflag == 1, displayflag == 3)
    % figure(1)
    % plot3(detphoton(1).p(:,1),detphoton(1).p(:,2),detphoton(1).p(:,3),'r.');
    % view([0 0 1])
    % title('photon detections, top down')
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
        drefmean = zeros(length(freq),1);
        for i = 1:length(freq)
            avg = 0;
            for j = 0:2
                avg = avg + drefmc((i*3)-j);
            end
            drefmean(i) = avg;
        end
        
        figure(pf)
        plot(freq, log10(drefmean))
        %ylim([-2 0])
        %xlim([0 0.5])
        ylabel("Log10 Rd")
        xlabel("Spatial freq. (mm^{-1})")
        %yticks([-2 -1.5 -1 -0.5 0])
        %yticklabels(["10^{-2}," "10^{-1.5}", "10^{-1}", "10^{-0.5}", "10^{0}"])
        grid on
        title(['L* = ', num2str(lstar), ' mm^{-1}'])
        pf=pf+1;
        
        figure(pf)
        fcw=fluence((n*3)-1).data*cfg.tstep;
        imagesc(log10(abs(squeeze(fcw(:,:,1)))))
        axis equal; colorbar
        %caxis([-5,-4])
        title(['flux ', num2str(freq(4)), ' freq phase 2']);
        pf=pf+1;

        figure(pf)
        fcw=fluence(n*3).data*cfg.tstep;
        imagesc(log10(abs(squeeze(fcw(:,:,1)))))
        axis equal; colorbar
        %caxis([-5,-4])
        title(['flux ', num2str(freq(n)), ' freq phase 3']);
        pf=pf+1;
    else
        figure(pf)
        plot(freq, dreflog)
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
        fcw=fluence((n*3)-1).dref*cfg.tstep;
        imagesc(log10(abs(squeeze(fcw(:,:,1)))))
        axis equal; colorbar
        %caxis([-12,-5.9])
        title(['flux ', num2str(freq(n)), ' freq phase 2']);
        pf=pf+1;

        figure(pf)
        fcw=fluence(n*3).dref*cfg.tstep;
        imagesc(log10(abs(squeeze(fcw(:,:,1)))))
        axis equal; colorbar
        %caxis([-12,-5.9])
        title(['flux ', num2str(freq(n)), ' freq phase 3']);
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
        hist(detphoton((n*3)-2).ppath(:,1),20); title('partial path tissue#1');
    end
    subplot(233);
    plot(squeeze(fluence((n*3)-2).data(30/scale,30/scale,10/scale,:)),'-o');title('TPSF at [30,30,10]');
    subplot(234);
    imagesc(squeeze(log(fluence((n*3)-2).data(xscaled/2,:,:))));title('fluence at x=70mm');
    subplot(235);
    imagesc(squeeze(log(fluence((n*3)-2).data(:,yscaled/2,:))));title('fluence at y=50mm');
    subplot(236);
    imagesc(squeeze(log(fluence((n*3)-2).data(:,:,zscaled/2))));title('fluence at z=10mm');
end

% figure(5)
% fcw=fluence(7).dref*cfg.tstep;
% imagesc(log10(abs(squeeze(fcw(:,:,1)))))
% axis equal; colorbar
% caxis([-5,-2])
% title('dref 0.10 freq');
% 
% figure(6)
% fcw=fluence(7).data*cfg.tstep;
% imagesc(log10(abs(squeeze(fcw(:,:,1)))))
% axis equal; colorbar
% caxis([-5,-2])
% title('flux 0.10 freq');
% figure(4)
% fcw=fluence(10).data*cfg.tstep;
% imagesc(log10(abs(squeeze(fcw(:,:,1)))))
% axis equal; colorbar
% title('0.15 freq');
% figure(5)
% fcw=fluence(13).data*cfg.tstep;
% imagesc(log10(abs(squeeze(fcw(:,:,1)))))
% axis equal; colorbar
% title('0.20 freq');
% figure(3)
% plot3(detphoton(4).p(:,1),detphoton(4).p(:,2),detphoton(4).p(:,3),'r.');
% view([0 0 1])
% figure(4)
% fcw=fluence(4).data*cfg.tstep;
% imagesc(log10(abs(squeeze(fcw(:,:,1)))))
% axis equal; colorbar
% title('0.05 freq');
% figure(5)
% plot3(detphoton(7).p(:,1),detphoton(7).p(:,2),detphoton(7).p(:,3),'r.');
% view([0 0 1])
% figure(6)
% fcw=fluence(7).data*cfg.tstep;
% imagesc(log10(abs(squeeze(fcw(:,:,1)))))
% axis equal; colorbar
% title('0.10 freq');
%figure(7)
mcxpreview(cfg); title('domain preview')
