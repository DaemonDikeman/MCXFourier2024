%% SFDI 2-layer white MLMC 
% Karissa Tilbury & Christian Crane 07_27_2021
close all
clear
clear cfg


%% CONTROL PANEL 
%info=mcxlab('gpuinfo');
%Total Volume Creation (mm)
x=100; %Total x dimension of box
y=60; %Total y dimension of box
z=60; %Total z dimension of box
cfg.unitinmm     = 0.2; % length of a voxel in mm can't be larger than layer thickness estimation (ex. 5.1mm thickness would require voxel size 0.1 or smaller)
cfg.nphoton      = 1e6;     % number of photons
cfg.srctype = 'fourier'; 

%Layer Creation
numLayer = 3; % number of tissue layers (not including top air layer (layer 1))
layerthickness(2)=1; %thickness of epidermis in mm (layer 1 is air)
layerthickness(3)=2; %thickness of dermis in mm
layer2=[0.4, 64.13, 0.82, 1.4];% Layer #2 Epidermis optical properties [mua,mus,g,n]
layer3=[1.2, 26.27, 0.82, 1.4]; % Layer #3 Dermis optical properties [mua, mus,g,n]
layer4=[.6,  5,  0.75,  1.44];% Layer #4 Adipose optical properties (remaining volume)

%Trajectory Detector Creation
Dx=20;  %Detector size in mm along x axis
Dy=20; %Detector size in mm along y axis
Dz=600; %Detector distance from box in mm
Dpos=[x/2 y/2];

%% CONSTANTS DO NOT CHANGE=================================================
cfg.gpuid        = 1;
cfg.autopilot    = 1;
cfg.isreflect    = 0;       % Don't reflect from outside the boundary
cfg.isrefint     = 1;       % Do consider internal reflections
cfg.maxdetphoton = 1e7;     % save up to this many photons (effects memory significantly)
cfg.issaveexit   = 0;       % Save the position and direction of detected photons (can't hurt)
cfg.issaveref    = 1;       % Also save diffuse reflectance
layerscaled(1)=1; %One voxel for air on top layer of volume
% Timing information
cfg.tstart = 0;       % Time the simulation starts
cfg.tstep  = 5e-9;   % Step size
cfg.tend   = 5e-9;    % When to end.  The output will have [tstart:tstep:tend] slices at each of the different time points
%Scaling to Voxel
scale  = cfg.unitinmm;
xscaled=x/scale;
yscaled=y/scale;
zscaled=z/scale;
layer2(1)=layer2(1)*scale;
layer2(2)=layer2(2)*scale;
layer3(1)=layer3(1)*scale;
layer3(2)=layer3(2)*scale;
layer4(1)=layer4(1)*scale;
layer4(2)=layer4(2)*scale;
% Layer Properties
air= [0,0,1,1];% air (formatted as [mua, mus, g, n]) Layer #1
cfg.prop = [    air;      
                layer2;      
                layer3;     
                layer4];     
layerscaled(2)=layerthickness(2)/scale;
layerscaled(3)=layerthickness(3)/scale;
cfg.vol = uint8(ones(xscaled,yscaled,zscaled))*numLayer; %Fills remaining space to requested box dimensions with layer 4
%cfg.vol(:,:,layerscaled(1)) = 0; %pad a layer of 0x to get diffuse reflectance
cfg.vol(:,:,1:(layerscaled(2))) = 1; %layer 2
cfg.vol(:,:,(1+layerscaled(2)):(1+layerscaled(2)+layerscaled(3))) = 2; % layer 3   
% simulation parameters & settings
cfg.srcpos = [0, 0, 0];           % source position
cfg.srcdir = [0,  0, 1];            % source direction
%cfg.detpos = [xscaled/2, yscaled/2, 0, yscaled/2];    % detector position and radius (x,y,z,r)*****************
cfg.bc='______001000';  % capture photons exiting from z=0 face
cfg.savedetflag='dpxv';
illumY           = yscaled;      % size of illumination y
% trajectory check simulation calculations
bz=Dz/scale;
bxc=Dpos(1)/scale;
byc=Dpos(2)/scale;
bx=Dx/scale;
by=Dy/scale;
bx0=bxc - bx/2;
bx1=bxc + bx/2;
by0=byc - by/2;
by1=byc + by/2;
%==========================================================================  

%% Model Running
idx = 0;
for fy = [0.05] % spatial frequencies ,0.1, 0.15, 0.2
    for phase = [0, 1/3, 2/3] % phase 
        
        kx = 0;
        ky = fy * illumY; % number of cycles
         % Approach to SFDI pattern projection
        cfg.srcparam1 = [xscaled, 0, 0, kx+phase]; % [30, 0, 0] is offset of the 2nd corner to the 1st
        cfg.srcparam2 = [0, illumY, 0, ky]; %  S=0.5*[1+M*cos(2*pi*(fx*x+fy*y)+phi0)], (0<=x,y,M<=1)                   

        
        idx = idx + 1;
        %[fluence(idx),detphoton(idx),vol(idx),seed(idx),trajectory(idx)]=mcxlab(cfg);
        [fluence(idx),detphoton(idx)]=mcxlab(cfg, 'cuda');
    end
end

%% Trajectory Prep
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

%% Trajectory calculator
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
%% Trajectory Display
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

%% Fluence+ Display
figure(1)
plot3(detphoton(1).p(:,1),detphoton(1).p(:,2),detphoton(1).p(:,3),'r.');
view([0 0 1])
title('photon detections, top down')
figure(2)
fcw=fluence(1).data*cfg.tstep;
imagesc(log10(abs(squeeze(fcw(:,:,1)))))
axis equal; colorbar
title('0.05 freq');
%figure(3)
%plot3(detphoton(4).p(:,1),detphoton(4).p(:,2),detphoton(4).p(:,3),'r.');
%view([0 0 1])
% figure(4)
% fcw=fluence(4).data*cfg.tstep;
% imagesc(log10(abs(squeeze(fcw(:,:,1)))))
% axis equal; colorbar
% title('0.05 freq');
%figure(5)
%plot3(detphoton(7).p(:,1),detphoton(7).p(:,2),detphoton(7).p(:,3),'r.');
%view([0 0 1])
% figure(6)
% fcw=fluence(7).data*cfg.tstep;
% imagesc(log10(abs(squeeze(fcw(:,:,1)))))
% axis equal; colorbar
% title('0.10 freq');
%figure(7)
%mcxpreview(cfg); title('domain preview')
