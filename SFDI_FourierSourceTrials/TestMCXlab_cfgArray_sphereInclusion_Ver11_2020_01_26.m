% fields with * are required; options in [] are the default values

set(groot,'defaultAxesLineStyleOrder','-|-.|:|--')

tic

clear;

fx = [0 0.3]; % mm^-1
fxstr = {'0','0.3'};
numfx = length(fx);

%Dimension of domain in voxels
xy = 500;
zz = 70;
%Generate voxel positions
[xi,yi,zi]=meshgrid(1:xy,1:xy,1:zz);
% sphere sit at bottom
dist=(xi-(xy/2)).^2+(yi-(xy/2)).^2+(zi-(21)).^2;

for ii = 1:numfx
    
    cfg(ii).nphoton=1e7;
    cfg(ii).maxdetphoton = 1e7; % save up to this many photons
    % %     cfg(ii).seed=randi(1e7);
    cfg(ii).seed=99999;
    
    cfg(ii).unitinmm = 0.1; %Each grid unit in mm
    cfg(ii).vol = ones(size(xi));
    cfg(ii).vol(dist<400) = 2; % 2 mm radius sphere
    %cfg.vol(:,:,60) = 0; % pad a layer of 0s to get diffuse reflectance
    cfg(ii).vol = uint8(cfg(ii).vol);
    
    cfg(ii).gpuid=1;
    cfg(ii).autopilot=1;
    
    cfg(ii).prop=[0 0 1 1; 0.005 8 0.9 1.4; 0.1 8 0.9 1.4];
    
    % % cfg.issaveexit=1; %Save the position and direction of detected photons (can't hurt)'
    % % cfg.issaveref=1; %Also save diffuse reflectance
        
    cfg(ii).srctype='fourier';
    cfg(ii).srcpos=[0 0 zz+1];
    cfg(ii).srcdir=[0 0 -1];
    cfg(ii).srcparam1=[xy 0 0 fx(ii)*cfg(ii).unitinmm*xy];
    cfg(ii).srcparam2=[0 xy 0 0];
    
    cfg(ii).tstart=0;
    cfg(ii).tend=1e-9;
    cfg(ii).tstep=1e-9;
    
    % % cfg.isreflect=0; %Don''t reflect from outside the boundary (assume photons can escape)
    % % cfg.isrefint=1; %Do consider internal reflections
    cfg(ii).bc = 'aaaaar';
    cfg(ii).outputtype = 'flux';
    cfg(ii).isspecular = 1; % 1-calculate specular reflection if source is outside, [0] no specular reflection
    % % cfg.isnormalized = 0; % [1]-normalize the output fluence to unitary source, 0-no reflection
    % % cfg.issrcrom0 = 1; % 1-first voxel is [0 0 0], [0]- first voxel is [1 1 1]
    
end


%% run mcx

[flux,detpt,vol,seeds,traj]=mcxlab(cfg);


%% plots

for ii = 1:numfx
    
% %     fluence(:,:,:,ii) = flux(ii).data * cfg(ii).tstep / ( flux(ii).stat.normalizer * (flux(ii).stat.energyabs/flux(ii).stat.energytot) );
    fluence(:,:,:,ii) = flux(ii).data * cfg(ii).tstep / flux(ii).stat.normalizer;
    ff = fluence(:,:,:,ii);
    
    for jj = 1:zz
        
        meanFluence(jj,ii) = mean2( ff( 0.1*xy:0.9*xy , 0.1*xy:0.9*xy , zz - (jj-1) ) );
        
    end
    
    % Plot entire volumes
    figure;
    hs=slice(log10(abs(double(fluence(:,:,:,ii)))),1,1,zz);
    caxis([-9,-7])
    set(hs,'linestyle','none');
    axis equal; colorbar;box on;
    title(['fluence ' num2str(fx(ii)) ' mm^{-1}']);
    
    % Plot xz slice
    figure;
    imagesc(log10(abs(double(squeeze(fluence(:,xy/2,:,ii))))));
    caxis([-8,-6.8])
    colorbar;box on;
    title(['fluence ' num2str(fx(ii)) ' mm^{-1}']);
    xlabel('depth (10 = 1 mm)')
    
end

depth = cfg(1).unitinmm:cfg(1).unitinmm:cfg(1).unitinmm*zz;
depth = depth';
depthMat = repmat(depth,1,numfx);

figure
plot(depthMat,meanFluence)
title('Fluence')
xlabel('Depth (mm)')
legend(fxstr)

% % norm1 = flux1.stat.normalizer;
% % norm2 = flux2.stat.normalizer;
% % 
% % fcw1 = flux1.data*cfg.tstep/(norm1* (flux1.stat.energyabs/flux1.stat.energytot) );
% % fcw2 = flux2.data*cfg.tstep/(norm2* (flux2.stat.energyabs/flux2.stat.energytot) );
% % 

% % %%Plot entire volumes
% % figure;
% % hs=slice(log10(abs(double(fluence(:,:,:,1)))),1,1,zz);
% % caxis([-9,-7])
% % set(hs,'linestyle','none');
% % axis equal; colorbar;box on;
% % title(['fluence ' num2str(fx(1)) ' mm^{-1}']);
% % 
% % 
% % %%Plot xz slice
% % figure;
% % % % hs=slice(log10(abs(double(fluence(:,:,:,1)))),250,[],[]);
% % imagesc(log10(abs(double(squeeze(fluence(:,250,:,1))))));
% % caxis([-8,-6.8])
% % colorbar;box on;
% % title(['fluence ' num2str(fx(1)) ' mm^{-1}']);
% % xlabel('depth (10 = 1 mm)')


% % % % 
% % % % figure;
% % % % hs=slice(abs(double(fcw2)),1,1,zz);
% % % % % % caxis([0,1e-5])
% % % % set(hs,'linestyle','none');
% % % % axis equal; colorbar;box on;
% % % % title(['a spatial frequency domain source ' num2str(fx2) ' mm^{-1}']);
% % 
% % %%Plot mean fluence(depth)
% % f1 = zeros(1,zz);
% % f2 = f1;
% % for ii = 1:zz
% %     f1(ii) = mean2(fcw1(:,:,zz - (ii-1)));
% %     f2(ii) = mean2(fcw2(:,:,zz - (ii-1)));
% % end
% % % % f1 = f1/f1(1);
% % % % f2 = f2/f2(2);
% % figure
% % plot(f1);
% % hold on
% % plot(f2);
% % title('normalized fluence')
% % xlabel('depth (mm)')
% % legend(num2str(fx1),num2str(fx2))

toc