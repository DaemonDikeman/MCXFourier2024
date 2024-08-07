close all
clear all

import numpy.*
import matplotlib.pyplot.*
import matplotlib.cm.ScalarMappable.*

% LUT1 = load('C:\Users\daemo\Documents\MATLAB\Generated_LUTs\my_LUT_Fourier_flipped_multilayer_botlayval_[0.001018;0.11308]_100x100.mat'); %variable top layer, static bottom layer
% LUT2 = load('C:\Users\daemo\Documents\MATLAB\Generated_LUTs\my_LUT_Fourier_SAC_multilayer_light_100x100.mat'); %normal multilayer LUT: static top layer, varying bottom layer
LUT1 = load('C:\Users\daemo\Documents\MATLAB\Generated_LUTs\my_LUT_Fourier_flipped_multilayer_botlayval_[0.002793;1.1688]_100x100.mat'); %variable top layer, static bottom layer
LUT2 = load('C:\Users\daemo\Documents\MATLAB\Generated_LUTs\my_LUT_Fourier_multilayer_toplayval_[0.0002793;0.11688]_100x100.mat'); %normal multilayer LUT: static top layer, varying bottom layer
% LUT1 = load('C:\Users\daemo\Documents\MATLAB\Generated_LUTs\my_LUT_Fourier_flipped_SAC_multilayer_tan_100x100.mat'); %variable top layer, static bottom layer
% LUT2 = load('C:\Users\daemo\Documents\MATLAB\Generated_LUTs\my_LUT_Fourier_SAC_multilayer_tan_100x100.mat'); %normal multilayer LUT: static top layer, varying bottom layer

contset = linspace(0,0.9,10);
contset2 = linspace(0, 0.9, 46);
wv = 659;
[xq, yq] = meshgrid(0.01:0.01:1, 0.01:0.01:1);
[xq2, yq2] = meshgrid(0.0095:0.0095:0.95, 0.0064:0.0064:0.64);


%% Comparing LUTs

LUT = LUT2;
f = figure(1);
f.Position = [10,10,600,400];
t = tiledlayout(2,2);
ax1 = nexttile;
hold on
contourf(LUT.Mua, LUT.Musp, LUT.M2(:,:),contset2, 'LineStyle','none')
contour(LUT.Mua, LUT.Musp, LUT.M2(:,:),contset, 'black','ShowText','on')
hold off
title('Static Top Layer AC')
caxis([0,0.9]);

ax2 = nexttile;
hold on
contourf(LUT.Mua, LUT.Musp, LUT.M1(:,:),contset2, 'LineStyle','none')
contour(LUT.Mua, LUT.Musp, LUT.M1(:,:),contset, 'black','ShowText','on')
hold off
title('Static Top Layer DC')
caxis([0,0.9]);

ax3 = nexttile;
hold on
contourf(LUT1.Mua, LUT1.Musp, LUT1.M2(:,:),contset2, 'LineStyle','none')
contour(LUT1.Mua, LUT1.Musp, LUT1.M2(:,:),contset, 'black','ShowText','on')
hold off
title('Dynamic Top Layer AC')
caxis([0,0.9]);

ax4 = nexttile;
hold on
contourf(LUT1.Mua, LUT1.Musp, LUT1.M1(:,:),contset2, 'LineStyle','none')
contour(LUT1.Mua, LUT1.Musp, LUT1.M1(:,:),contset, 'black','ShowText','on')
hold off
title('Dynamic Top Layer DC')
caxis([0,0.9]);

cb = colorbar;
cb.Layout.Tile = 'east';
cb.Label.String = 'Diffuse Reflectance';
cb.Label.FontSize = 12;

t.TileSpacing = 'compact';
t.Padding = 'compact';
title(t,'Diffuse Reflectance vs. Absorbance and Scattering')
xlabel(t,'\mu_a (mm^{-1})')
ylabel(t,'\mu_s'' (mm^{-1})')

%% Basic LUT inversions

absset1 = linspace(0,0.2,11);
absset2 = linspace(0,0.2,51);
absset3 = linspace(0,10,11);
absset4 = linspace(0,10,51);
scaset1 = linspace(0,5,11);
scaset2 = linspace(0,5,51);
scaset3 = linspace(0,30,11);
scaset4 = linspace(0,30,51);

clear op_fit_maps

g = 2;
for i=1:length(wv)

    disp(['Processing ' num2str(wv) ' nm...'])
    
    LUT = LUT1; % Dynamic Top Layer
    op_fit_maps1(:,:,i,1)=griddata(LUT.M1(:), LUT.M2(:), LUT.Mua(:), xq2, yq2, 'linear');
    op_fit_maps1(:,:,i,2)=griddata(LUT.M1(:), LUT.M2(:), LUT.Musp(:), xq2, yq2,'linear');

    LUT = LUT2; % Statis Top Layer
    op_fit_maps2(:,:,i,1)=griddata(LUT.M1(:), LUT.M2(:), LUT.Mua(:), xq2, yq2, 'linear');
    op_fit_maps2(:,:,i,2)=griddata(LUT.M1(:), LUT.M2(:), LUT.Musp(:), xq2, yq2,'linear');

    f = figure(g);
    f.Position = [10,10,900,600];
    t = tiledlayout(2,2);
    
    nexttile
    hold on
    contourf(xq2(1,:),yq2(:,1),op_fit_maps2(:,:,i,1),absset2,'LineStyle','none');
    contour(xq2(1,:),yq2(:,1),op_fit_maps2(:,:,i,1),absset1,'black','ShowText','on');
    caxis([0,0.2]);
    colorbar
    hold off
    title('Static Top Layer \mu_a (mm^{-1})')

    nexttile
    hold on
    [m,n] = contourf(xq2(1,:),yq2(:,1),op_fit_maps2(:,:,i,2),scaset2,'LineStyle','none');
    [a,b] = contour(xq2(1,:),yq2(:,1),op_fit_maps2(:,:,i,2),scaset1,'black','ShowText','on');
    colorbar
    caxis([0,5]);
    hold off
    title('Static Top Layer \mu_s'' (mm^{-1})')

    nexttile
    hold on
    contourf(xq2(1,:),yq2(:,1),op_fit_maps1(:,:,i,1),absset4,'LineStyle','none');
    contour(xq2(1,:),yq2(:,1),op_fit_maps1(:,:,i,1),absset3,'black','ShowText','on');
    caxis([0,10]);
    colorbar
    hold off
    title('Dynamic Top Layer \mu_a (mm^{-1})')
    
    nexttile
    hold on
    contourf(xq2(1,:),yq2(:,1),op_fit_maps1(:,:,i,2),scaset4,'LineStyle','none');
    contour(xq2(1,:),yq2(:,1),op_fit_maps1(:,:,i,2),scaset3,'black','ShowText','on');
    caxis([0,30]);
    colorbar
    hold off
    title('Dynamic Top Layer \mu_s'' (mm^{-1})')
    
    t.TileSpacing = 'compact';
    t.Padding = 'compact';
    title(t,'Diffuse Reflectance vs. Absorbance and Scattering')
    xlabel(t,'R_d (f_x = 0 mm^{-1})')
    ylabel(t,'R_d (f_x = 0.1 mm^{-1})')
    
    g = g + 1;
end


%% Compared LUT Inversions
% op_fit_maps3 = op_fit_maps2(:,:,:,:) - op_fit_maps1(:,:,:,:);
% 
% op_fit_maps4 = op_fit_maps2(:,:,:,:)./op_fit_maps1(:,:,:,:);

% bpavMua     = 0.001018;
% bpavMusp    = 0.628212;

for i=1:length(wv)
    op_fit_maps3(:,:,i,1) = griddata(LUT2.M1(:), LUT2.M2(:), LUT2.Mua(:), LUT1.M1, LUT1.M2, 'linear');
    op_fit_maps3(:,:,i,2) = griddata(LUT2.M1(:), LUT2.M2(:), LUT2.Musp(:), LUT1.M1, LUT1.M2, 'linear');
end

sz = size(op_fit_maps3);
bpavMua  = zeros(sz) + 0.002793;
bpavMusp = zeros(sz) + 1.1688;
% bpavMua  = zeros(sz) + 0.001018;
% bpavMusp = zeros(sz) + 0.628212;
% bpavMua  = zeros(sz) + 3.26464;
% bpavMusp = zeros(sz) + 6.413;

op_fit_maps4(:,:,:,1) = op_fit_maps3(:,:,:,1) - bpavMua(:,:,:,1);
op_fit_maps4(:,:,:,2) = op_fit_maps3(:,:,:,2) - bpavMusp(:,:,:,1);

%% Absolute Btm Lyr OP Error by Top Lyr OP Shift
% Skin Musp:        6.413
% Light Skin Mua:   0.66766
% Tan Skin Mua:     3.26464
% Dark Skin Mua:    9.10785

% LUT3.MuaD = LUT1.Mua(:,:) - 0.66766;
% LUT3.MuaD = LUT1.Mua(:,:) - 3.26464;
% LUT3.MuspD = LUT1.Musp(:,:) - 6.413;

LUT3.MuaD = LUT1.Mua(:,:) - 0.002793; % DOSI
LUT3.MuspD = LUT1.Musp(:,:) - 1.1688;

[xq3, yq3] = meshgrid(-0.1:0.1:9.8, -1:0.2:19);

% for i =1:length(wv)
%     op_fit_maps5(:,:,i,1) = griddata(LUT3.MuaD(:), LUT3.MuspD(:),reshape(op_fit_maps4(:,:,i,1),[],1),xq3,yq3);
%     op_fit_maps5(:,:,i,2) = griddata(LUT3.MuaD(:), LUT3.MuspD(:),reshape(op_fit_maps4(:,:,i,2),[],1),xq3,yq3);
% end

abssub1 = linspace(-0.06,0.18,9);
abssub2 = linspace(-0.05,0.2,41);

scasub1 = linspace(-2,5,15);
scasub2 = linspace(-2,5,71);

% for i =1:length(wv)
%     op_fit_maps5(:,:,i,1) = griddata(LUT3.MuaD(:), LUT3.MuspD(:),LUT2.Mua(:),xq3,yq3);
%     op_fit_maps5(:,:,i,2) = griddata(LUT3.MuaD(:), LUT3.MuspD(:),LUT2.Musp(:),xq3,yq3);
% end

for i=1:length(wv)
    f = figure(g);
    f.Position = [10,10,1200,800];
    t = tiledlayout(1,2);

    nexttile
    hold on
    contourf(LUT3.MuaD(1,:),LUT3.MuspD(:,1),op_fit_maps4(:,:,i,1),abssub2,'LineStyle','none');
    contour(LUT3.MuaD(1,:),LUT3.MuspD(:,1),op_fit_maps4(:,:,i,1),abssub1,'black','ShowText','on');
    caxis([-0.05,0.2]);
    colorbar
    hold off
    title('Absolute error in returned bottom layer \mu_a (mm^{-1})')

    nexttile
    hold on
    [m,n] = contourf(LUT3.MuaD(1,:),LUT3.MuspD(:,1),op_fit_maps4(:,:,i,2),scasub2,'LineStyle','none');
    [a,b] = contour(LUT3.MuaD(1,:),LUT3.MuspD(:,1),op_fit_maps4(:,:,i,2),scasub1,'black','ShowText','on');
    colorbar
    caxis([-2,5]);
    hold off
    title('Absolute error in returned bottom layer \mu_s'' (mm^{-1})')
    
    t.TileSpacing = 'compact';
    t.Padding = 'compact';
    title(t,'Error Propagation From Inaccurate Skin OP Assumptions')
    xlabel(t,'Top Layer Mua Change')
    ylabel(t,'Top layer Musp Change')
    g = g + 1;
end


%% Absolute Btm Lyr OP Error by Rd Vals

absrel1 = linspace(-.06,0.18,9);
absrel2 = linspace(-.05,0.2,51);

scarel1 = linspace(-2,5,15);
scarel2 = linspace(-2,5,71);

for i=1:length(wv)
    f = figure(g);
    f.Position = [10,10,1200,800];
    t = tiledlayout(1,2);

    nexttile
    hold on
    contourf(LUT.M1,LUT.M2,op_fit_maps4(:,:,i,1),absrel2,'LineStyle','none');
    contour(LUT.M1,LUT.M2,op_fit_maps4(:,:,i,1),absrel1,'black','ShowText','on');
    caxis([-0.05,0.2]);
    colorbar
    hold off
    title('Absolute error in returned bottom layer \mu_a (mm^{-1})')

    nexttile
    hold on
    [m,n] = contourf(LUT.M1,LUT.M2,op_fit_maps4(:,:,i,2),scarel2,'LineStyle','none');
    [a,b] = contour(LUT.M1,LUT.M2,op_fit_maps4(:,:,i,2),scarel1,'black','ShowText','on');
    colorbar
    caxis([-2,5]);
    hold off
    title('Absolute error in returned bottom layer \mu_s'' (mm^{-1})')
    
    t.TileSpacing = 'compact';
    t.Padding = 'compact';
    title(t,'Extracted OP error vs. R_d value')
    xlabel(t,'R_d (f_x = 0 mm^{-1})')
    ylabel(t,'R_d (f_x = 0.1 mm^{-1})')
    g = g + 1;
end

%% Relative Btm Lyr OP Error by Top Lyr OP Shift

rabsrel1 = linspace(0,70,8);
rabsrel2 = linspace(0,70,36);

rscarel1 = linspace(-1,5,13);
rscarel2 = linspace(-1,5,61);

op_fit_maps5(:,:,:,1) = op_fit_maps4(:,:,:,1) ./ bpavMua(:,:,:,1);
op_fit_maps5(:,:,:,2) = op_fit_maps4(:,:,:,2) ./ bpavMusp(:,:,:,1);

for i=1:length(wv)
    f = figure(g);
    f.Position = [10,10,1200,800];
    t = tiledlayout(1,2);

    nexttile
    hold on
    contourf(LUT3.MuaD(1,:),LUT3.MuspD(:,1),op_fit_maps5(:,:,i,1),rabsrel2,'LineStyle','none');
    contour(LUT3.MuaD(1,:),LUT3.MuspD(:,1),op_fit_maps5(:,:,i,1),rabsrel1,'black','ShowText','on');
    caxis([0,72]);
    colorbar
    hold off
    title('Relative error in returned bottom layer \mu_a (mm^{-1})')

    nexttile
    hold on
    [m,n] = contourf(LUT3.MuaD(1,:),LUT3.MuspD(:,1),op_fit_maps5(:,:,i,2),rscarel2,'LineStyle','none');
    [a,b] = contour(LUT3.MuaD(1,:),LUT3.MuspD(:,1),op_fit_maps5(:,:,i,2),rscarel1,'black','ShowText','on');
    colorbar
    caxis([-1,5]);
    hold off
    title('Relative error in returned bottom layer \mu_s'' (mm^{-1})')
    
    t.TileSpacing = 'compact';
    t.Padding = 'compact';
    title(t,'Error Propagation From Inaccurate Skin OP Assumptions')
    xlabel(t,'Top Layer Mua Change')
    ylabel(t,'Top layer Musp Change')
    g = g + 1;
end