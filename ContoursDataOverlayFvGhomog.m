close all
clear all

import numpy.*
import matplotlib.pyplot.*
import matplotlib.cm.ScalarMappable.*

%LUT1 = load('C:\Users\daemo\Documents\MATLAB\Generated_LUTs\my_LUT_8032022.mat');
%LUT2 = load('C:\Users\daemo\Documents\MATLAB\Generated_LUTs\my_LUT_7272022.mat');
%LUT1 = load('C:\Users\daemo\Documents\MATLAB\Generated_LUTs\my_LUT_Fourier_multilayer_toplayval_[03;024375]_20x20.mat');
%LUT2 = load('C:\Users\daemo\Documents\MATLAB\Generated_LUTs\my_LUT_Fourier_multilayer_toplayval_[0015;024375]_20x20.mat');
%LUT1 = load('C:\Users\daemo\Documents\MATLAB\Generated_LUTs\my_LUT_Fourier_SAC_multilayer_dark_100x100.mat');
%LUT2 = load('C:\Users\daemo\Documents\MATLAB\Generated_LUTs\my_LUT_Gardner_SAC_multilayer_light_100x100.mat');
%LUT1 = load('C:\Users\daemo\Documents\MATLAB\Generated_LUTs\my_LUT_Fourier_SAC_multilayer_light_100x100.mat');
%LUT2 = load('C:\Users\daemo\Documents\MATLAB\Generated_LUTs\my_LUT_Gardner_SAC_multilayer_dark_100x100.mat');
%LUT2 = load('C:\Users\daemo\Documents\MATLAB\Generated_LUTs\my_LUT_Fourier_SAC_multilayer_tan_100x100.mat');
LUT2 = load('C:\Users\daemo\Documents\MATLAB\Generated_LUTs\my_LUT_Fourier_homogenous_100x100.mat');
LUT1 = load('C:\Users\daemo\Documents\MATLAB\SFDI Analysis Code 5-11-20\LUTs\Gardner_homogeneous_Large_Fx.mat');

contset = linspace(0,0.9,10);
contset2 = linspace(0, 0.9, 46);
wv = 659;
[xq, yq] = meshgrid(0.01:0.01:1.00, 0.009:0.009:0.90);

CarGardmcxMua = [0.0040,0.0015,0.0245,0.0083,0.1064];
CarGardmcxMusp = [3.4547,3.9897,2.3753,2.4239,1.9394];
CarGardhomogMua = [0.005855473,0.023570119,0.060782356,0.036852787,0.108724719,0.001053252];
CarGardhomogMusp = [0.910441792,0.85052589,0.681855078,0.625509067,0.614984873,0.537075917];
CarFourhomogMua = [0.017625017,0.071274647,0.167294103,0.104784842,0.199813919,0.001122827];
CarFourhomogMusp = [1.072762156,1.024691874,0.819354912,0.750393706,0.64469535,0.63475942];

CarGardhomogRd0 = [0.661659273,0.454093399,0.269750823,0.333765742,0.174230169,0.785160221];
CarGardhomogRd1 = [0.382892035,0.309907368,0.198389603,0.216937922,0.137250882,0.252560959];
CarFourhomogRd0 = [0.661633277,0.454472423,0.269886738,0.334279583,0.176331452,0.784648624];
CarFourhomogRd1 = [0.382915569,0.310110184,0.198432518,0.217352956,0.138855877,0.252079626];

CarlayerLabels = ['Skin2b';'Skin3b';'Skin4b';'Skin5b';'Skin6b'];
CarhomogLabels = ['Skin2';'Skin3';'Skin4';'Skin5';'Skin6';'bpav4'];

%% Comparing LUTs

LUT = LUT1.LUT;
f = figure(1);
f.Position = [10,10,600,400];
t = tiledlayout(2,2);
ax1 = nexttile;
hold on
contourf(LUT2.Mua, LUT2.Musp, LUT2.M2(:,:),contset2, 'LineStyle','none')
contour(LUT2.Mua, LUT2.Musp, LUT2.M2(:,:),contset, 'black','ShowText','on')
scatter(CarFourhomogMua,CarFourhomogMusp,[],"red","filled")
labelpoints(CarFourhomogMua,CarFourhomogMusp,CarhomogLabels,'N',0.1)
hold off
title('Fourier AC')
caxis([0,0.9]);

ax2 = nexttile;
hold on
contourf(LUT2.Mua, LUT2.Musp, LUT2.M1(:,:),contset2, 'LineStyle','none')
contour(LUT2.Mua, LUT2.Musp, LUT2.M1(:,:),contset, 'black','ShowText','on')
%scatter(CarGardMua,CarGardMusp,[],"red","filled")
scatter(CarFourhomogMua,CarFourhomogMusp,[],"red","filled")
labelpoints(CarFourhomogMua,CarFourhomogMusp,CarhomogLabels,'N',0.1)
hold off
title('Fourier DC')
caxis([0,0.9]);

ax3 = nexttile;
hold on
contourf(LUT.Mua, LUT.Musp, LUT.M2(:,:),contset2, 'LineStyle','none')
contour(LUT.Mua, LUT.Musp, LUT.M2(:,:),contset, 'black','ShowText','on')
scatter(CarGardhomogMua,CarGardhomogMusp,[],"red","filled")
labelpoints(CarGardhomogMua,CarGardhomogMusp,CarhomogLabels,'N',0.1)
hold off
title('Gardner AC')
caxis([0,0.9]);

ax4 = nexttile;
hold on
contourf(LUT.Mua, LUT.Musp, LUT.M1(:,:),contset2, 'LineStyle','none')
contour(LUT.Mua, LUT.Musp, LUT.M1(:,:),contset, 'black','ShowText','on')
scatter(CarGardhomogMua,CarGardhomogMusp,[],"red","filled")
labelpoints(CarGardhomogMua,CarGardhomogMusp,CarhomogLabels,'N',0.1)
hold off
title('Gardner DC')
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
scaset1 = linspace(0,5,11);
scaset2 = linspace(0,5,51);

clear op_fit_maps

g = 2;
for i=1:length(wv)

    disp(['Processing ' num2str(wv) ' nm...'])
    
    LUT = LUT1.LUT;
    op_fit_maps1(:,:,i,1)=griddata(LUT.M1(:), LUT.M2(:), LUT.Mua(:), xq, yq, 'linear');
    op_fit_maps1(:,:,i,2)=griddata(LUT.M1(:), LUT.M2(:), LUT.Musp(:), xq, yq,'linear');

    LUT = LUT2;
    op_fit_maps2(:,:,i,1)=griddata(LUT.M1(:), LUT.M2(:), LUT.Mua(:), xq, yq, 'linear');
    op_fit_maps2(:,:,i,2)=griddata(LUT.M1(:), LUT.M2(:), LUT.Musp(:), xq, yq,'linear');

    f = figure(g);
    f.Position = [10,10,900,600];
    t = tiledlayout(2,2);

    nexttile
    hold on
    contourf(xq(1,:),yq(:,1),op_fit_maps1(:,:,i,1),absset2,'LineStyle','none');
    contour(xq(1,:),yq(:,1),op_fit_maps1(:,:,i,1),absset1,'black','ShowText','on');
    scatter(CarGardhomogRd0,CarGardhomogRd1,[],"red","filled")
    labelpoints(CarGardhomogRd0,CarGardhomogRd1,CarhomogLabels,'N',0.1)
    caxis([0,0.2]);
    colorbar
    hold off
    title('Gardner \mu_a (mm{-1})')


    nexttile
    hold on
    contourf(xq(1,:),yq(:,1),op_fit_maps1(:,:,i,2),scaset2,'LineStyle','none');
    contour(xq(1,:),yq(:,1),op_fit_maps1(:,:,i,2),scaset1,'black','ShowText','on');
    scatter(CarGardhomogRd0,CarGardhomogRd1,[],"red","filled")
    labelpoints(CarGardhomogRd0,CarGardhomogRd1,CarhomogLabels,'N',0.1)
    caxis([0,5]);
    colorbar
    hold off
    title('Gardner \mu_s'' (mm{-1})')
    
    nexttile
    hold on
    contourf(xq(1,:),yq(:,1),op_fit_maps2(:,:,i,1),absset2,'LineStyle','none');
    contour(xq(1,:),yq(:,1),op_fit_maps2(:,:,i,1),absset1,'black','ShowText','on');
    scatter(CarFourhomogRd0,CarFourhomogRd1,[],"red","filled")
    labelpoints(CarFourhomogRd0,CarFourhomogRd1,CarhomogLabels,'N',0.1)
    caxis([0,0.2]);
    colorbar
    hold off
    title('Fourier \mu_a (mm{-1})')

    nexttile
    hold on
    [m,n] = contourf(xq(1,:),yq(:,1),op_fit_maps2(:,:,i,2),scaset2,'LineStyle','none');
    [a,b] = contour(xq(1,:),yq(:,1),op_fit_maps2(:,:,i,2),scaset1,'black','ShowText','on');
    scatter(CarFourhomogRd0,CarFourhomogRd1,[],"red","filled")
    labelpoints(CarFourhomogRd0,CarFourhomogRd1,CarhomogLabels,'N',0.1)
    colorbar
    caxis([0,5]);
    hold off
    title('Fourier \mu_s'' (mm{-1})')
    
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
% 
% abssub1 = linspace(-0.15,0.05,11);
% absrel1 = linspace(0,3,11);
% scasub1 = linspace(-0.5,3,15);
% scarel1 = linspace(0,2,11);
% 
% abssub2 = linspace(-0.15,0.05,51);
% absrel2 = linspace(0,3,51);
% scasub2 = linspace(-0.5,3,71);
% scarel2 = linspace(0,2,51);
% 
% for i=1:length(wv)    
%     f = figure(g);
%     f.Position = [10,10,1200,800];
%     t = tiledlayout(2,2);
%     
%     nexttile
%     hold on
%     contourf(xq(1,:),yq(:,1),op_fit_maps3(:,:,i,1),abssub2,'LineStyle','none');
%     contour(xq(1,:),yq(:,1),op_fit_maps3(:,:,i,1),abssub1,'black','ShowText','on');
%     scatter(CarGardmcxRd0,CarGardmcxRd1,[],"red","filled")
%     labelpoints(CarGardmcxRd0,CarGardmcxRd1,CarlayerLabels,'N',0.1)
%     caxis([-0.05,0.05]);
%     colorbar
%     hold off
%     title('Absolute difference in \mu_a (mm{-1})')
% 
%     nexttile
%     hold on
%     contourf(xq(1,:),yq(:,1),op_fit_maps3(:,:,i,2),scasub2,'LineStyle','none');
%     contour(xq(1,:),yq(:,1),op_fit_maps3(:,:,i,2),scasub1,'black','ShowText','on');
%     scatter(CarGardmcxRd0,CarGardmcxRd1,[],"red","filled")
%     labelpoints(CarGardmcxRd0,CarGardmcxRd1,CarlayerLabels,'N',0.1)
%     caxis([-0.5,3]);
%     colorbar
%     hold off
%     title('Absolute difference in \mu_s'' (mm{-1})')
%     
%     nexttile
%     hold on
%     contourf(xq(1,:),yq(:,1),op_fit_maps4(:,:,i,1),absrel2,'LineStyle','none');
%     contour(xq(1,:),yq(:,1),op_fit_maps4(:,:,i,1),absrel1,'black','ShowText','on');
%     scatter(CarGardmcxRd0,CarGardmcxRd1,[],"red","filled")
%     labelpoints(CarGardmcxRd0,CarGardmcxRd1,CarlayerLabels,'N',0.1)
%     caxis([0,3]);
%     colorbar
%     hold off
%     title('Relative change in \mu_a (mm{-1})')
% 
%     nexttile
%     hold on
%     [m,n] = contourf(xq(1,:),yq(:,1),op_fit_maps4(:,:,i,2),scarel2,'LineStyle','none');
%     [a,b] = contour(xq(1,:),yq(:,1),op_fit_maps4(:,:,i,2),scarel1,'black','ShowText','on');
%     scatter(CarGardmcxRd0,CarGardmcxRd1,[],"red","filled")
%     labelpoints(CarGardmcxRd0,CarGardmcxRd1,CarlayerLabels,'N',0.1)
%     colorbar
%     caxis([0,2]);
%     hold off
%     title('Relative change in \mu_s'' (mm{-1})')
%     
%     t.TileSpacing = 'compact';
%     t.Padding = 'compact';
%     title(t,'Changes in R_d vs. Absorbance and Scattering')
%     xlabel(t,'R_d (f_x = 0 mm^{-1})')
%     ylabel(t,'R_d (f_x = 0.1 mm^{-1})')
%     g = g + 1;
% end