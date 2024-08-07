close all
clear all

import numpy.*
import matplotlib.pyplot.*
import matplotlib.cm.ScalarMappable.*

%LUT1 = load('C:\Users\daemo\Documents\MATLAB\Generated_LUTs\my_LUT_8032022.mat');
%LUT2 = load('C:\Users\daemo\Documents\MATLAB\Generated_LUTs\my_LUT_7272022.mat');
%LUT1 = load('C:\Users\daemo\Documents\MATLAB\Generated_LUTs\my_LUT_Fourier_multilayer_toplayval_[03;024375]_20x20.mat');
%LUT2 = load('C:\Users\daemo\Documents\MATLAB\Generated_LUTs\my_LUT_Fourier_multilayer_toplayval_[0015;024375]_20x20.mat');
LUT1 = load('C:\Users\daemo\Documents\MATLAB\Generated_LUTs\my_LUT_Fourier_SAC_multilayer_dark_100x100.mat');
LUT2 = load('C:\Users\daemo\Documents\MATLAB\Generated_LUTs\my_LUT_Fourier_SAC_multilayer_light_100x100.mat');
%LUT2 = load('C:\Users\daemo\Documents\MATLAB\Generated_LUTs\my_LUT_Gardner_SAC_multilayer_dark_100x100.mat');

contset = linspace(0,0.9,10);
contset2 = linspace(0, 0.9, 46);
wv = 659;
[xq, yq] = meshgrid(0.001:0.001:0.1, 0.0007:0.0007:0.07);
[xq2, yq2] = meshgrid(0.008:0.008:0.80, 0.005:0.005:0.50);

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
title('Light skin AC')
caxis([0,0.9]);

ax2 = nexttile;
hold on
contourf(LUT.Mua, LUT.Musp, LUT.M1(:,:),contset2, 'LineStyle','none')
contour(LUT.Mua, LUT.Musp, LUT.M1(:,:),contset, 'black','ShowText','on')
hold off
title('Light skin DC')
caxis([0,0.9]);
cb = colorbar;
cb.Layout.Tile = 'east';
cb.Label.String = 'Diffuse Reflectance';
cb.Label.FontSize = 11;

ax3 = nexttile;
hold on
contourf(LUT1.Mua, LUT1.Musp, LUT1.M2(:,:),contset2./10, 'LineStyle','none')
contour(LUT1.Mua, LUT1.Musp, LUT1.M2(:,:),contset./10, 'black','ShowText','on')
hold off
title('Dark skin AC')
caxis([0,0.09]);

ax4 = nexttile;
hold on
contourf(LUT1.Mua, LUT1.Musp, LUT1.M1(:,:),contset2./10, 'LineStyle','none')
contour(LUT1.Mua, LUT1.Musp, LUT1.M1(:,:),contset./10, 'black','ShowText','on')
hold off
title('Dark skin DC')
caxis([0,0.09]);

cb = colorbar;
cb.Layout.Tile = 'east';
cb.Label.String = 'Diffuse Reflectance';
cb.Label.FontSize = 11;

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
    
    LUT = LUT1; % Dark
    op_fit_maps1(:,:,i,1)=griddata(LUT.M1(:), LUT.M2(:), LUT.Mua(:), xq, yq, 'linear');
    op_fit_maps1(:,:,i,2)=griddata(LUT.M1(:), LUT.M2(:), LUT.Musp(:), xq, yq,'linear');

    LUT = LUT2; % Light
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
    cb = colorbar;
    cb.Label.String = '\mu_a (mm^{-1})';
    cb.Label.FontSize = 12;
    hold off
    title('Light \mu_a (mm^{-1})')

    nexttile
    hold on
    [m,n] = contourf(xq2(1,:),yq2(:,1),op_fit_maps2(:,:,i,2),scaset2,'LineStyle','none');
    [a,b] = contour(xq2(1,:),yq2(:,1),op_fit_maps2(:,:,i,2),scaset1,'black','ShowText','on');
    cb = colorbar;
    cb.Label.String = '\mu_s'' (mm^{-1})';
    cb.Label.FontSize = 12;
    caxis([0,5]);
    hold off
    title('Light \mu_s'' (mm^{-1})')
    
    nexttile
    hold on
    contourf(xq(1,:),yq(:,1),op_fit_maps1(:,:,i,1),absset2,'LineStyle','none');
    contour(xq(1,:),yq(:,1),op_fit_maps1(:,:,i,1),absset1,'black','ShowText','on');
    caxis([0,0.2]);
    cb = colorbar;
    cb.Label.String = '\mu_a (mm^{-1})';
    cb.Label.FontSize = 12;
    hold off
    title('Dark \mu_a (mm^{-1})')

    nexttile
    hold on
    contourf(xq(1,:),yq(:,1),op_fit_maps1(:,:,i,2),scaset2,'LineStyle','none');
    contour(xq(1,:),yq(:,1),op_fit_maps1(:,:,i,2),scaset1,'black','ShowText','on');
    caxis([0,5]);
    cb = colorbar;
    cb.Label.String = '\mu_s'' (mm^{-1})';
    cb.Label.FontSize = 12;
    hold off
    title('Dark \mu_s'' (mm^{-1})')
    
    t.TileSpacing = 'compact';
    t.Padding = 'compact';
    title(t,'Diffuse Reflectance vs. Absorbance and Scattering')
    xlabel(t,'R_d (f_x = 0 mm^{-1})')
    ylabel(t,'R_d (f_x = 0.1 mm^{-1})')
    
    g = g + 1;
end


%% Compared LUT Inversions
op_fit_maps3 = op_fit_maps2(:,:,:,:) - op_fit_maps1(:,:,:,:);

op_fit_maps4 = op_fit_maps2(:,:,:,:)./op_fit_maps1(:,:,:,:);

abssub1 = linspace(0.0,0.2,11);
absrel1 = linspace(0,5,11);
scasub1 = linspace(0,1,11);
scarel1 = linspace(1,2,11);

abssub2 = linspace(0.0,0.2,51);
absrel2 = linspace(0,5,51);
scasub2 = linspace(0,1,51);
scarel2 = linspace(1,2,51);

for i=1:length(wv)
    f = figure(g);
    f.Position = [10,10,1200,800];
    t = tiledlayout(2,2);
    
    nexttile
    hold on
    contourf(xq(1,:),yq(:,1),op_fit_maps3(:,:,i,1),abssub2,'LineStyle','none');
    contour(xq(1,:),yq(:,1),op_fit_maps3(:,:,i,1),abssub1,'black','ShowText','on');
    caxis([0,0.2]);
    colorbar
    hold off
    title('Absolute difference in \mu_a (mm{-1})')

    nexttile
    hold on
    contourf(xq(1,:),yq(:,1),op_fit_maps3(:,:,i,2),scasub2,'LineStyle','none');
    contour(xq(1,:),yq(:,1),op_fit_maps3(:,:,i,2),scasub1,'black','ShowText','on');
    caxis([0,1]);
    colorbar
    hold off
    title('Absolute difference in \mu_s'' (mm{-1})')
    
    nexttile
    hold on
    contourf(xq(1,:),yq(:,1),op_fit_maps4(:,:,i,1),absrel2,'LineStyle','none');
    contour(xq(1,:),yq(:,1),op_fit_maps4(:,:,i,1),absrel1,'black','ShowText','on');
    caxis([1,5]);
    colorbar
    hold off
    title('Relative change in \mu_a (mm{-1})')

    nexttile
    hold on
    [m,n] = contourf(xq(1,:),yq(:,1),op_fit_maps4(:,:,i,2),scarel2,'LineStyle','none');
    [a,b] = contour(xq(1,:),yq(:,1),op_fit_maps4(:,:,i,2),scarel1,'black','ShowText','on');
    colorbar
    caxis([1,2]);
    hold off
    title('Relative change in \mu_s'' (mm{-1})')
    
    t.TileSpacing = 'compact';
    t.Padding = 'compact';
    title(t,'Changes in R_d vs. Absorbance and Scattering')
    xlabel(t,'R_d (f_x = 0 mm^{-1})')
    ylabel(t,'R_d (f_x = 0.1 mm^{-1})')
    g = g + 1;
end

%% How different are the patterns?
f = figure(g);
t = tiledlayout(2,2);

nexttile
hold on
contourf(LUT.Mua, LUT.Musp, LUT.M2(:,:)./LUT1.M2(:,:),'LineStyle','none')
contour(LUT.Mua, LUT.Musp, LUT.M2(:,:)./LUT1.M2(:,:),'black','ShowText','on');
colorbar
hold off
title('Relative Rd values for 0.1mm^{-1}')

nexttile
hold on
contourf(LUT.Mua, LUT.Musp, LUT.M1(:,:)./LUT1.M1(:,:),'LineStyle','none')
contour(LUT.Mua, LUT.Musp, LUT.M1(:,:)./LUT1.M1(:,:),'black','ShowText','on');
colorbar
hold off
title('Relative Rd values for 0mm^{-1}')

S1 = std(reshape(LUT.M2./LUT1.M2,1,[]));
M1 = mean(reshape(LUT.M2./LUT1.M2,1,[]));
S0 = std(reshape(LUT.M1./LUT1.M1,1,[]));
M0 = mean(reshape(LUT.M1./LUT1.M1,1,[]));

nexttile
hold on
contourf(LUT.Mua, LUT.Musp, (LUT.M2(:,:)./LUT1.M2(:,:)-M1)/S1,'LineStyle','none')
contour(LUT.Mua, LUT.Musp, (LUT.M2(:,:)./LUT1.M2(:,:)-M1)/S1,'black','ShowText','on');
colorbar
hold off
title('Relative Rd value standard deviations for 0.1mm^{-1}')

nexttile
hold on
contourf(LUT.Mua, LUT.Musp, (LUT.M1(:,:)./LUT1.M1(:,:)-M0)/S0,'LineStyle','none')
contour(LUT.Mua, LUT.Musp, (LUT.M1(:,:)./LUT1.M1(:,:)-M0)/S0,'black','ShowText','on');
colorbar
hold off
title('Relative Rd value standard deviations for 0mm^{-1}')

t.TileSpacing = 'compact';
t.Padding = 'compact';
title(t,'Relative Diffuse Reflectance values from simulation Light skin vs Dark skin')
xlabel(t,'\mu_a (mm^{-1})')
ylabel(t,'\mu_s'' (mm^{-1})')

%% Axis Shift - LUT inversions

absset1 = linspace(0,0.2,11);
absset2 = linspace(0,0.2,51);
scaset1 = linspace(0,5,11);
scaset2 = linspace(0,5,51);

clear op_fit_maps

g = 5;
for i=1:length(wv)

    disp(['Processing ' num2str(wv) ' nm...'])
    
    LUT = LUT1; % Dark
    op_fit_maps1(:,:,i,1)=griddata(LUT.M1(:), LUT.M2(:), LUT.Mua(:), xq2, yq2, 'linear');
    op_fit_maps1(:,:,i,2)=griddata(LUT.M1(:), LUT.M2(:), LUT.Musp(:), xq2, yq2,'linear');

    LUT = LUT2; % Light
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
    cb = colorbar;
    cb.Label.String = '\mu_a (mm^{-1})';
    cb.Label.FontSize = 12;
    hold off
    title('Light \mu_a (mm^{-1})')

    nexttile
    hold on
    [m,n] = contourf(xq2(1,:),yq2(:,1),op_fit_maps2(:,:,i,2),scaset2,'LineStyle','none');
    [a,b] = contour(xq2(1,:),yq2(:,1),op_fit_maps2(:,:,i,2),scaset1,'black','ShowText','on');
    cb = colorbar;
    cb.Label.String = '\mu_s'' (mm^{-1})';
    cb.Label.FontSize = 12;
    caxis([0,5]);
    hold off
    title('Light \mu_s'' (mm^{-1})')
    
    nexttile
    hold on
    contourf(xq2(1,:),yq2(:,1),op_fit_maps1(:,:,i,1),absset2,'LineStyle','none');
    contour(xq2(1,:),yq2(:,1),op_fit_maps1(:,:,i,1),absset1,'black','ShowText','on');
    caxis([0,0.2]);
    cb = colorbar;
    cb.Label.String = '\mu_a (mm^{-1})';
    cb.Label.FontSize = 12;
    hold off
    title('Dark \mu_a (mm^{-1})')

    nexttile
    hold on
    contourf(xq2(1,:),yq2(:,1),op_fit_maps1(:,:,i,2),scaset2,'LineStyle','none');
    contour(xq2(1,:),yq2(:,1),op_fit_maps1(:,:,i,2),scaset1,'black','ShowText','on');
    caxis([0,5]);
    cb = colorbar;
    cb.Label.String = '\mu_s'' (mm^{-1})';
    cb.Label.FontSize = 12;
    hold off
    title('Dark \mu_s'' (mm^{-1})')
    
    t.TileSpacing = 'compact';
    t.Padding = 'compact';
    title(t,'Diffuse Reflectance vs. Absorbance and Scattering')
    xlabel(t,'R_d (f_x = 0 mm^{-1})')
    ylabel(t,'R_d (f_x = 0.1 mm^{-1})')
    
    g = g + 1;
end