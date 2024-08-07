close all
clear all

import numpy.*
import matplotlib.pyplot.*
import matplotlib.cm.ScalarMappable.*

LUT1 = load('C:\Users\daemo\Documents\MATLAB\Generated_LUTs\my_LUT_Fourier_homogenous_100x100.mat');
LUT2 = load('C:\Users\daemo\Documents\MATLAB\Generated_LUTs\my_LUT_Gardner_MCX_homogenous_100x100.mat');
LUT3 = load('C:\Users\daemo\Documents\MATLAB\SFDI Analysis Code 5-11-20\LUTs\LUT_Homogeneous_Large_1e6.mat');

contset = linspace(0,0.1,11);
contset2 = linspace(0, 0.1, 51);
wv = 659;
xq = ;
yq = ;

%% Comparing LUTs
% 
% LUT = LUT2.LUT;
% f = figure(1);
% f.Position = [10,10,600,400];
% t = tiledlayout(2,2);
% ax1 = nexttile;
% hold on
% contourf(LUT.Mua, LUT.Musp, LUT.M2(:,:),contset2, 'LineStyle','none')
% contour(LUT.Mua, LUT.Musp, LUT.M2(:,:),contset, 'black','ShowText','on')
% hold off
% title('Gardner AC')
% caxis([0,0.09]);
% 
% ax2 = nexttile;
% hold on
% contourf(LUT.Mua, LUT.Musp, LUT.M1(:,:),contset2, 'LineStyle','none')
% contour(LUT.Mua, LUT.Musp, LUT.M1(:,:),contset, 'black','ShowText','on')
% hold off
% title('Gardner DC')
% caxis([0,0.09]);
% 
% ax3 = nexttile;
% hold on
% contourf(LUT1.Mua, LUT1.Musp, LUT1.M2(:,:),contset2, 'LineStyle','none')
% contour(LUT1.Mua, LUT1.Musp, LUT1.M2(:,:),contset, 'black','ShowText','on')
% hold off
% title('Fourier AC')
% caxis([0,0.09]);
% 
% ax4 = nexttile;
% hold on
% contourf(LUT1.Mua, LUT1.Musp, LUT1.M1(:,:),contset2, 'LineStyle','none')
% contour(LUT1.Mua, LUT1.Musp, LUT1.M1(:,:),contset, 'black','ShowText','on')
% hold off
% title('Fourier DC')
% caxis([0,0.09]);
% 
% cb = colorbar;
% cb.Layout.Tile = 'east';
% 
% t.TileSpacing = 'compact';
% t.Padding = 'compact';
% title(t,'Diffuse Reflectance vs. Absorbance and Scattering')
% xlabel(t,'\mu_a (mm^{-1})')
% ylabel(t,'\mu_s'' (mm^{-1})')
% 
% 
% 

%% Basic LUT inversions

absset1 = linspace(0,0.2,11);
absset2 = linspace(0,0.2,51);
scaset1 = linspace(0,5,11);
scaset2 = linspace(0,5,51);

clear op_fit_maps

g = 2;
for i=1:length(wv)

    disp(['Processing ' num2str(wv) ' nm...'])
    
    LUT = LUT1;
    op_fit_maps1(:,:,i,1)=griddata(LUT.M1(:), LUT.M2(:), LUT.Mua(:), xq, yq, 'linear');
    op_fit_maps1(:,:,i,2)=griddata(LUT.M1(:), LUT.M2(:), LUT.Musp(:), xq, yq,'linear');

    LUT = LUT2.LUT;
    op_fit_maps2(:,:,i,1)=griddata(LUT.M1(:), LUT.M2(:), LUT.Mua(:), xq, yq, 'linear');
    op_fit_maps2(:,:,i,2)=griddata(LUT.M1(:), LUT.M2(:), LUT.Musp(:), xq, yq,'linear');

    f = figure(g);
    f.Position = [10,10,900,600];
    t = tiledlayout(2,2);

    nexttile
    hold on
    contourf(xq(1,:),yq(:,1),op_fit_maps1(:,:,i,1),absset2,'LineStyle','none');
    contour(xq(1,:),yq(:,1),op_fit_maps1(:,:,i,1),absset1,'black','ShowText','on');
    caxis([0,0.2]);
    colorbar
    hold off
    title('Fourier \mu_a (mm{-1})')


    nexttile
    hold on
    contourf(xq(1,:),yq(:,1),op_fit_maps1(:,:,i,2),scaset2,'LineStyle','none');
    contour(xq(1,:),yq(:,1),op_fit_maps1(:,:,i,2),scaset1,'black','ShowText','on');
    caxis([0,5]);
    colorbar
    hold off
    title('Fourier \mu_s'' (mm{-1})')
    
    nexttile
    hold on
    contourf(xq(1,:),yq(:,1),op_fit_maps2(:,:,i,1),absset2,'LineStyle','none');
    contour(xq(1,:),yq(:,1),op_fit_maps2(:,:,i,1),absset1,'black','ShowText','on');
    caxis([0,0.2]);
    colorbar
    hold off
    title('Gardner \mu_a (mm{-1})')

    nexttile
    hold on
    [m,n] = contourf(xq(1,:),yq(:,1),op_fit_maps2(:,:,i,2),scaset2,'LineStyle','none');
    [a,b] = contour(xq(1,:),yq(:,1),op_fit_maps2(:,:,i,2),scaset1,'black','ShowText','on');
    colorbar
    caxis([0,5]);
    hold off
    title('Gardner \mu_s'' (mm{-1})')
    
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

abssub1 = linspace(-0.12,0.02,8);
absrel1 = linspace(0,2.5,11);
scasub1 = linspace(-0.5,2.5,13);
scarel1 = linspace(0,3,13);

abssub2 = linspace(-0.12,0.02,36);
absrel2 = linspace(0,2.5,51);
scasub2 = linspace(-0.5,2.5,61);
scarel2 = linspace(0,3,61);

for i=1:length(wv)    
    f = figure(g);
    f.Position = [10,10,1200,800];
    t = tiledlayout(2,2);
    
    nexttile
    hold on
    contourf(xq(1,:),yq(:,1),op_fit_maps3(:,:,i,1),abssub2,'LineStyle','none');
    contour(xq(1,:),yq(:,1),op_fit_maps3(:,:,i,1),abssub1,'black','ShowText','on');
    caxis([-0.15,0.05]);
    colorbar
    hold off
    title('Absolute difference in \mu_a (mm{-1})')

    nexttile
    hold on
    contourf(xq(1,:),yq(:,1),op_fit_maps3(:,:,i,2),scasub2,'LineStyle','none');
    contour(xq(1,:),yq(:,1),op_fit_maps3(:,:,i,2),scasub1,'black','ShowText','on');
    caxis([-0.5,2.5]);
    colorbar
    hold off
    title('Absolute difference in \mu_s'' (mm{-1})')
    
    nexttile
    hold on
    contourf(xq(1,:),yq(:,1),op_fit_maps4(:,:,i,1),absrel2,'LineStyle','none');
    contour(xq(1,:),yq(:,1),op_fit_maps4(:,:,i,1),absrel1,'black','ShowText','on');
    caxis([0,2.5]);
    colorbar
    hold off
    title('Relative change in \mu_a (mm{-1})')

    nexttile
    hold on
    [m,n] = contourf(xq(1,:),yq(:,1),op_fit_maps4(:,:,i,2),scarel2,'LineStyle','none');
    [a,b] = contour(xq(1,:),yq(:,1),op_fit_maps4(:,:,i,2),scarel1,'black','ShowText','on');
    colorbar
    caxis([0,3]);
    hold off
    title('Relative change in \mu_s'' (mm{-1})')
    
    t.TileSpacing = 'compact';
    t.Padding = 'compact';
    title(t,'Changes in R_d vs. Absorbance and Scattering')
    xlabel(t,'R_d (f_x = 0 mm^{-1})')
    ylabel(t,'R_d (f_x = 0.1 mm^{-1})')
    g = g + 1;
end