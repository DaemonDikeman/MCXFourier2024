%% Values to Graph
clear all

freq = [0, 0.05/3, 0.10/3, 0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5];

curfold = pwd;

e = dir(curfold +"\\SFDI_Fourier_LUT_0*.csv");
names = {e.name};
for i=1:length(e)
    L(:,:,i) = readmatrix(string(names(i)));
end

ee = dir(curfold +"\\SFDI_Fourier_LUT_2*.csv");
names = {ee.name};
for i=1:length(ee)
    L1(:,:,i) = readmatrix(string(names(i)));
end

%% Heatmap
close all

figure(1)
A = L(:,:,1);

h = heatmap(0.005:0.005:0.05, 0.5:0.5:5, log10(A), 'CellLabelColor', 'none');

h.Title = 'SFDI Fourier Look-Up Table Freq 0 Log';
h.XLabel = 'Absorbance';
h.YLabel = 'Scattering';
h.GridVisible = 'off';

figure(2)
B = L(:,:,2);

g = heatmap(0.005:0.005:0.05, 0.5:0.5:5, log10(B), 'CellLabelColor', 'none');

g.Title = 'SFDI Fourier Look-Up Table Freq 0.1 Log';
g.XLabel = 'Absorbance';
g.YLabel = 'Scattering';
g.GridVisible = 'off';

figure(3)
C = L1(:,:,1);

j = heatmap(0.00:0.005:0.2, 0.15:0.15:3, log10(C), 'CellLabelColor', 'none');

j.Title = 'SFDI Fourier Look-Up Table Multilayer Freq 0 Log';
j.XLabel = 'Absorbance';
j.YLabel = 'Scattering';
j.GridVisible = 'off';

figure(4)
D = L1(:,:,2);

f = heatmap(0.00:0.005:0.2, 0.15:0.15:3, log10(D), 'CellLabelColor', 'none');

f.Title = 'SFDI Fourier Look-Up Table Multilauyer Freq 0.1 Log';
f.XLabel = 'Absorbance';
f.YLabel = 'Scattering';
f.GridVisible = 'off';

load('C:\Users\daemo\Documents\MATLAB\Generated_LUTs\my_LUTa.mat')

figure(5)
v = heatmap(LUT.Mua(1,:),LUT.Musp(:,1),log10(LUT.M1));

v.Title = 'SFDI Gardner Look-Up Table Freq 0 Log';
v.XLabel = 'Absorbance';
v.YLabel = 'Scattering';
v.GridVisible = 'off';

figure(6)
w = heatmap(LUT.Mua(1,:),LUT.Musp(:,1),log10(LUT.M2));

w.Title = 'SFDI Gardner Look-Up Table Freq 0.1 Log';
w.XLabel = 'Absorbance';
w.YLabel = 'Scattering';
w.GridVisible = 'off';
