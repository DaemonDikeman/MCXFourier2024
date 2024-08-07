%% Values to Graph
clear
close all

freq = [0, 0.05/3, 0.10/3, 0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5];
freq2 = [0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.10];

curfold = pwd;

a = dir(curfold +"\\SFDI_dref_nphoton_e9" + "\\SFDI_drefval_*.csv");
names = {a.name};
for i=1:length(a)
    L(i,:) = readmatrix(string(names(i)));
end

aa = dir(curfold +"\\SFDI_dref_nphoton_e7_g09" + "\\SFDI_drefval_*.csv");
names = {aa.name};
for i=1:length(aa)
    L2(i,:) = readmatrix(string(names(i)));
end

aaa = dir(curfold +"\\SFDI_dref_nphoton_e5_g09" +"\\SFDI_drefval_*.csv");
names = {aaa.name};
for i=1:length(aaa)
    L3(i,:) = readmatrix(string(names(i)));
end

b = dir(curfold + "\\SFDI_drefcalc" + "\\SFDI_drefcalc_g_0.90_dual_axis_nphoton_1.0e+07*.csv");
names = {b.name};
for i=1:length(b)
    M((i*4-3):(i*4),:) = readmatrix(string(names(i)));
end

bb = dir(curfold + "\\SFDI_drefcalc" + "\\SFDI_drefcalc_g_0.90_dual_axis_nphoton_1.0e+05*.csv");
names = {bb.name};
for i=1:length(bb)
    M2((i*4-3):(i*4),:) = readmatrix(string(names(i)));
end

c = dir(curfold + "\\SFDI_dref_nphoton_4e7_05mm_vx" + "\\SFDI_drefval_*.csv");
names = {c.name};
for i=1:length(c)
    N(i,:) = readmatrix(string(names(i)));
end

cc = dir(curfold + "\\SFDI_dref_nphoton_e7_05mm_vx" + "\\SFDI_drefval_*.csv");
names = {cc.name};
for i=1:length(cc)
    N2(i,:) = readmatrix(string(names(i)));
end

d = dir(curfold + "\\SFDI_drefcalc" + "\\SFDI_drefcalc_g_0.90_single_axis_nphoton_1.0e+07*.csv");
names = {d.name};
for i=1:length(d)
    P((i*4-3):(i*4),:) = readmatrix(string(names(i)));
end

dd = dir(curfold + "\\SFDI_drefcalc" + "\\SFDI_drefcalc_g_0.90_single_axis_nphoton_1.0e+05*.csv");
names = {dd.name};
for i=1:length(dd)
    P2((i*4-3):(i*4),:) = readmatrix(string(names(i)));
end

ddd = dir(curfold +"\\SFDI_dref_nphoton_e7_g071" + "\\SFDI_drefval_*.csv");
names = {ddd.name};
for i=1:length(ddd)
    P3(i,:) = readmatrix(string(names(i)));
end

% e = dir(curfold + "\\SFDI_focused_freq" + "\\SFDI_drefval_g_0.71_*.csv");
% names = {e.name};
% for i=1:length(e)
%     Q(i,:) = readmatrix(string(names(i)));
% end
% 
% ee = dir(curfold + "\\SFDI_focused_freq" + "\\SFDI_drefcalc_g_adjust_*.csv");
% names = {ee.name};
% for i=1:length(ee)
%     Q2((i*4-3):(i*4),:) = readmatrix(string(names(i)));
% end

Mlog = log10(M);
Mlog2 = log10(M2);
Llog = log10(L);
Llog2 = log10(L2);
Llog3 = log10(L3);
Nlog = log10(N);
Nlog2 = log10(N2);
Plog = log10(P);
Plog2 = log10(P2);
Plog3 = log10(P3);
% Qlog = log10(Q);
% Qlog2 = log10(Q2);

%% Graph
close all

wide = 2;
thin = 1;

newcolors =    [0.9290    0.6040    0.1250
                     0    0.4470    0.7410
                0.1250    0.8600    0.1250
                0.4940    0.1840    0.5560];

% ['Diffuse Reflectance vs. Spatial freq. (mm^{-1}) at diff. L* values']

%% Fourier vs Pencil

figure(1) % e5 fourier source vs. e5 pencil source
hold on
for i = 1:4
    plot(freq, Llog3(i,:), '-o', 'Linewidth', wide)
end
for i = 1:4
    plot(freq, Mlog2(i,:), '--*', 'Linewidth', thin)
end
% for i = 1:length(d)
%     plot(freq, Llog3(i,:), '-o', 'Linewidth', wide) % increase line width
%     plot(freq, Mlog(i,:), '--*', 'Linewidth', thin)
% end
colororder(newcolors)
ylim([-2.1 0])
xlim([0 0.5])
ylabel("Log10 Rd")
xlabel("Spatial freq. (mm^{-1})")
yticks([-2 -1.5 -1 -0.5 0])
yticklabels(["10^{-2}", "10^{-1.5}", "10^{-1}", "10^{-0.5}", "10^{0}"])
legend("Fourier L*=0.5","Fourier L*=1.0","Fourier L*=2.0","Fourier L*=4.0",...
       "Pencil  L*=0.5","Pencil  L*=1.0","Pencil  L*=2.0","Pencil  L*=4.0")
grid on
[t,s] = title(['Diffuse Reflectance vs. Spatial Freq. (mm^{-1}) at diff. L* values'],...
              ['1mm voxels, Fourier e5 vs. Pencil e5 Sources']);


figure(2) % not log plot; e5 fourier source vs. e5 pencil source
hold on
for i = 1:4
    plot(freq, L3(i,:), '-o', 'Linewidth', wide)
end
for i = 1:4
    plot(freq, M2(i,:), '--*', 'Linewidth', thin)
end
colororder(newcolors)
% ylim([-2 0])
% xlim([0 0.5])
ylabel("Rd")
xlabel("Spatial freq. (mm^{-1})")
% yticks([-2 -1.5 -1 -0.5 0])
% yticklabels(["10^{-2}", "10^{-1.5}", "10^{-1}", "10^{-0.5}", "10^{0}"])
legend("Fourier L*=0.5","Fourier L*=1.0","Fourier L*=2.0","Fourier L*=4.0",...
       "Pencil  L*=0.5","Pencil  L*=1.0","Pencil  L*=2.0","Pencil  L*=4.0")
grid on
[t,s] = title(['Diffuse Reflectance vs. Spatial Freq. (mm^{-1}) at diff. L* values'],...
              ['1mm voxels, Fourier e5 vs. Pencil e5 Sources']);

%% Photon count differences

figure(3) % e5 fourier source vs. e7 fourier source vs. e9 fourier source
hold on
for i = 1:4
    plot(freq, Llog3(i,:), '-o', 'Linewidth', wide)
end
for i = 1:4
    plot(freq, Llog2(i,:), '--*', 'Linewidth', thin)
end
for i = 1:4
    plot(freq, Llog(i,:), '-.', 'Linewidth', wide)
end

colororder(newcolors)
ylim([-2.1 0])
xlim([0 0.5])
x0=280;
y0=320;
width=750; % 1000
height=600; % 800
set(gcf,'position',[x0,y0,width,height])
ax = gca;
ax.FontSize = 14; % 18
ylabel("Log10 Rd",'FontSize',14) % 18
xlabel("Spatial freq. (mm^{-1})",'FontSize',14) % 18
yticks([-2 -1.5 -1 -0.5 0])
yticklabels(["10^{-2}", "10^{-1.5}", "10^{-1}", "10^{-0.5}", "10^{0}"])
% legend("F. e5 L*=0.5","F. e5 L*=1.0","F. e5 L*=2.0","F. e5 L*=4.0",...
%        "F. e7 L*=0.5","F. e7 L*=1.0","F. e7 L*=2.0","F. e7 L*=4.0",...
%        "F. e9 L*=0.5","F. e9 L*=1.0","F. e9 L*=2.0","F. e9 L*=4.0")
legend({"1e5 photons",'','','',...
       "1e7 photons",'','','',...
       "1e9 photons",'','',''},...
       'FontSize',14) % 18
grid on
hold off
% [t,s] = title(['Diffuse Reflectance vs. Spatial freq. (mm^{-1}) w/ varying L* values and photon counts'],...
%               ['1mm voxels, g=0.9 e5 v. e7 v. e9 Fourier sources']);
[t] = title(['MCX Fourier R_d vs. Spatial freq. (mm^{-1}) at diff. photon counts'],'FontSize',16); % 18

figure(4) % e5 pencil source vs. e7 pencil source
hold on
for i = 1:4
    plot(freq, Plog2(i,:), '-o', 'Linewidth', wide)
end
for i = 1:4
    plot(freq, Plog(i,:), '--*', 'Linewidth', thin)
end
colororder(newcolors)
ylim([-2.1 0])
xlim([0 0.5])
x0=1280;
y0=320;
width=750; % 1000
height=600; % 800
set(gcf,'position',[x0,y0,width,height])
ax = gca;
ax.FontSize = 14; % 18
ylabel("Log10 Rd",'FontSize',14) % 18
xlabel("Spatial freq. (mm^{-1})",'FontSize',14) % 18
yticks([-2 -1.5 -1 -0.5 0])
yticklabels(["10^{-2}", "10^{-1.5}", "10^{-1}", "10^{-0.5}", "10^{0}"])
% legend("P. e5 L*=0.5","P. e5 L*=1.0","P. e5 L*=2.0","P. e5 L*=4.0",...
%        "P. e7 L*=0.5","P. e7 L*=1.0","P. e7 L*=2.0","P. e7 L*=4.0")
legend({"1e5 photons",'','','',...
       "1e7 photons",'','',''},...
       'FontSize',14) % 'FontSize',18)
grid on
hold off
% [t,s] = title(['Diffuse Reflectance vs. Spatial freq. (mm^{-1}) w/ varying L* values and photon counts'],...
%               ['1mm voxels, g=0.9 Pencil e5 vs. Pencil e7 Sources']);
[t] = title(['MCX Gardner R_d vs. Spatial freq. (mm^{-1}) at diff. photon counts'],'FontSize',16); % 18

figure(5) % e7 & e9 fourier source vs. e7 single and dual pencil source (gardner)
hold on
% for i = 1:4
%     plot(freq, Llog3(i,:), '-o', 'Linewidth', wide) % Fourier e5
% end
for i = 1:4
    plot(freq, Llog2(i,:), '-o', 'Linewidth', wide) % Fourier e7
end
for i = 1:4
    plot(freq, Llog(i,:), '--*', 'Linewidth', thin) % Fourier e9 -- prev: ['-.', 'Linewidth', wide]
end
for i = 1:4
    plot(freq, Plog(i,:), ':x', 'Linewidth', wide) % Gardner e7 single axis
end
% for i = 1:4
%     plot(freq, Mlog(i,:), '-.', 'Linewidth', wide) % Gardner e7 dual axis -- prev: ['-.diamond', 'Linewidth', thin]
% end
colororder(newcolors)
% ylim([-2 0])
% xlim([0 0.5])
ylabel("Log10 Rd")
xlabel("Spatial freq. (mm^{-1})")
yticks([-2 -1.5 -1 -0.5 0])
yticklabels(["10^{-2}", "10^{-1.5}", "10^{-1}", "10^{-0.5}", "10^{0}"])
legend("F. e7 L*=0.5","F. e7 L*=1.0","F. e7 L*=2.0","F. e7 L*=4.0",...
       "F. e9 L*=0.5","F. e9 L*=1.0","F. e9 L*=2.0","F. e9 L*=4.0",...
       "P. e7 L*=0.5","P. e7 L*=1.0","P. e7 L*=2.0","P. e7 L*=4.0")%,...
%        "P. e7 2a L*=0.5","P. e7 2a L*=1.0","P. e7 2a L*=2.0","P. e7 2a L*=4.0")
grid on
[t,s] = title(['Diffuse Reflectance vs. Spatial freq. (mm^{-1}) w/ varying L* values, photon counts, and simulation methods'],...
              ['1mm voxels, g=0.9 Fourier e7, & e9 vs. Pencil e7 singleSources']);

%% Voxel size differences
          
figure(6) %e7 fourier source w/ 1mm voxels vs. 4e7 fourier source w/ 0.5mm voxels
hold on
for i = 1:4
    plot(freq, Plog3(i,:), '-o', 'Linewidth', wide)
end
for i = 1:4
    plot(freq, Nlog2(i,:), '--*', 'Linewidth', thin)
end
for i = 1:4
    plot(freq, Nlog(i,:), '-.', 'Linewidth', wide)
end
colororder(newcolors)
% ylim([-2 0])
% xlim([0 0.5])
ylabel("Log10 Rd")
xlabel("Spatial freq. (mm^{-1})")
yticks([-2 -1.5 -1 -0.5 0])
yticklabels(["10^{-2}", "10^{-1.5}", "10^{-1}", "10^{-0.5}", "10^{0}"])
legend("1.0mm vox. 1e7 photons",'','','',...
       "0.5mm vox. 1e7 photons",'','','',...
       "0.5mm vox. 4e7 photons",'','','')
% legend("1.0mm vox. 1e7 L*=0.5","1.0mm vox. 1e7 L*=1.0","1.0mm vox. 1e7 L*=2.0","1.0mm vox. 1e7 L*=4.0",...
%        "0.5mm vox. 1e7 L*=0.5","0.5mm vox. 1e7 L*=1.0","0.5mm vox. 1e7 L*=2.0","0.5mm vox. 1e7 L*=4.0",...
%        "0.5mm vox. 4e7 L*=0.5","0.5mm vox. 4e7 L*=1.0","0.5mm vox. 4e7 L*=2.0","0.5mm vox. 4e7 L*=4.0")
grid on
[t,s] = title(['MCX Fourier R_d vs. Spatial freq. (mm^{-1}) at diff. voxel sizes & photon counts'],...
              ['Fourier Source g=0.71 1e7 w/1.0mm vox. vs. 1e7 w/0.5mm vox. vs 4e7 w/0.5mm vox.']);
          

%% Dual vs Single Axis

figure(7) % e5 pencil sources vs. e7 pencil sources
hold on
for i = 1:4
    plot(freq, Mlog2(i,:), '-o', 'Linewidth', wide) % Dual Axis e5
end
for i = 1:4
    plot(freq, Mlog(i,:), '--*', 'Linewidth', thin) % Dual Axis e7
end
for i = 1:4
    plot(freq, Plog2(i,:), '-.diamond', 'Linewidth', wide) % Single Axis e5
end
for i = 1:4
    plot(freq, Plog(i,:), ':x', 'Linewidth', wide) % Single Axis e7
end
colororder(newcolors)
ylabel("Log10 Rd")
xlabel("Spatial freq. (mm^{-1})")
yticks([-2 -1.5 -1 -0.5 0])
yticklabels(["10^{-2}", "10^{-1.5}", "10^{-1}", "10^{-0.5}", "10^{0}"])
legend("P. e5 1a L*=0.5","P. e5 1a L*=1.0","P. e5 1a L*=2.0","P. e5 1a L*=4.0",...
       "P. e7 1a L*=0.5","P. e7 1a L*=1.0","P. e7 1a L*=2.0","P. e7 1a L*=4.0",...
       "P. e5 2a L*=0.5","P. e5 2a L*=1.0","P. e5 2a L*=2.0","P. e5 2a L*=4.0",...
       "P. e7 2a L*=0.5","P. e7 2a L*=1.0","P. e7 2a L*=2.0","P. e7 2a L*=4.0")
grid on
[t,s] = title(['Diffuse Reflectance vs. Spatial freq. (mm^{-1}) w/ varying L* values and photon counts'],...
              ['1mm voxels, g=0.9 Single vs Dual Axis Pencil e5 vs. Pencil e7 Sources']);
          
%% Dummy Plot to generate Legend
figure(8)
hold on
for i = 1:4
    plot(freq,Mlog(i,:),'-', 'Linewidth', wide)
end
colororder(newcolors)
legend({"L*=0.5","L*=1.0","L*=2.0","L*=4.0"},'Fontsize',18)

