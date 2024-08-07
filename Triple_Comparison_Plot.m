%% Values to Graph
clear all

freq = [0, 0.05/3, 0.10/3, 0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5];

curfold = pwd;

e = dir(curfold +"\\SFDI_dref_nphoton_e7_g071" + "\\SFDI_drefval_g_0.71_*.csv");
names = {e.name};
for i=1:length(e)
    P(i,:) = readmatrix(string(names(i)));
end

ee = dir(curfold +"\\SFDI_dref_nphoton_e8_multilayer" + "\\SFDI_multilayer_drefval_2layermodel*.csv");
names = {ee.name};
for i=1:length(ee)
    P1(i,:) = readmatrix(string(names(i)));
end

h = dir(curfold + "\\SFDI_drefcalc" + "\\SFDI_drefcalc_g_adjust_dual_axis*.csv");
names = {h.name};
for i=1:length(h)
    Q((i*4-3):(i*4),:) = readmatrix(string(names(i)));
end

hh = dir(curfold + "\\SFDI_drefcalc" + "\\SFDI_drefcalc_g_0.71_single*1.0e+07*.csv"); % single axis Gardner MCX
names = {hh.name};
for i=1:length(hh)
    Q1((i*4-3):(i*4),:) = readmatrix(string(names(i)));
end

hhh = dir(curfold+ "\\SFDI_drefcalc\\SFDI_drefcalc_g_0.71_dual*_nphoton_1.0e+07*.csv"); % dual axis Gardner MCX
names = {hhh.name};
for i=1:length(hhh)
    Q2((i*4-3):(i*4),:) = readmatrix(string(names(i)));
end

d = dir(curfold + "\\SFDI_MCCL_Gardner_Results.csv");
names = {d.name};
for i=1:length(d)
    R((i*4-3):(i*4),:) = readmatrix(string(names(i)));
end

dd = dir(curfold + "\\SFDI_dref_nphoton_e7_test\\SFDI_approx*.csv");
names = {dd.name};
for i=1:length(dd)
    S(i,:) = readmatrix(string(names(i)));
end

Plog = log10(P);
P1log = log10(P1);
Qlog = log10(Q);
Q1log = log10(Q1);
Q2log = log10(Q2);
Rlog = log10(abs(R));
Slog = log10(S);

% [0.9290    0.6040    0.1250
%      0    0.4470    0.7410
% 1.0000    0.8600    0.1250
% 0.4940    0.1840    0.5560
% 1.0000    0.5000    0.8000
% 0.2330    0.3370    0.0940
% 0.3010    0.7450    0.9330
% 0.6350    0.0780    0.1840]
%% Plot
close all

wide = 2;
thin = 1;

newcolors =    [0.9290    0.6040    0.1250
                     0    0.4470    0.7410
                0.1250    0.8600    0.1250
                0.4940    0.1840    0.5560];
            
%% Three-Way

figure(1) % e7 fourier source w/ g=0.71 vs e7 pencil source w/ g=0.71
hold on
for i = 1:4
    plot(freq, Plog(i,:), '--*', 'Linewidth', thin)
end
for i = 1:4
    plot(freq, Q1log(i,:), '-o', 'Linewidth', wide)
end
for i = 1:4
    plot(0:0.01:0.50, Rlog(i,:), '-.', 'Linewidth', wide)
end
hold off
colororder(newcolors)
ylim([-2 0])
xlim([0 0.5])

ylabel("Log10 R_d")
xlabel("Spatial freq. (mm^{-1})")
yticks([-2 -1.5 -1 -0.5 0])
yticklabels(["10^{-2}", "10^{-1.5}", "10^{-1}", "10^{-0.5}", "10^{0}"])
legend( "MCX Fourier L*=0.5","MCX Fourier L*=1.0","MCX Fourier L*=2.0","MCX Fourier L*=4.0",...
        "MCX Gardner L*=0.5","MCX Gardner L*=1.0","MCX Gardner L*=2.0","MCX Gardner L*=4.0",...
        "MCCL Gardner L*=0.5","MCCL Gardner L*=1.0","MCCL Gardner L*=2.0","MCCL Gardner L*=4.0")
grid on
[t,s] = title(['Diffuse Reflectance vs. Spatial freq. (mm^{-1}) at diff. L* values'],...
              ['1mm voxels, 1e7 photons g=0.71 Fourier source vs Pencil source']);

%%% Three-Way non-Log
figure(2) % e7 fourier source w/ g=0.71 vs e7 pencil source w/ g=0.71
hold on
for i = 1:4
    plot(freq, P(i,:), '--*', 'Linewidth', thin)
end
for i = 1:4
    plot(freq, Q1(i,:), '-o', 'Linewidth', wide)
end
for i = 1:4
    plot(0:0.01:0.50, R(i,:), '-.', 'Linewidth', wide)
end
hold off
colororder(newcolors)
ylim([0 0.75])
xlim([0 0.5])

ylabel("R_d")
xlabel("Spatial freq. (mm^{-1})")
%yticks([-2 -1.5 -1 -0.5 0])
%yticklabels(["10^{-2}", "10^{-1.5}", "10^{-1}", "10^{-0.5}", "10^{0}"])
legend( "MCX Fourier L*=0.5","MCX Fourier L*=1.0","MCX Fourier L*=2.0","MCX Fourier L*=4.0",...
        "MCX Gardner L*=0.5","MCX Gardner L*=1.0","MCX Gardner L*=2.0","MCX Gardner L*=4.0",...
        "MCCL Gardner L*=0.5","MCCL Gardner L*=1.0","MCCL Gardner L*=2.0","MCCL Gardner L*=4.0")
grid on
[t,s] = title(['Diffuse Reflectance vs. Spatial freq. (mm^{-1}) at diff. L* values'],...
              ['1mm voxels, 1e7 photons g=0.71 Fourier source vs Pencil source']);
          
%% Two/Three-Way Close-Up
figure(3) % e7 fourier source w/ g=0.71 vs e7 pencil source w/ g=0.71
hold on
for i = 1:4
    plot(freq, Plog(i,:), '--*', 'Linewidth', thin)
end
for i = 1:4
    plot(freq, Qlog(i,:), '-o', 'Linewidth', wide)
end
for i = 1:4
    plot(0:0.01:0.50, Rlog(i,:), '-.', 'Linewidth', wide)
end
hold off
colororder(newcolors)
ylim([-1.2 0])
xlim([0 0.1])

ylabel("Log10 R_d")
xlabel("Spatial freq. (mm^{-1})")
% yticks([-2 -1.5 -1 -0.5 0]);
% yticklabels(["10^{-2}", "10^{-1.5}", "10^{-1}", "10^{-0.5}", "10^{0}"])
yticks([-1 -0.5 0]);
yticklabels(["10^{-1}", "10^{-0.5}", "10^{0}"])
legend( "MCX Fourier L*=0.5","MCX Fourier L*=1.0","MCX Fourier L*=2.0","MCX Fourier L*=4.0",...
        "MCX Gardner L*=0.5","MCX Gardner L*=1.0","MCX Gardner L*=2.0","MCX Gardner L*=4.0",...
        "MCCL Gardner L*=0.5","MCCL Gardner L*=1.0","MCCL Gardner L*=2.0","MCCL Gardner L*=4.0")
grid on
[t,s] = title(['Multilayer Diffuse Reflectance vs. Spatial freq. (mm^{-1}) at diff. L* values'],...
              ['1mm voxels, g=0.71 e7 Fourier source vs e7 Pencil source vs MCCL Gardner']);

%% Single vs Dual Axis?
figure(5) % comparing dual axis vs. single axis gardner method vs fourier
hold on
for i = 1:4
    plot(freq, Q1log(i,:), '--*', 'Linewidth', thin)
end
for i = 1:4
    plot(freq, Q2log(i,:), '-o', 'Linewidth', wide)
end
for i = 1:4
    plot(freq, Plog(i,:), '-.', 'Linewidth', wide)
end
hold off
colororder(newcolors);
ylim([-2.1 0]);
xlim([0 0.5]);

ylabel("Log10 R_d");
xlabel("Spatial freq. (mm^{-1})");
yticks([-2 -1.5 -1 -0.5 0]);
yticklabels(["10^{-2}", "10^{-1.5}", "10^{-1}", "10^{-0.5}", "10^{0}"]);
legend( "MCX Dual Axis L*=0.5","MCX Dual Axis L*=1.0","MCX Dual Axis L*=2.0","MCX Dual Axis L*=4.0",...
        "MCX Single Axis L*=0.5","MCX Single Axis L*=1.0","MCX Single Axis L*=2.0","MCX Single Axis L*=4.0",...
        "MCX Fourier L*=0.5","MCX Fourier L*=1.0","MCX Fourier L*=2.0","MCX Fourier L*=4.0")
grid on
[t,s] = title(['Multilayer Diffuse Reflectance vs. Spatial freq. (mm^{-1}) at diff. L* values'],...
              ['1mm voxels, e7 g=0.71 Dual Axis Gardner vs Single Axis Gardner vs Fourier']);
          
%% Four-Way
figure(7) % e7 fourier source w/ g=0.71 vs e7 pencil source w/ g=0.71 vs MCCL vs DiffApprox
hold on
for i = 1:4
    plot(freq, Plog(i,:), '--*', 'Linewidth', thin)
end
for i = 1:4
    plot(freq, Q1log(i,:), '-o', 'Linewidth', wide)
end
for i = 1:4
    plot(0:0.01:0.50, Rlog(i,:), '-.', 'Linewidth', wide)
end
for i = 1:4
    plot(freq, Slog(i,:), ':x','Linewidth', wide)
end
hold off
colororder(newcolors)
ylim([-2.1 0])
xlim([0 0.5])

ylabel("Log10 R_d")
xlabel("Spatial freq. (mm^{-1})")
yticks([-2 -1.5 -1 -0.5 0])
yticks([-1 -0.5 0])
yticklabels(["10^{-2}", "10^{-1.5}", "10^{-1}", "10^{-0.5}", "10^{0}"])
yticklabels(["10^{-1}", "10^{-0.5}", "10^{0}"])
legend( "MCX Fourier L*=0.5","MCX Fourier L*=1.0","MCX Fourier L*=2.0","MCX Fourier L*=4.0",...
        "MCX Gardner L*=0.5","MCX Gardner L*=1.0","MCX Gardner L*=2.0","MCX Gardner L*=4.0",...
        "MCCL Gardner L*=0.5","MCCL Gardner L*=1.0","MCCL Gardner L*=2.0","MCCL Gardner L*=4.0",...
        "DiffApprox L*=0.5","DiffApprox L*=1.0","DiffApprox L*=2.0","DiffApprox L*=4.0")
grid on
[t,s] = title(['Diffuse Reflectance vs. Spatial freq. (mm^{-1}) w/ Varying L* Values and Simulation Methods'],...
              ['1mm voxels, 1e7 photons g=0.71 Fourier vs Gardner vs Differential Approximation']);
          
          
%% Three-Way w/ Diff Approx
f = figure(9); % e7 fourier source w/ g=0.0.71 vs e7 pencil source w/ g=0.71
%f.Position = [10,10,1200,800];
hold on
for i = 1:4
    plot(freq, Plog(i,:), '--*', 'Linewidth', thin)
end
for i = 1:4
    plot(freq, Q1log(i,:), '-o', 'Linewidth', wide)
end
for i = 1:4
    plot(freq, Slog(i,:), '-.','Linewidth', wide)
end
hold off
colororder(newcolors)
ylim([-2.1 0])
xlim([0 0.5])

ylabel("Log10 R_d")
xlabel("Spatial freq. (mm^{-1})")
yticks([-2 -1.5 -1 -0.5 0])
yticklabels(["10^{-2}", "10^{-1.5}", "10^{-1}", "10^{-0.5}", "10^{0}"])
legend( "MCX Fourier L*=0.5","MCX Fourier L*=1.0","MCX Fourier L*=2.0","MCX Fourier L*=4.0",...
        "MCX Gardner L*=0.5","MCX Gardner L*=1.0","MCX Gardner L*=2.0","MCX Gardner L*=4.0",...
        "DiffApprox L*=0.5","DiffApprox L*=1.0","DiffApprox L*=2.0","DiffApprox L*=4.0")
grid on
[t,s] = title(['Diffuse Reflectance vs. Spatial freq. (mm^{-1}) at diff. L* values'],...
              ['Fourier src vs Pencil src vs Diffusion Approximation']);

%%% Non-Log
f = figure(10); % e7 fourier source w/ g=0.0.71 vs e7 pencil source w/ g=0.71
%f.Position = [10,10,1200,800];
hold on
for i = 1:4
    plot(freq, P(i,:), '--*', 'Linewidth', thin)
end
for i = 1:4
    plot(freq, Q1(i,:), '-o', 'Linewidth', wide)
end
for i = 1:4
    plot(freq, S(i,:), '-.','Linewidth', wide)
end
hold off
colororder(newcolors)
ylim([0 0.75])
xlim([0 0.5])

ylabel("R_d")
xlabel("Spatial freq. (mm^{-1})")
%yticks([-2 -1.5 -1 -0.5 0])
%yticklabels(["10^{-2}", "10^{-1.5}", "10^{-1}", "10^{-0.5}", "10^{0}"])
legend( "MCX Fourier L*=0.5","MCX Fourier L*=1.0","MCX Fourier L*=2.0","MCX Fourier L*=4.0",...
        "MCX Gardner L*=0.5","MCX Gardner L*=1.0","MCX Gardner L*=2.0","MCX Gardner L*=4.0",...
        "DiffApprox L*=0.5","DiffApprox L*=1.0","DiffApprox L*=2.0","DiffApprox L*=4.0")
grid on
[t,s] = title(['Diffuse Reflectance vs. Spatial freq. (mm^{-1}) at diff. L* values'],...
              ['Fourier src vs Pencil src vs Diffusion Approximation']);
          
%% Three-Way w/Dual Axis?
f = figure(11); % e7 fourier source w/ g=0.71 vs e7 pencil source w/ g=0.71
%f.Position = [10,10,1200,800];
hold on
for i = 1:4
    plot(freq, (P(i,:)./Q1(i,:)), '--*', 'Linewidth', thin)
end
% for i = 1:4
%     plot(freq, Q2log(i,:), '-o', 'Linewidth', wide)
% end
hold off
colororder(newcolors)
%ylim([1 1.2])
xlim([0 0.5])
hold off
ylabel("R_d")
xlabel("Spatial freq. (mm^{-1})")
%yticks([-2 -1.5 -1 -0.5 0])
%yticklabels(["10^{-2}", "10^{-1.5}", "10^{-1}", "10^{-0.5}", "10^{0}"])
legend( "F/G L*=0.5","F/G L*=1.0","F/G L*=2.0","F/G L*=4.0")
grid on
[t,s] = title(['Diffuse Reflectance vs. Spatial freq. (mm^{-1}) at diff. L* values'],...
              ['MCX Fourier src relative to MCX Gardner']);
          
%% MCCL vs MCX Gardner Presentation
f = figure(12); % e7 mcx gardner vs mccl gardner

hold on

for i = 1:4
    plot(freq, Q1log(i,:), '-o', 'Linewidth', wide)
end
for i = 1:4
    plot(0:0.01:0.50, Rlog(i,:), '-.', 'Linewidth', wide)
end
hold off
colororder(newcolors)
ylim([-2 0])
xlim([0 0.5])

ylabel("Log10 R_d")
xlabel("Spatial freq. (mm^{-1})")
yticks([-2 -1.5 -1 -0.5 0])
yticklabels(["10^{-2}", "10^{-1.5}", "10^{-1}", "10^{-0.5}", "10^{0}"])
legend( "MCX Gardner L*=0.5","MCX Gardner L*=1.0","MCX Gardner L*=2.0","MCX Gardner L*=4.0",...
        "MCCL Gardner L*=0.5","MCCL Gardner L*=1.0","MCCL Gardner L*=2.0","MCCL Gardner L*=4.0")
grid on
[t,s] = title(['MCX vs. MCCL: Diffuse Reflectance vs. Spatial freq. (mm^{-1}) at diff. L* values'],...
              ['1mm voxels, 1e7 photons g=0.71']);

%% MCX Fourier v Gardner Presentation
f = figure(13); % e7 mcx gardner vs mccl gardner

hold on

for i = 1:4
    plot(freq, Plog(i,:), '--*', 'Linewidth', thin)
end
for i = 1:4
    plot(freq, Q1log(i,:), '-o', 'Linewidth', wide)
end
hold off
colororder(newcolors)
ylim([-2 0])
xlim([0 0.5])

ylabel("Log10 R_d")
xlabel("Spatial freq. (mm^{-1})")
yticks([-2 -1.5 -1 -0.5 0])
yticklabels(["10^{-2}", "10^{-1.5}", "10^{-1}", "10^{-0.5}", "10^{0}"])
legend( "MCX Fourier L*=0.5","MCX Fourier L*=1.0","MCX Fourier L*=2.0","MCX Fourier L*=4.0",...
        "MCX Gardner L*=0.5","MCX Gardner L*=1.0","MCX Gardner L*=2.0","MCX Gardner L*=4.0")
grid on
[t,s] = title(['Fourier vs Gardner: Diffuse Reflectance vs. Spatial freq. (mm^{-1}) at diff. L* values'],...
              ['1mm voxels, 1e7 photons g=0.71']);
