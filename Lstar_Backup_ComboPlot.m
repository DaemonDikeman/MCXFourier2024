%% Values to Graph
clear all

freq = [0, 0.05/3, 0.10/3, 0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5];
freq2 = [0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.10];

curfold = pwd;

d = dir(curfold +"\\SFDI_dref_nphoton_e9" + "\\SFDI_drefval_*.csv");
names = {d.name};
for i=1:length(d)
    L(i,:) = readmatrix(string(names(i)));
end

dd = dir(curfold +"\\SFDI_dref_nphoton_e7" + "\\SFDI_drefval_*.csv");
names = {dd.name};
for i=1:length(dd)
    L2(i,:) = readmatrix(string(names(i)));
end

ddd = dir(curfold +"\\SFDI_dref_nphoton_e5_g09" +"\\SFDI_drefval_*.csv");
names = {ddd.name};
for i=1:length(ddd)
    L3(i,:) = readmatrix(string(names(i)));
end
e = dir(curfold +"\\SFDI_drefcalc" + "\\SFDI_drefcalc_nphoton_1.0e+07*.csv");
names = {e.name};
for i=1:length(e)
    M((i*4-3):(i*4),:) = readmatrix(string(names(i)));
end

ee = dir(curfold +"\\SFDI_drefcalc" + "\\SFDI_drefcalc_nphoton_1.0e+05*.csv");
names = {ee.name};
for i=1:length(ee)
    M2((i*4-3):(i*4),:) = readmatrix(string(names(i)));
end

f = dir(curfold +"\\SFDI_dref_nphoton_4e7_half_mm_voxel" + "\\SFDI_drefval_*.csv");
names = {f.name};
for i=1:length(f)
    N(i,:) = readmatrix(string(names(i)));
end

ff = dir(curfold +"\\SFDI_dref_nphoton_e7_half_mm_voxel" + "\\SFDI_drefval_*.csv");
names = {ff.name};
for i=1:length(ff)
    N2(i,:) = readmatrix(string(names(i)));
end

h = dir(curfold +"\\SFDI_dref_nphoton_e7_g077" + "\\SFDI_drefval_*.csv");
names = {h.name};
for i=1:length(h)
    P(i,:) = readmatrix(string(names(i)));
end

hh = dir(curfold + "\\SFDI_drefcalc_g_adjust_*.csv");
names = {hh.name};
for i=1:length(hh)
    Q((i*4-3):(i*4),:) = readmatrix(string(names(i)));
end

k = dir(curfold +"\\SFDI_focused_freq" + "\\SFDI_drefval_g_0.71_*.csv");
names = {k.name};
for i=1:length(k)
    P2(i,:) = readmatrix(string(names(i)));
end

kk = dir(curfold +"\\SFDI_focused_freq" + "\\SFDI_drefcalc_g_adjust_*.csv");
names = {kk.name};
for i=1:length(kk)
    Q2((i*4-3):(i*4),:) = readmatrix(string(names(i)));
end

hhh = dir(curfold +"\\SFDI_dref_nphoton_e7_g071" + "\\SFDI_drefval_*.csv");
names = {hhh.name};
for i=1:length(hhh)
    P3(i,:) = readmatrix(string(names(i)));
end

Mlog = log10(M);
Mlog2 = log10(M2);
Llog = log10(L);
Llog2 = log10(L2);
Llog3 = log10(L3);
Nlog = log10(N);
Nlog2 = log10(N2);
Plog = log10(P);
Qlog = log10(Q);
Plog2 = log10(P2);
Qlog2 = log10(Q2);
Plog3 = log10(P3);

%% Graph
close all

wide = 2;
thin = 1;

newcolors =    [0.9290    0.6040    0.1250
                     0    0.4470    0.7410
                0.1250    0.8600    0.1250
                0.4940    0.1840    0.5560];

% ['Diffuse Reflectance vs. Spatial freq. (mm^{-1}) at diff. L* values']

%% g 0.9 vs g 0.77 
figure(8) % e7 fourier source w/ g=0.9 vs e7 fourier source w/ g=0.77
hold on
for i = 1:length(d)
    plot(freq, Llog2(i,:), '-o', 'Linewidth', wide)
    plot(freq, Plog(i,:), '--*', 'Linewidth', thin)
end
colororder(newcolors)
ylim([-2.1 0])
xlim([0 0.5])

ylabel("Log10 Rd")
xlabel("Spatial freq. (mm^{-1})")
yticks([-2 -1.5 -1 -0.5 0])
yticklabels(["10^{-2}", "10^{-1.5}", "10^{-1}", "10^{-0.5}", "10^{0}"])
legend("g=0.9 L*=0.5","g=0.77 L*=0.5",...
       "g=0.9 L*=1.0","g=0.77 L*=1.0",...
       "g=0.9 L*=2.0","g=0.77 L*=2.0",...
       "g=0.9 L*=4.0","g=0.77 L*=4.0")
grid on
[t,s] = title(['Diffuse Reflectance vs. Spatial freq. (mm^{-1}) at diff. L* values'],...
              ['1mm voxels, Fourier 1e7 source vs g=0.9 vs. g=0.77']);
          
figure(9) % e7 fourier source w/ g=0.9 vs e7 fourier source w/ g=0.77
hold on
for i = 1:length(d)
    plot(freq, L2(i,:), '-o', 'Linewidth', wide)
    plot(freq, P(i,:), '--*', 'Linewidth', thin)
end
colororder(newcolors)
% ylim([-2.1 0])
% xlim([0 0.5])

ylabel("Rd")
xlabel("Spatial freq. (mm^{-1})")
% yticks([-2 -1.5 -1 -0.5 0])
% yticklabels(["10^{-2}", "10^{-1.5}", "10^{-1}", "10^{-0.5}", "10^{0}"])
legend("g=0.9 L*=0.5","g=0.77 L*=0.5",...
       "g=0.9 L*=1.0","g=0.77 L*=1.0",...
       "g=0.9 L*=2.0","g=0.77 L*=2.0",...
       "g=0.9 L*=4.0","g=0.77 L*=4.0")
grid on
[t,s] = title(['Diffuse Reflectance vs. Spatial freq. (mm^{-1}) at diff. L* values'],...
              ['1mm voxels, Fourier 1e7 source vs g=0.9 vs. g=0.77']);
          
          
figure(10) % e7 pencil source w/ g=0.9 vs e7 pencil source w/ g=0.77
hold on
for i = 1:length(d)
    plot(freq, Mlog(i,:), '-o', 'Linewidth', wide)
    plot(freq, Qlog(i,:), '--*', 'Linewidth', thin)
end
colororder(newcolors)
ylim([-2.1 0])
xlim([0 0.5])

ylabel("Log10 Rd")
xlabel("Spatial freq. (mm^{-1})")
yticks([-2 -1.5 -1 -0.5 0])
yticklabels(["10^{-2}", "10^{-1.5}", "10^{-1}", "10^{-0.5}", "10^{0}"])
legend("g=0.9 L*=0.5","g=0.77 L*=0.5",...
       "g=0.9 L*=1.0","g=0.77 L*=1.0",...
       "g=0.9 L*=2.0","g=0.77 L*=2.0",...
       "g=0.9 L*=4.0","g=0.77 L*=4.0")
grid on
[t,s] = title(['Diffuse Reflectance vs. Spatial freq. (mm^{-1}) at diff. L* values'],...
              ['1mm voxels, Pencil 1e7 source g=0.9 vs. g=0.77']);

%% Comparing Pencil to fourier w/ g=0.71
figure(11) % e7 fourier source w/ g=0.71 vs e7 pencil source w/ g=0.71
hold on
for i = 1:length(d)
    plot(freq, Plog(i,:), '-o', 'Linewidth', wide)
    plot(freq, Qlog(i,:), '--*', 'Linewidth', thin)
end
colororder(newcolors)
ylim([-2.1 0])
xlim([0 0.5])

ylabel("Log10 Rd")
xlabel("Spatial freq. (mm^{-1})")
yticks([-2 -1.5 -1 -0.5 0])
yticklabels(["10^{-2}", "10^{-1.5}", "10^{-1}", "10^{-0.5}", "10^{0}"])
legend("Fourier L*=0.5","Pencil L*=0.5",...
       "Fourier L*=1.0","Pencil L*=1.0",...
       "Fourier L*=2.0","Pencil L*=2.0",...
       "Fourier L*=4.0","Pencil L*=4.0")
grid on
[t,s] = title(['Diffuse Reflectance vs. Spatial freq. (mm^{-1}) at diff. L* values'],...
              ['1mm voxels, 1e7 photons g=0.71 Fourier source vs Pencil source']);
          
figure(12) % e7 fourier source w/ g=0.77 vs e7 pencil source w/ g=0.77
hold on
for i = 1:length(d)
    plot(freq2, Plog2(i,:), '-o', 'Linewidth', wide)
    plot(freq2, Qlog2(i,:), '--*', 'Linewidth', thin)
end
colororder(newcolors)
ylim([-1.3 0])
xlim([0 0.1])

ylabel("Log10 Rd")
xlabel("Spatial freq. (mm^{-1})")
yticks([-1.25, -1 -0.75 -0.5 -0.25 0])
yticklabels(["10^{-1.25}", "10^{-1}", "10^{-0.75}", "10^{-0.5}", "10^{-0.25}", "10^{0}"])
legend("Fourier L*=0.5","Pencil L*=0.5",...
       "Fourier L*=1.0","Pencil L*=1.0",...
       "Fourier L*=2.0","Pencil L*=2.0",...
       "Fourier L*=4.0","Pencil L*=4.0")
grid on
[t,s] = title(['Diffuse Reflectance vs. Spatial freq. (mm^{-1}) at diff. L* values'],...
              ['1mm voxels, narrow freq. range, 1e7 photons g=0.71 Fourier source vs Pencil source']);