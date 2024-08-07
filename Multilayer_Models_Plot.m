%% Values to Graph
clear all

freq = [0, 0.05/3, 0.10/3, 0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5];

curfold = pwd;

d = dir(curfold +"\\SFDI_multilayer" + "\\SFDI_multilayer_*.csv");
names = {d.name};
for i=1:length(d)
    L(i,:) = readmatrix(string(names(i)));
end

Llog = log10(L);


%% Graph
close all

w1 = 2;
w2 = 1;

newcolors =    [0.9290    0.6040    0.1250
                     0    0.4470    0.7410
                1.0000    0.8600    0.1250
                0.4940    0.1840    0.5560
                1.0000    0.5000    0.8000
                0.2330    0.3370    0.0940
                0.3010    0.7450    0.9330
                0.6350    0.0780    0.1840];

% ['Diffuse Reflectance vs. Spatial freq. (mm^{-1}) at diff. L* values']

figure(1) % e9 fourier source vs. e5 pencil source
hold on
for i = 1:length(d)
    if i < 4
        linestyle = '-o';
    elseif i == 10
        linestyle = '--*';
    else
        linestyle = '-.|';
    end
    plot(freq, Llog(i,:), linestyle, 'Linewidth', w2) % increase line width
end
colororder(newcolors)
ylim([-2.1 0])
xlim([0 0.5])
ylabel("Log10 Rd")
xlabel("Spatial freq. (mm^{-1})")
yticks([-2 -1.5 -1 -0.5 0])
yticklabels(["10^{-2}", "10^{-1.5}", "10^{-1}", "10^{-0.5}", "10^{0}"])
legend("1 Layer (A)","1 Layer (B)","1 Layer (C)",...
       "2 Layer (A+B)","2 Layer (AB+)","2 Layer (A+C)","2 Layer (AC+)",...
       "2 Layer (B+C)","2 Layer (BC+)","3 Layer (ABC)")
grid on
[t,s] = title(['Diffuse Reflectance vs. Spatial Freq. (mm^{-1})'],...
              ['1mm voxels for different multilayer models']);