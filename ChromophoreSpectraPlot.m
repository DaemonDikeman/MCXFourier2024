%% Importing data sets

clear all
% Importing data to graph
HemRaw  	=   importdata('Hb_Spectra_Prahl.txt') ; %Oxyhemoglobin - dl
nm_range    =   600:2:1200 ;
MelaninData =   (1.70*(10^12))*(nm_range.^(-3.48)) ; %Melanin - calculate? µa = (1.70*(10^12))*(nm_range.^-3.48)
LipidData   =   importdata('fat_spectra.txt'); %Lipids - dl
WaterData1  =   importdata('kou93b_ice_water_spectra.txt'); %H20 lower range
%WaterData2  =   import ; %H20 upper range - dl

%% Processing Hemoglobin/Deoxyhemoglobin data

% µa = (2.303) e (x g/liter)/(64,500 g Hb/mole)
% in whole blood, x ~ 150

OxHemData   = HemRaw.data(:,2).*((2.303*150)/64500);
DeOxHemData = HemRaw.data(:,3).*((2.303*150)/64500);

%% Plotting

figure()
semilogy(HemRaw.data(176:end,1),OxHemData(176:end),'-');
hold on
semilogy(HemRaw.data(176:end,1),DeOxHemData(176:end),'-');
semilogy(nm_range,MelaninData,'-');
semilogy(LipidData.data(172:end,1),LipidData.data(172:end,2),'-');
semilogy(WaterData1.data(226:end,1),WaterData1.data(226:end,2),'-');
%plot(WaterData2(:,1),WaterData2(:,2),'-');
title('Chromophore Absorption Spectrums')
xlabel('Wavelength (nm)')
ylabel('Absorption coefficient (cm^{-1})'   )
legend('O_2Hb','HHb','Melanin','Lipids','H_2O')
hold off