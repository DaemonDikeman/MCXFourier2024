%% Define import parameters to create a data set
fx = [0,0.05,0.1,0.15,0.2];
wavelengthRange = [650,998];%[650,852]
format = 'MI';
[dataSet,phantomID] = importDataSet(fx,wavelengthRange,format);
%% Demodulate
dataSet = demodDataSet(dataSet);
%% Generate a mask for analysis and corrections
threshold = 0.25; %between 0 and 1
fxID = 1; % profile frequency index (1 or 2)
heightFactor = -0.97; % scalar that converts phase shift to height
phaseSmoothStd = 1; % standard deviation for Gaussian smoothing kernel (phase map)
% correctionFactorFileName = 'HCFactors'; % file that contains height correction matrix
angleSmoothStd = 1; % standard deviation for Gaussian smoothing kernel (angle map)
% dataSet = genCorrectionMask(dataSet,phantomID,threshold);
% dataSet = genHeightMap(dataSet,phantomID,fxID,heightFactor,phaseSmoothStd);
% dataSet = doHeightCorrection(dataSet,phantomID,correctionFactorFileName);
% dataSet = genAngleMap(dataSet,phantomID);
% dataSet = doAngleCorrection(dataSet,phantomID,angleSmoothStd);
%% Calibrate
%phantomOPFile = 'CarlosBigBoy_Optical_Properties';
phantomOPFile = 'MI_Phantom_Optical_Properties';
useCorrected = false; %use height/angle corrections to calculate Rd
dataSet = calibDataSet(dataSet,phantomID,phantomOPFile,useCorrected);
%% Use 2-fx LUT method
fx = [0,0.1];
%lutFile = 'my_LUT_Fourier_homogenous_100x100';
lutFile = 'my_LUT_Gardner_MCX_homogenous_100x100';
%lutFile = 'my_LUT_Fourier_SAC_multilayer_light_100x100';
%lutFile = 'LUT_Homogeneous_Large_1e6.mat';
dataSet = lutFit_new(dataSet,phantomID,fx,lutFile);
%dataSet = lutFit(dataSet,phantomID,fx);
%% Use 5-fx neural network to fit optical properties
% dataSet = NNFit_5fx(dataSet,phantomID);
%% Fit parameters from Rd space
% dataSet = NNFit_Param(dataSet,phantomID,'A');
% dataSet = NNFit_Param(dataSet,phantomID,'B');
% dataSet = NNFit_Param(dataSet,phantomID,'HbO2');
% dataSet = NNFit_Param(dataSet,phantomID,'HHb');
%% Crop Data/Image set
lambdas = dataSet{1}.Wavelengths;
%dataSet = CropDataSetTwoPoints(dataSet, phantomID, lambdas);
dataSet = CropDataSet(dataSet, phantomID, 100);
dataSet = cropImage(dataSet,phantomID,lambdas);
%% Find Average Optical Properties
dataSet = ROIaverageMua(dataSet,phantomID,lambdas);
%% Plot
k = 2;
for i = 1:(length(dataSet)-1)
    figure(i)
    t = tiledlayout(2,1);
    ax1 = nexttile;
    imagesc(dataSet{i}.MuaCrop(:,:,3))
    title('Mua')
    ax2 = nexttile;
    imagesc(dataSet{i}.MuspCrop(:,:,3))
    title('Musp')
end
%imagesc(dataSet{2}.HHb);

%% Rd Values
for k = 1:(length(dataSet)-1)
    figure
    imagesc(dataSet{k}.Rd(:,:,2,1))
    title('Rd')
    figure
    imagesc(dataSet{k}.Rd(:,:,2,3))
    title('Rd')
end

%% Compiled Averaged Values
L = 2; % 691nm
d = 1;
avmua.homog = zeros(1,d);
avmua.layer = zeros(1,length(dataSet)-(d+1));
avmusp.homog = zeros(1,d);
avmusp.layer = zeros(1,length(dataSet)-(d+1));
for i = 1:d
    avmua.homog(i) = dataSet{i}.avgMua(L);
    avmusp.homog(i) = dataSet{i}.avgMusp(L);
end
for i = 1:(length(dataSet)-(d+1))
    avmua.layer(i) = dataSet{i+d}.avgMua(L);
    avmusp.layer(i) = dataSet{i+d}.avgMusp(L);
end

avRd0.homog = zeros(1,d);
avRd0.layer = zeros(1,length(dataSet)-(d+1));
avRd1.homog = zeros(1,d);
avRd1.layer = zeros(1,length(dataSet)-(d+1));
for i = 1:d
    avRd0.homog(i) = dataSet{i}.avgRd0(L);
    avRd1.homog(i) = dataSet{i}.avgRd1(L);
end
for i = 1:(length(dataSet)-(d+1))
    avRd0.layer(i) = dataSet{i+d}.avgRd0(L);
    avRd1.layer(i) = dataSet{i+d}.avgRd1(L);
end

