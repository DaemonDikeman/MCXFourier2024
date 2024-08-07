%% Calibrate Data Set
% This function calculates diffuse reflectance (Rd) for a demodulated data
% set when supplied the phantom ID and the file containing the known
% phantom optical properties. Optical properties are extracted with matlabs
% linear interpolation function for the aquisition wavelengths in the data
% set. If the OP file is sparsly sampled or not in the wavelength range,
% interp1 will still try to determine the optical properties for the
% phantom. Make sure the OP file represents the aquisition wavelengths
% well.
% By default, this function uses a homogeneous LUT to determine expected
% diffuse reflectance of the phantom. The LUT contains
% Fx = [0,0.05,0.1,0.15,0.2]. Load a different LUT if the data set contains
% a different Fx.
function dataSet = calibDataSet(dataSet,phantomID,phantomOPFile,useCorrected)
    sampleIdx = 1:length(dataSet);
    sampleIdx(phantomID) = [];
    % get size information from the data set
    [Ny,Nx,Nfreqs,~,Nwavelengths] = size(dataSet{phantomID}.Raw);
    
    % get phantom optical properties
    disp('Getting Phantom Optical Properties ...');
    phantomFile = load(phantomOPFile);
    phantomOPs = phantomFile.phantomOPs;
    phantomOPs = [interp1(phantomOPs(:,1),phantomOPs(:,2),...
        dataSet{phantomID}.Wavelengths,...
        'pchip','extrap');interp1(phantomOPs(:,1),...
        phantomOPs(:,3),dataSet{phantomID}.Wavelengths,'pchip','extrap')]';
    Mua = phantomOPs(:,1);
    Musp = phantomOPs(:,2);
    
    % load the LUT... This can be changed out for any LUT with the same
    % format
    disp('Loading Lookup Table ...');
    LUT_H = load('LUT_Homogeneous');
    
    % Use 2-D interpolation to run the Forward model for calibration
    disp('Running Forward Model ...');
    ForwardRd = nan(Nfreqs,Nwavelengths);
    for i = 1:Nfreqs
        fx = dataSet{phantomID}.Freqs(i);
        idx = LUT_H.LUT.Fx == fx;
        if all(~idx)
            error(['LUT doesn''t contain spatial frequency' num2str(fx)]);
        end
        ForwardRd(i,:) = griddata(LUT_H.LUT.Mua,LUT_H.LUT.Musp,...
            LUT_H.LUT.Rd(:,:,idx),Mua,Musp);
    end
    if useCorrected == false
    % Calibrate the data set
    disp('Calibrating ...');
    for i = sampleIdx
        %Calibrate
        dataSet{i}.Rd = nan(Ny,Nx,Nfreqs,Nwavelengths);
        for j = 1:length(dataSet{phantomID}.Freqs)
            for k = 1:Nwavelengths
                dataSet{i}.Rd(:,:,j,k) = (dataSet{i}.Demodulated(:,:,j,k)...
                    ./dataSet{phantomID}.Demodulated(:,:,j,k))*...
                    ForwardRd(j,k);
            end
        end
    end
    end

    if useCorrected == true
    % Calibrate the data set
    disp('Calibrating ...');
    for i = sampleIdx
        %Calibrate
        dataSet{i}.Rd = nan(Ny,Nx,Nfreqs,Nwavelengths);
        for j = 1:length(dataSet{phantomID}.Freqs)
            for k = 1:Nwavelengths
                dataSet{i}.Rd(:,:,j,k) = (dataSet{i}.AngleCorrected(:,:,j,k)...
                    ./dataSet{phantomID}.Demodulated(:,:,j,k))*...
                    ForwardRd(j,k);
            end
        end
    end
    end
    disp('Done');
end %function