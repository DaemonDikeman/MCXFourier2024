%% import sfdi data and format the data set structure
function [dataSet,phantomID] = importDataSetMI(Freqs,WavelengthsBounds,path)
%Check format string
    fileIdx = 1; % 1 is MI system
    try % to open a directory to get scan files
        directory = path;
        files = dir(directory);
    catch % error and return null values
        dataSet = NaN;
        phantomID = NaN;
        return
    end
    % reformat filenames to an array
    names = {files.name};
    filenames = names(3:end); %ignore the first two entries
    dataSet = cell(1,length(filenames)); % define dataSet
    for i = 1:length(filenames)
        disp(['Loading ',num2str(i), ' out of ', num2str(length(dataSet)),' ... ']);
        filepath = [directory '\' filenames{i}];
        % load data from raw files
        [rawData,rawProfile,exposureTimes,wavelengths,freqs] = getMIRaw(filepath,fileIdx);
        % import only specified spatial frequencies
        idxf = zeros(1,length(Freqs));
        for k = 1:length(Freqs)
            idxf(k) = find(freqs==Freqs(k));
            if isempty(idxf(k))
                error('Spatial frequency not found');
            end
        end
        % find indices of specified wavelength range
        idxw = wavelengths > WavelengthsBounds(1) & wavelengths < WavelengthsBounds(2);
        nWavelengths = sum(idxw);
    
        nr = size(rawData,1);
        nc = size(rawData,2);
        %Raw data structure is row, column, spatial frequency, phase,
        %wavelength)
        Raw = zeros(nr,nc,length(idxf),3,nWavelengths);
        % Reshape Data Structure
        for j = 1:length(idxf)
            Raw(:,:,j,:,:) = rawData(:,:,[idxf(j)*3-2,idxf(j)*3-1,idxf(j)*3],idxw);
        end
        dataSet{i}.Raw = Raw;
        dataSet{i}.Profile = rawProfile;
        dataSet{i}.Name = filenames(i);
        dataSet{i}.Wavelengths = wavelengths(idxw);
        dataSet{i}.Exposure = 1000*exposureTimes(idxw);
        dataSet{i}.Freqs = Freqs;
        dataSet{i}.Phases = [0 2*pi/3 4*pi/3];
    end
    phantomID = input('Enter calibration phantom ID from the list:    ');
end 