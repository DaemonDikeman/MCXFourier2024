%% import sfdi data and format the data set structure
function [dataSet,phantomID] = importDataSet(Freqs,WavelengthsBounds,format)
%Check format string
if strcmp(format,'MI')
    fileIdx = 1; % 1 is MI system
    try % to open a directory to get scan files
        directory = uigetdir('Data_Sets','Select an MI data set');
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
    
    %Check that the data set to be imported is from openSFDI system   
    elseif strcmp(format,'OS')
        try
        directory = uigetdir('Data_Sets','Select an OpenSFDI data set');
        files = dir([directory '\*.mat']);
        catch
            dataSet = NaN;
            phantomID = NaN;
           return 
        end
        filenames = {files.name};
        dataSet = cell(length(filenames),1);
        for i = 1:length(filenames)
            disp(['Loading ',num2str(i), ' out of ', num2str(length(dataSet)),' ... ']);
            file = load([directory '\' filenames{i}]);
            file.S.Name = erase(filenames{i},'.mat');
            dataSet{i} = file.S;
        end
else
    error('Incorrect format type.')
end
    Names = string(filenames);
    for i = 1:length(Names)
        disp(strcat(num2str(i),': ',Names(i)));
    end
        phantomID = input('Enter calibration phantom ID from the list:    ');
end %function