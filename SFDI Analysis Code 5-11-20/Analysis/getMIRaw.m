%% GET MI RAW
% Inputs
    % directory (str): directory that contains the raw file with SFDI scan data
    % fileIndex (int): fileIndex specifies which raw file to import within the directory

% Outputs
    % data (4D int): pixel data from SFDI scan
    % exposureTimes (1D int): camera exposure times for each wavelength (ms)
    % freqs (1D double): spatial frequncies used for scan (1/mm)

function [rawData,rawProfile,exposureTimes,wavelengths,freqs] = getMIRaw(directory,fileIndex)
%Default camera binning for MI system is 2
bin = 2;
%bin = 1;

%Read raw data from raw file as a 1D array
fnames = dir(fullfile(directory,'*_raw'));
fid = fopen(fullfile(directory,fnames(fileIndex).name),'r');
rawData = fread(fid,'single');
fclose(fid);

%Read raw profile data from raw file as a 1D array
fnames = dir(fullfile(directory,'*_profile'));
fid = fopen(fullfile(directory,fnames(fileIndex).name),'r');
rawProfile = fread(fid,'single');
fclose(fid);

%Gather scan parameters from directory:
fnames = dir(fullfile(directory,'*_parameters.xml'));
fid = fopen(fullfile(directory,fnames(fileIndex).name),'r');
scanParams = fscanf(fid,'%s');
fclose(fid);

%Camera exposure times for each wavelength (1D int)
exposureTimes = tag2array(scanParams,'ExposureTime');
%Scan Wavelengths (1D int)
wavelengths = tag2array(scanParams,'Wavelengths');
%Scan frequencies (1D double)
freqs = tag2array(scanParams,'SpatialFrequencies');

%Reshape linear image data to 4D array (Ny,Nx,Nphases*Nfrequencies,Nwavelengths)
rawData = reshape(rawData,1392/bin,1040/bin,3*length(freqs),length(wavelengths));
%Reshape linear profile image data to 3D array (Ny,Nx,nPhases,NFrequencies)
%rawProfile = reshape(rawProfile,1392/bin,1040/bin,3,3);
rawProfile = reshape(rawProfile,1392/bin,1040/bin,3,[]);
end