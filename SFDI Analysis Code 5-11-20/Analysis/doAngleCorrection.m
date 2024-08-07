%% Angle corrections for a data set with height corrections
function dataSet = doAngleCorrection(dataSet,phantomID,smoothStd)
    idx = 1:length(dataSet);
    idx(phantomID) = [];
    for i = idx
        AngleMap = imgaussfilt(cos(dataSet{i}.AngleMap),smoothStd);
        for j = 1:length(dataSet{i}.Freqs)
            for k = 1:length(dataSet{i}.Wavelengths)   
                dataSet{i}.AngleCorrected(:,:,j,k) = dataSet{i}.HeightCorrected(:,:,j,k)./AngleMap;
            end
        end
    end
end %function