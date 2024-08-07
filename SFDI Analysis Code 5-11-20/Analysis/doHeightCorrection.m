%% Run height correction on data set with height maps
function dataSet = doHeightCorrection(dataSet,phantomID,correctionFactorFileName)
    CF = load(correctionFactorFileName);
    idx = 1:length(dataSet);
    idx(phantomID) = [];
    for i = idx
        for j = 1:length(dataSet{i}.Freqs)
            fx = CF.HCFactors.Freqs == dataSet{i}.Freqs(j);
            if all(~fx)
                error(['correction factor not found at fx = ' num2str(dataSet{i}.Freqs(j))]);
            end
            for k = 1:length(dataSet{i}.Wavelengths)   
                wv = CF.HCFactors.Wavelengths == dataSet{i}.Wavelengths(k);
                if all(~wv)
                    error(['correction factor not found at lambda = ' num2str(dataSet{i}.Wavelengths(k))]);
                end
                HCF = CF.HCFactors.Factors(fx,wv);
                dataSet{i}.HeightCorrected(:,:,j,k) = dataSet{i}.Demodulated(:,:,j,k) + HCF*dataSet{i}.HeightMap;
            end
        end
    end
end %function