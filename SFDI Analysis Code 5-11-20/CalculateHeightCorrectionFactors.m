%% Calculate Height Correction Factors
[dataSet,phantomID] = importDataSet([0,0.05,0.1,0.15,0.2],[0,2000],'MI');
dataSet = demodDataSet(dataSet);
%% Get phase maps
profile1 = getProfile(dataSet{phantomID}.Profile,dataSet{3}.Profile);
profile2 = getProfile(dataSet{phantomID}.Profile,dataSet{5}.Profile);
%% Apply height factor to get absolute height
s = -0.9745;
H1 = mean2(s*profile1(250:450,180:340));
H2 = mean2(s*profile2(250:450,180:340));
%% Estimate height correction factors
nFx = length(dataSet{phantomID}.Freqs);
nWvs = length(dataSet{phantomID}.Wavelengths);
MAC1= nan(nFx,nWvs);
MAC2 = nan(nFx,nWvs);
for fx = 1:nFx
    for wv = 1:nWvs
        MAC1(fx,wv) = mean2(dataSet{3}.Demodulated(250:450,180:340,fx,wv));
        MAC2(fx,wv) = mean2(dataSet{5}.Demodulated(250:450,180:340,fx,wv));
    end
end
HCFactors.Factors = (MAC2 - MAC1)./(H2 - H1);
HCFactors.Wavelengths = dataSet{phantomID}.Wavelengths;
HCFactors.Freqs = dataSet{phantomID}.Freqs;
%% Save height correction factors
fileName = 'HCFactors';
save(fileName,HCFactors);