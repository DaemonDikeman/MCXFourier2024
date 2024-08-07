%% Generate heigth map for a data set with a correction mask
function dataSet = genHeightMap(dataSet,phantomID,fxID,heightFactor,smoothStd)
    N = ceil(sqrt(length(dataSet)-1));
    % reference phase map
    s = dataSet{phantomID}.Profile;
    I1 = s(:,:,1,fxID);
    I2 = s(:,:,2,fxID);
    I3 = s(:,:,3,fxID);
    I1 = imgaussfilt(I1,smoothStd);
    I2 = imgaussfilt(I2,smoothStd);
    I3 = imgaussfilt(I3,smoothStd);
    wrapppedPhase = atan2((sqrt(3)*(I1-I3)),(2*I2 - I1 - I3));
    phase = unwrap(wrapppedPhase')';
    dataSet{phantomID}.PhaseMap = phase;

    % height maps
    figure;
    idx = 1:length(dataSet);
    idx(phantomID) = [];
    count = 0;
    for i = idx
        count = count + 1;
        s = dataSet{i}.Profile;
        I1 = s(:,:,1,fxID);
        I2 = s(:,:,2,fxID);
        I3 = s(:,:,3,fxID);
        I1 = imgaussfilt(I1,smoothStd);
        I2 = imgaussfilt(I2,smoothStd);
        I3 = imgaussfilt(I3,smoothStd);
        wrapppedPhase = atan2((sqrt(3)*(I1-I3)),(2*I2 - I1 - I3));
        wrapppedPhase(~dataSet{i}.Mask) = nan;
        wrapppedPhase = fillmissing(wrapppedPhase,'nearest');
        phase = unwrap(wrapppedPhase')';
        phase(~dataSet{i}.Mask) = nan;
        dataSet{i}.HeightMap = heightFactor*(phase-dataSet{phantomID}.PhaseMap);
        subplot(N,N,count);
        imagesc(dataSet{i}.HeightMap,[-5,5]);
        title(['Height Map ',num2str(count)]);
    end
end %function