%% Generate angle maps for data set with height maps
function dataSet = genAngleMap(dataSet,phantomID)
    N = ceil(sqrt(length(dataSet)-1));
    figure;
    idx = 1:length(dataSet);
    idx(phantomID) = [];
    count = 0;
    for i = idx
        count = count + 1;
        dataSet{i}.AngleMap = atan(imgradient(dataSet{i}.HeightMap));
        subplot(N,N,count);
        imagesc(dataSet{i}.AngleMap);
        title(['AngleMap ',num2str(count)]);
    end
end %function