%% Calculates ROI Average Mua
function dataSet = ROIaverageScatteringProps(dataSet,phantomID,lambdas);


sampleIdx = 2:length(dataSet); % In this case 1 is calibration but will this always hold? 


avgMua = zeros(length(sampleIdx),length(lambdas));
stdevMua = zeros(length(sampleIdx), length(lambdas));
avgMusp = zeros(length(sampleIdx),length(lambdas));
stdevMusp = zeros(1, length(lambdas));

for i = sampleIdx;
%     if sampleIdx(phantomID) == []
%     dataSet{sampleIdx}.Mua = NaN;
%     dataSet{sampleIdx}.Musp = NaN;
%     end
    for lambda= 1:length(lambdas);
        MuaCrop = dataSet{i}.MuaCrop(:,:,lambda);
        MuspCrop = dataSet{i}.MuspCrop(:,:,lambda);
        dataSet{i}.avgMua(lambda) = mean(mean(MuaCrop,'omitnan'),'omitnan');
        dataSet{i}.stdevMua(lambda) = nanstd(nanstd(MuaCrop));
        dataSet{i}.avgMusp(lambda) = mean(mean(MuspCrop,'omitnan'),'omitnan');
        dataSet{i}.stdevMusp(lambda) = nanstd(nanstd(MuspCrop));
    end
end
