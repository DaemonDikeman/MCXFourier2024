%% Calculates ROI Average Mua
function dataSet = ROIaverageMua(dataSet,phantomID,lambdas)

sampleIdx = 1:length(dataSet); % for whatever reason, I get errors with 1
sampleIdx(phantomID) = [];

avgMua = zeros(length(sampleIdx),length(lambdas));
stdevMua = zeros(length(sampleIdx),length(lambdas));
avgMusp = zeros(length(sampleIdx),length(lambdas));
stdevMusp = zeros(length(sampleIdx),length(lambdas));
avgRd0 = zeros(length(sampleIdx),length(lambdas));
stdevRd0 = zeros(length(sampleIdx),length(lambdas));
avgRd1 = zeros(length(sampleIdx),length(lambdas));
stdevRd1 = zeros(length(sampleIdx),length(lambdas));

for i = sampleIdx
    for lambda = 1:length(lambdas)
        MuaCrop = dataSet{i}.MuaCrop(:,:,lambda);
        MuspCrop = dataSet{i}.MuspCrop(:,:,lambda);
        Rd0Crop = dataSet{i}.Rd0Crop(:,:,lambda);
        Rd1Crop = dataSet{i}.Rd1Crop(:,:,lambda);
        dataSet{i}.avgMua(lambda) = mean(mean(MuaCrop,'omitnan'),'omitnan');
        dataSet{i}.stdevMua(lambda) = std(std(MuaCrop, 'omitnan'), 'omitnan');
        dataSet{i}.avgMusp(lambda) = mean(mean(MuspCrop,'omitnan'),'omitnan');
        dataSet{i}.stdevMusp(lambda) = std(std(MuspCrop, 'omitnan'), 'omitnan');
        dataSet{i}.avgRd0(lambda) = mean(mean(Rd0Crop,'omitnan'),'omitnan');
        dataSet{i}.stdevRd0(lambda) = std(std(Rd0Crop, 'omitnan'), 'omitnan');
        dataSet{i}.avgRd1(lambda) = mean(mean(Rd1Crop,'omitnan'),'omitnan');
        dataSet{i}.stdevRd1(lambda) = std(std(Rd1Crop, 'omitnan'), 'omitnan');
    end
end
