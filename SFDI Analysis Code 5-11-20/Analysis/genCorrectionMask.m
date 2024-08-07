%% Generate a mask to run height and angle corrections
function dataSet = genCorrectionMask(dataSet,phantomID,threshold)
    figure;
    N = ceil(sqrt(length(dataSet)-1));
    idx = 1:length(dataSet);
    idx(phantomID) = [];
    count = 0;
    for i = idx
        count = count+1;
        planar = dataSet{i}.Raw(:,:,1,1,1);
        mx = max(max(planar));
        mn = min(min(planar));
        planar = (planar - mn)./(mx-mn);
        mask = planar < threshold;
        mask = imfill(mask,'holes');
        mask = imfill(~mask,'holes');
        dataSet{i}.Mask = mask;
        subplot(N,N,count);
        imagesc(mask);
        title(['Correction Mask ',num2str(count)]);
    end
    
end %function