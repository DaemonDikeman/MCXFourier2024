%% Crop Images using Point ROI
function dataSet = cropImage1(dataSet,phantomID,lambdas)
sampleIdx = 1:length(dataSet);
    sampleIdx(phantomID) = [];
    for i = sampleIdx
         for j = 1:length(lambdas)
            ImageInputMua = dataSet{i}.Mua(:,:,j);
            ImageInputMusp = dataSet{i}.Musp(:,:,j);
            RawImage971 = dataSet{i}.Raw(:,:,1,1,9);
            dataSet{i}.MuaCrop1(:,:,j) = imcrop(ImageInputMua,[dataSet{i}.ROIpos1]);
            dataSet{i}.MuspCrop1(:,:,j) = imcrop(ImageInputMusp,[dataSet{i}.ROIpos1]);
            dataSet{i}.Raw971Crop1(:,:,9) = imcrop(RawImage971, [dataSet{i}.ROIpos1]);
         end
    end
    