%% Crop Images using Point ROI
function dataSet = cropImage2(dataSet,phantomID,lambdas)
sampleIdx = 1:length(dataSet);
    sampleIdx(phantomID) = [];
    for i = sampleIdx
         for j = 1:length(lambdas)
            ImageInputMua = dataSet{i}.MuaCrop1(:,:,j);
            ImageInputMusp = dataSet{i}.MuspCrop1(:,:,j);
            RawImage971 = dataSet{i}.Raw971Crop1(:,:,9);
            dataSet{i}.MuaCrop2(:,:,j) = imcrop(ImageInputMua,[dataSet{i}.ROIpos2]);
            dataSet{i}.MuspCrop2(:,:,j) = imcrop(ImageInputMusp,[dataSet{i}.ROIpos2]);
            dataSet{i}.Raw971Crop2(:,:,9) = imcrop(RawImage971, [dataSet{i}.ROIpos2]);
         end
    end
    