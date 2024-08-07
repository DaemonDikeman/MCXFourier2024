%% Crop Images using Point ROI
function dataSet = cropImage(dataSet,phantomID,lambdas)
sampleIdx = 1:length(dataSet);
    sampleIdx(phantomID) = [];
    for i = sampleIdx
         for j = 1:length(lambdas)
            ImageInputMua = dataSet{i}.Mua(:,:,j);
            ImageInputMusp = dataSet{i}.Musp(:,:,j);
            ImageInputRd00 = dataSet{i}.Rd(:,:,1,j);
            ImageInputRd01 = dataSet{i}.Rd(:,:,3,j);
            dataSet{i}.MuaCrop(:,:,j) = imcrop(ImageInputMua,[dataSet{i}.ROIpos]);
            dataSet{i}.MuspCrop(:,:,j) = imcrop(ImageInputMusp,[dataSet{i}.ROIpos]);
            dataSet{i}.Rd0Crop(:,:,j) = imcrop(ImageInputRd00, [dataSet{i}.ROIpos]);
            dataSet{i}.Rd1Crop(:,:,j) = imcrop(ImageInputRd01, [dataSet{i}.ROIpos]);
         end
    end
    