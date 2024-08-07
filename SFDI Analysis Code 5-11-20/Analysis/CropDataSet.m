%% Select Point ROI
function dataSet = CropDataSet(dataSet, phantomID, ROIsize, lambdas)
    sampleIdx = 1:length(dataSet);
    sampleIdx(phantomID) = [];
    figure()
    N = ceil(sqrt(length(dataSet)-1));
    count = 0;
    message = msgbox('Select a point in the middle of desired ROI in each image');
    for i = sampleIdx
        count = count+1;
        Image = dataSet{i}.Musp(:,:,1);
        subplot(N,N,count);
        imagesc(Image);
        title(['Musp Image to Draw ROI' ,num2str(count)]);
        pointROI = drawpoint();
        xposROI = round(pointROI.Position(1));
        yposROI = round(pointROI.Position(2));
        xposROImin = round(xposROI-ROIsize/2);
        yposROImin = round(yposROI-ROIsize/2);
        dataSet{i}.ROIpos = [xposROImin, yposROImin, ROIsize, ROIsize];  
    
%          for j = 1:length(lambdas)
%             ImageInputMua = dataSet{i}.Mua(:,:,j);
%             ImageInputMusp = dataSet{i}.Musp(:,:,j);
%             dataSet{i).MuaCrop(:,:,j) = imcrop(ImageInputMua,[dataSet{i}.ROIpos]);
%             dataSet{i}.MuspCrop(:,:,j) = imcrop(ImageInputMusp,[dataSet{i}.ROIpos]);
%          end
            
        end
        

    
    