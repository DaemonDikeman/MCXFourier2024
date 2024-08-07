%% Select Point ROI
function dataSet = SelectROI2(dataSet, phantomID, ROIsize)
    sampleIdx = 1:length(dataSet);
    sampleIdx(phantomID) = [];
    figure()
    N = ceil(sqrt(length(dataSet)-1));
    count = 0;
    message = msgbox('Select a point in the middle of desired ROI in each image');
    for i = sampleIdx
        count = count+1;
        Image = dataSet{i}.MuaCrop1(:,:,4);
        subplot(N,N,count);
        imagesc(Image);
        title(['Mua Image to Draw ROI' ,num2str(count)]);
        pointROI = drawpoint();
        xposROI = round(pointROI.Position(1));
        yposROI = round(pointROI.Position(2));
        xposROImin = round(xposROI-ROIsize/2);
        yposROImin = round(yposROI-ROIsize/2);
        dataSet{i}.ROIpos2 = [xposROImin, yposROImin, ROIsize-1, ROIsize-1];  
    
%          for j = 1:length(lambdas)
%             ImageInputMua = dataSet{i}.Mua(:,:,j);
%             ImageInputMusp = dataSet{i}.Musp(:,:,j);
%             dataSet{i).MuaCrop(:,:,j) = imcrop(ImageInputMua,[dataSet{i}.ROIpos]);
%             dataSet{i}.MuspCrop(:,:,j) = imcrop(ImageInputMusp,[dataSet{i}.ROIpos]);
%          end
            
        end
        

    
    