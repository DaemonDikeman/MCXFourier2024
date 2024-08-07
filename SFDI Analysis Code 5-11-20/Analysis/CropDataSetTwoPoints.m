%% Select Point ROI
function dataSet = CropDataSetTwoPoints(dataSet, phantomID, lambdas)
    sampleIdx = 1:length(dataSet);
    sampleIdx(phantomID) = [];
    figure()
    N = ceil(sqrt(length(dataSet)-1));
    count = 0;
    message = msgbox('Select a point in the middle of desired ROI in each image');
    for i = sampleIdx
        count = count+1;
        Image = dataSet{i}.Mua(:,:,4);
        subplot(N,N,count);
        imagesc(Image);
        title(['Mua Image to Draw ROI' ,num2str(count)]);
        pointROI = drawpoint();
        pointROI2 = drawpoint();
        xposROImax = round(pointROI.Position(1));
        yposROImax = round(pointROI.Position(2));
        xposROImin = round(pointROI2.Position(1));
        yposROImin = round(pointROI2.Position(2));
        dataSet{i}.ROIpos = [min([xposROImin xposROImax]), min([yposROImin yposROImax]), abs(xposROImax-xposROImin), abs(yposROImax-yposROImin)];  
    
%          for j = 1:length(lambdas)
%             ImageInputMua = dataSet{i}.Mua(:,:,j);
%             ImageInputMusp = dataSet{i}.Musp(:,:,j);
%             dataSet{i).MuaCrop(:,:,j) = imcrop(ImageInputMua,[dataSet{i}.ROIpos]);
%             dataSet{i}.MuspCrop(:,:,j) = imcrop(ImageInputMusp,[dataSet{i}.ROIpos]);
%          end
            
        end
        

    
    