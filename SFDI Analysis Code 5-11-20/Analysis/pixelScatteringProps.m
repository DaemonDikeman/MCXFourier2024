%% Calculate pixel optical properties within selected ROI

function dataSet = pixelScatteringProps(dataSet,phantomID,lambdas,lambdasSelect, Ny, Nx);

sampleIdx = 1:length(dataSet);
sampleIdx(phantomID) = [];

%depthLambdas = length(lambdas);

count = 0;

    for i = sampleIdx   
        
        disp('Fitting A and bMie')
        count = count +1;
        Musp = dataSet{i}.MuspCrop1;
        MuspSelect = Musp(:,:,4:end);
        musp_measured = reshape(MuspSelect,[],length(lambdasSelect));
        NaNcount = sum(isnan(reshape(musp_measured,[],1)))/numel(musp_measured);
        musp_measuredSize = size(musp_measured);
        numMusp_measured = musp_measuredSize(1);
        numPixels = Nx * Ny;
        
         A = zeros(1,numPixels);
         bMie = zeros(1,numPixels);
         
         disp(['Working on ', num2str(count), ' of ', num2str(length(sampleIdx)),' ...']);
          
          for k = 1:size(musp_measured,1)
              musp_indexed = musp_measured(k,:);
            %size(musp_measured);  
               if any(isnan(musp_indexed))
                   
                    A(k) = NaN;
                    bMie(k) = NaN;
            
               else
                   [A(k),bMie(k),fit(k)] = scatFit1(lambdasSelect,musp_indexed);
        
               end
          end
                Aval = reshape(A, [Ny,Nx]);
                bMieVal = reshape(bMie, [Ny,Nx]);
                
                dataSet{i}.A = Aval;
                dataSet{i}.bMie = bMieVal;
    end
    disp('Done')
end