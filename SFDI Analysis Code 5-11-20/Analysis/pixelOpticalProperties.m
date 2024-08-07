%% Calculate pixel optical properties within selected ROI

function dataSet = pixelScatteringProps(dataSet,phantomID,lambdas,lambdasSelect, Ny, Nx);

sampleIdx = 2:length(dataSet);

%depthLambdas = length(lambdas);

count = 0;

    for i = sampleIdx   
        
        disp('Fitting A and bMie')
        count = count +1;
       
        Musp = dataSet{i}.MuspCrop;
        
        MuspSelect = Musp(:,:,4:end);
      
     
        musp_measured = reshape(MuspSelect,[],length(lambdasSelect));
          
        musp_measuredSize = size(musp_measured);
        numMusp_measured = musp_measuredSize(1);
          
        disp(['Working on ', num2str(count), ' of ', num2str(length(sampleIdx)),' ...']);
          
          for k = 1:length(numMusp_measured);
              
              A = zeros(1,numMusp_measured);
              bMie = zeros(1,numMusp_measured);
              
              % This if statement is far too restrictive, I want to just
              % ignore those pixels. Right now it's making the entire image
              % NaN. Otherwise it works. 
              
               if any(isnan(musp_measured(k,:)));
                   
                    A(k) = NaN;
                    bMie(k) = NaN;
            
               else
              
                    musp_measured = musp_measured(k,:);
              
        
                   [a(k),bMie(k),fit(k)] = scatFit1(lambdasSelect,musp_measured);
                
               end
          end
          
          
                Aval = reshape(A, [Ny,Nx]);
                bMieVal = reshape(bMie, [Ny,Nx]);
                
                dataSet{i}.A = Aval;
                dataSet{i}.bMie = bMieVal;
    end
end

    







