%% Calculate pixel optical properties within selected ROI

function dataSet = pixel2Chromophores(dataSet,phantomID,lambdas, lambdasSelect, HbO2spectrum,Hhbspectrum, Nx, Ny);

sampleIdx = 1:length(dataSet);
sampleIdx(phantomID) = [];


MuaMeasured = zeros(Ny,Nx,length(lambdasSelect));

%depthLambdas = length(lambdas);
count = 0;
for i = sampleIdx    
    
    disp('Fitting Oxy and Deoxy Chromophores')
    count = count +1;
        
        Mua = dataSet{i}.MuaCrop;
        
         MuaSelect = Mua(:,:,4:end);
         MuaMeasured = reshape(MuaSelect,[],length(lambdasSelect));
         mua_measuredSize = size(MuaMeasured);
         numMua_measured = mua_measuredSize(1);
         NumPixels = Ny*Nx;
         HbO2 = zeros(1,NumPixels);
         Hhb = zeros(1,NumPixels);
     
          
       disp(['Working on ', num2str(count), ' of ', num2str(length(sampleIdx)),' ...']);
         
       for k = 1:size(MuaMeasured,1);
           
           Mua_indexed = MuaMeasured(k,:);
           
           size(MuaMeasured);
               if any(isnan(Mua_indexed));
              
               % HbO2 has 1 value at the 1st pixel, the rest is NaN. A is
              % all NaN. It is something wrong with the for loop k =
              % 1:length(numMusp_measured)
                   
                    HbO2(k) = NaN;
                    Hhb(k) = NaN;
                   
            
               else
              
                    MuaMeasured = MuaMeasured(k,:);
              
        
                   [coeff, fit] = spectrumFit3(HbO2spectrum,Hhbspectrum,Mua_indexed);
                    
                       HbO2(k) = coeff(1,1);
                       Hhb(k) = coeff(1,2);
                       
               end
       end
          
          HbO2 = reshape(HbO2, [Ny,Nx]);
          Hhb = reshape(Hhb, [Ny,Nx]); 
          
          dataSet{i}.HbO2 = HbO2;
          dataSet{i}.Hhb = Hhb;
          
end
end


