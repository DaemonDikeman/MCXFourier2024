%% ROI average 2 Chromophore Extraction 
function dataSet = ROIaverage2Chromophores(dataSet,phantomID,lambdas,HbO2spectrum, Hhbspectrum, Ny, Nx);

sampleIdx = 1:length(dataSet); 
sampleIdx(phantomID) = [];

count = 0;

for i = sampleIdx
    count = count + 1;
    Mua = dataSet{i}.MuaCrop;
    MuaSelect = Mua(:,:,4:end);
    mua_measured = reshape(MuaSelect, [], length(lambdasSelect));
    NaNcount = sum(isnan(reshape(mua_measured,[],1)))/numel(mua_measured);
    mua_measuredSize = size(mua_measured);
    numMua_measured = mua_measuredSize(1);
    numPixels = Nx * Ny;
    
    HbO2 = zeros(1,numPixels);
    Hhb = zeros(1,numPixels); 
       
    disp('Fitting for Tissue Chromophores')
    
    for k = 1:size(mua_measured,1);
              mua_indexed = mua_measured(k,:);
            %size(mua_measured);  
               if any(isnan(mua_indexed));
                   
                    HbO2(k) = NaN;
                    Hhb(k) = NaN;
            
               else
                [coeff, fit] = spectrumFit2(HbO2spectrum, Hhbspectrum, mua_indexed);
                dataSet{k}.HbO2val = coeff(1,1);
                dataSet{k}.HHbval = coeff(1,2);
               end
               
        
    end

end
