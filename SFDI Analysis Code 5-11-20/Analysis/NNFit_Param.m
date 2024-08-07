%% LUT interpolation
function dataSet = NNFit_Param(dataSet,phantomID,parameter)
%error checking
    fx = dataSet{phantomID}.Freqs;
    lambdas = dataSet{phantomID}.Wavelengths;
    Lambdas = [621   659   691   731   851];
    Fx = [0,0.05,0.1,0.15,0.2];
    try
        if ~all(fx == Fx)
            error('fx must be [0, 0.05, 0.1, 0.15, 0.2]')
        end %if
    catch
        error('wrong number of fx specified')
    end %try
    
    idxs = nan(1,length(Lambdas));
    for i = 1:length(Lambdas)
        idx = find(lambdas==Lambdas(i));
        if ~isempty(idx)
            idxs(i) = idx;
        else
            error(['wavelength ' num2str(Lambdas(i)) ' not found in data set']);
        end
    end
    
    
    
    
    
    
    
    % processing
    ids = 1:length(dataSet);
    ids(phantomID) = [];
    count = 0;
    for i = ids
        count = count + 1;
        disp(['Working on ', num2str(count), ' of ', num2str(length(ids)),' ...']);
        Rd = permute(dataSet{i}.Rd(:,:,:,idxs),[1,2,4,3]);
        if strcmp(parameter,'A')
            output = Anet(reshape(Rd,[],25)');
            badParam = 0 > output(1,:) | 4 < output(1,:);
            output(badParam) = nan;
            dataSet{i}.A = reshape(output,size(dataSet{i}.Rd(:,:,1,1)));
        elseif strcmp(parameter,'B')
            output = Bnet(reshape(Rd,[],25)');
            badParam = 0 > output(1,:) | 4 < output(1,:);
            output(badParam) = nan;
            dataSet{i}.B = reshape(output,size(dataSet{i}.Rd(:,:,1,1)));
        elseif strcmp(parameter,'HbO2')
            output = HbO2Net(reshape(Rd,[],25)');
            badParam = 0 > output(1,:) | 0.1 < output(1,:);
            output(badParam) = nan;
            dataSet{i}.HbO2 = reshape(output,size(dataSet{i}.Rd(:,:,1,1)));
        elseif strcmp(parameter,'HHb')
            output = HHbNet(reshape(Rd,[],25)');
            badParam = 0 > output(1,:) | 0.1 < output(1,:);
            output(badParam) = nan;
            dataSet{i}.HHb = reshape(output,size(dataSet{i}.Rd(:,:,1,1)));
        end

    end
    disp('Done');
end %function