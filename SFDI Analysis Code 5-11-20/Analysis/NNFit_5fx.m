%% LUT interpolation
function dataSet = NNFit_5fx(dataSet,phantomID)
%error checking
    fx = dataSet{phantomID}.Freqs;
    Fx = [0,0.05,0.1,0.15,0.2];
    try
        if ~all(fx == Fx)
            error('fx must be [0, 0.05, 0.1, 0.15, 0.2]')
        end %if
    catch
        error('wrong number of fx specified')
    end %try
    
    % processing
    ids = 1:length(dataSet);
    ids(phantomID) = [];
    count = 0;
    for i = ids
        count = count + 1;
        disp(['Working on ', num2str(count), ' of ', num2str(length(ids)),' ...']);
        Rd = permute(dataSet{i}.Rd,[1,2,4,3]);
        
        input = [reshape(Rd(:,:,:,1),[],1),reshape(Rd(:,:,:,2),[],1),...
            reshape(Rd(:,:,:,3),[],1),reshape(Rd(:,:,:,4),[],1),...
            reshape(Rd(:,:,:,5),[],1)];
        
        output = OPNet_5fx(input'); % adapted from zhao et al. 2018
        
        badMua = 0.005 > output(1,:) | 0.2 < output(1,:);
        badMusp = 0.1 > output(2,:) | 3.58 < output(2,:);
        output(1,badMua) = nan;
        output(2,badMusp) = nan;
        dataSet{i}.Mua = reshape(output(1,:),size(Rd(:,:,:,1)));
        dataSet{i}.Musp = reshape(output(2,:),size(Rd(:,:,:,1)));
    end
    disp('Done');
end %function