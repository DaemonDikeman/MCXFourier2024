%% LUT interpolation
function dataSet = lutFit(dataSet,phantomID,fx)
    if length(fx) ~= 2
        error('more than 2 fx specified in data set')
    end
    LUT = load('LUT_Homogeneous');
    idxs = nan(1,2);
    for i = 1:length(fx)
        idx = find(LUT.LUT.Fx==fx(i));
        if ~isempty(idx)
            idxs(i) = idx;
        else
            error(['spatial frequency ' num2str(fx(i)) ' not found in LUT']);
        end
    end
    ids = 1:length(dataSet);
    ids(phantomID) = [];
    count = 0;
    for i = ids
        count = count + 1;
        disp(['Working on ', num2str(count), ' of ', num2str(length(ids)),' ...']);
        Rd = permute(dataSet{i}.Rd,[1,2,4,3]);
        Mua = griddata(LUT.LUT.Rd(:,:,idxs(1)),LUT.LUT.Rd(:,:,idxs(2)),...
            LUT.LUT.Mua,reshape(Rd(:,:,:,idxs(1)),[],1),reshape(Rd(:,:,:,idxs(2)),[],1));
        Musp = griddata(LUT.LUT.Rd(:,:,idxs(1)),LUT.LUT.Rd(:,:,idxs(2)),...
            LUT.LUT.Musp,reshape(Rd(:,:,:,idxs(1)),[],1),reshape(Rd(:,:,:,idxs(2)),[],1));
        dataSet{i}.Mua = reshape(Mua,size(Rd(:,:,:,1)));
        dataSet{i}.Musp = reshape(Musp,size(Rd(:,:,:,1)));
    end
    disp('Done');
end %function