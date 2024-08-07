%% LUT interpolation
function dataSet = lutFit_new(dataSet,phantomID,fx,lutFile)
    if length(fx) ~= 2
        error('more than 2 fx specified in data set')
    end
    LUTa = load(lutFile);
    if isfield(LUTa,'LUT')
        LUT = LUTa.LUT;
    else
        LUT = LUTa;
    end
    if isfield(LUT,'Fx')
        idxs = nan(1,2);
        for i = 1:length(fx)
            idx = find(abs(LUT.Fx-fx(i)) < 0.001);
            if ~isempty(idx)
                idxs(i) = idx;
            else
                error(['spatial frequency ' num2str(fx(i)) ' not found in LUT']);
            end
        end
%     elseif isfield(LUT,'M3')
%         idxs = nan(1,2);
%         for i = 1:length(fx)
%             prompt = "What field should be used for frequency " + num2str(fx(i)) + "?";
%             field = input(prompt,"s");
%             idx = find(LUT.(field));
%             if ~isempty(idx)
%                 idxs(i) = idx;
%             else
%                 error(['field ' field ' not found in LUT']);
%             end
%         end
    elseif fx ~= [0,0.1]
        disp(['Warning, LUT does not specify frequencies; frequencies ' num2str(fx(1)) ' and ' num2str(fx(2)) ' are atypical']);
        idxs = nan(1,2);
        for i = 1:length(fx)
            idx = find(dataSet{1}.Freqs==fx(i));
            if ~isempty(idx)
                idxs(i) = idx;
            else
                error(['spatial frequency ' num2str(fx(i)) ' not found in LUT']);
            end
        end
    else
        idxs = nan(1,2);
        for i = 1:length(fx)
            idx = find(dataSet{1}.Freqs==fx(i));
            if ~isempty(idx)
                idxs(i) = idx;
            else
                error(['spatial frequency ' num2str(fx(i)) ' not found in LUT']);
            end
        end
    end
    
    ids = 1:length(dataSet);
    ids(phantomID) = [];
    count = 0;
    for i = ids
        count = count + 1;
        disp(['Working on ', num2str(count), ' of ', num2str(length(ids)),' ...']);
        Rd = permute(dataSet{i}.Rd,[1,2,4,3]);
        if isfield(LUT,'Rd')
            Mua = griddata(double(LUT.Rd(:,:,idxs(1))),double(LUT.Rd(:,:,idxs(2))),double(LUT.Mua),...
                reshape(Rd(:,:,:,idxs(1)),[],1),reshape(Rd(:,:,:,idxs(2)),[],1));
            Musp = griddata(double(LUT.Rd(:,:,idxs(1))),double(LUT.Rd(:,:,idxs(2))),double(LUT.Musp),...
                reshape(Rd(:,:,:,idxs(1)),[],1),reshape(Rd(:,:,:,idxs(2)),[],1));
        elseif isfield(LUT,'M4')
            Mua = griddata(LUT.M1(:,:),LUT.M4(:,:),LUT.Mua,...
                reshape(Rd(:,:,:,idxs(1)),[],1),reshape(Rd(:,:,:,idxs(2)),[],1));
            Musp = griddata(LUT.M1(:,:),LUT.M4(:,:),LUT.Musp,...
                reshape(Rd(:,:,:,idxs(1)),[],1),reshape(Rd(:,:,:,idxs(2)),[],1));
        else
            Mua = griddata(LUT.M1(:,:),LUT.M2(:,:),LUT.Mua,...
                reshape(Rd(:,:,:,idxs(1)),[],1),reshape(Rd(:,:,:,idxs(2)),[],1));
            Musp = griddata(LUT.M1(:,:),LUT.M2(:,:),LUT.Musp,...
                reshape(Rd(:,:,:,idxs(1)),[],1),reshape(Rd(:,:,:,idxs(2)),[],1));
        end
        dataSet{i}.Mua = reshape(Mua,size(Rd(:,:,:,1)));
        dataSet{i}.Musp = reshape(Musp,size(Rd(:,:,:,1)));
    end
    disp('Done');
end %function