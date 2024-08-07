%% Demodualte Data Set
% This function demodulates the data set and corrects the pixel values
% according to camera exposure time.
function dataSet = demodDataSet(dataSet)
    % get size information from the data set
    [Ny,Nx,Nfreqs,~,Nwavelengths] = size(dataSet{1}.Raw);
    disp('Demodulating ...');
    for i = 1:length(dataSet)
        %Demodulate and correct for exposure times
        Demodulated = demodAC(reshape(dataSet{i}.Raw(:,:,:,1,:),1,[]),...
            reshape(dataSet{i}.Raw(:,:,:,2,:),1,[]),...
            reshape(dataSet{i}.Raw(:,:,:,3,:),1,[]));
        dataSet{i}.Demodulated = reshape(Demodulated,Ny,Nx,Nfreqs,Nwavelengths);
        for k = 1:Nwavelengths
            dataSet{i}.Demodulated(:,:,:,k) = dataSet{i}.Demodulated(:,:,:,k)...
                /dataSet{i}.Exposure(k);
        end
        n = 10;
        pause(n)
        disp(['Data Set ', num2str(i), ' Demodulated. Waiting ', num2str(n), ' seconds for memory management purposes.']);
    end
    disp('Done');
end