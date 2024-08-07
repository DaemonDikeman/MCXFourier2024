%% Generate training and testing data for wavelength-specific tissue parameters
% Define the wavelengths that will be used to collect Rd data
lambdas = [621   659   691   731   851];
fx = [0,0.05,0.1,0.15,0.2];
% Spectra we want to include in the inverse model. The absorption
% coefficients depend on the chromophores we include the model.
spectra = {'HbO2','HHb'};

%Number of random samples we want to generate. The data set will be trimmed
% to exlude parameters that lead to optical properties that are
% out-of-range in the forward model.
N = 100000;

% define the optical property bounds of the forward model
% absorption is only defined between [0.001,0.2] mm^-1
% scattering is only defined between [0.2,5] mm^-1
muaBounds = [0.0005,0.2];
muspBounds = [0.2,3.58];

% Concentration bounds for each chromophore spectrum
CBounds = [0,  0.1;... %lower and upper bounds for HbO2 [mM]
           0,  0.1];   %lower and upper bounds for HHb [mM]


%Scattering Parameters from mu_s' = A*(lambda/500nm)^(-B)
ABounds = [0,4]; % Upper and lower bounds for A
BBounds = [0,4]; % Upper and lower bounds for B

%==========================================================================
%===============================CODE======================================
%=========================================================================
% We can generate random values for each parameter with rand(), which
% returns uniform distribution of values on the open interval (0,1). Values
% are mapped to the intervals we want with... 
% uniformDist(lower,upper,numSamples)
A = uniformDist(ABounds(1),ABounds(2),N); %Random A values
B = uniformDist(BBounds(1),BBounds(2),N); %Random B values

% Here we're generating the absorption parameters from the 'CBounds' matrix
absorption_parameters = zeros(length(spectra),N); %is a preallocated matrix
for i = 1:length(spectra) %for each spectrum we want to fit for...
    %generate the random values and store them in 'absorption_parameters'
    absorption_parameters(i,:) = uniformDist(CBounds(i,1),CBounds(i,2),N);
end

% store all input parameters in params. We will use this to train the 
% inverse model later
params = [absorption_parameters;A;B];

% Now that we have the randomized parameters, we need to convert them to
% optical properties with a wavelength dependent model.
% For absorption, we simply pull mu_a from the superposition of absorption
% spectra at the different concentrations.
%
% For scattering, we calculate mu_s' from the scattering model:
% Mu_s' = A*(lambda/500nm)^(-B)

% We need to keep track of the OPs that fall out of range with a list of
% 'badValues'
badValues = false(1,N);

% preallocate the input array for the 2 optical properties. This will be
% used to feed a forward model to generate diffuse reflectance
input = zeros(length(lambdas),2,N);
for i = 1:N % for each random input parameter combination
    badMuaFlag = false; %set badMuaFlag to false so we can catch bad values
    absorption_spectrum = zeros(length(lambdas),1); % allocate the net mu_a spectrum
    for j = 1:length(spectra) %for each spectrum
        % get absorption values
        mu_a = getChromoData(lambdas,spectra{j})';
        % superimpose absorption spectra
        absorption_spectrum = absorption_spectrum + absorption_parameters(j,i)*mu_a;
        % check that mua values are within bounds
        if any(absorption_spectrum < muaBounds(1) | absorption_spectrum > muaBounds(2))
            %if one of them isn't
            badMuaFlag = true; %set the bad mua value flag true
            break %out of the inner spectrum loop
        end
    end
    % if we catch a bad mua value...
    if badMuaFlag == true
        % don't calculate scattering and record the bad value location
        badValues(i) = true;
    else %if the mua values are good...
        % get mu_s' from scattering model
        scattering_spectrum = getScattering1(lambdas,A(i),B(i))';
        % check for out-of-range musp values
        if any(scattering_spectrum < muspBounds(1) | scattering_spectrum > muspBounds(2))
            % record bad scattering value location
            badValues(i) = true;
        else %if all mua and musp values are in range...
            % load mua and musp into the 'input' tensor
            input(:,:,i) = [absorption_spectrum,scattering_spectrum];
        end
    end
    progressbar(i/N);
end
% Now we need to remove the out-of-range OPs from the input array and
% corresponding parameter matrix
input = input(:,:,~badValues);
params = params(:,~badValues)';
% show how the % of samples we removed
disp(['Removed ',num2str(100*sum(badValues)/N),'% of the data set' ]);
% get Rd values
Rd = zeros(size(input,3),length(fx),length(lambdas));
input = permute(input,[1,3,2]);
LUT = load('LUT_Homogeneous');
for i = 1:length(fx)
    Rd(:,:,i) = griddata(LUT.LUT.Mua,LUT.LUT.Musp,LUT.LUT.Rd(:,:,i),input(:,:,1),input(:,:,2))';
    progressbar(i/length(fx));
end

paramDataSet.Rd = Rd;
paramDataSet.Params = params;
paramDataSet.Fx = fx;
paramDataSet.Lambdas = lambdas;
paramDataSet.Spectra = [spectra,'A','B'];
%% Train Network on generated data
emptyNet = fitnet([10,10,10],'trainbr');
input = reshape(paramDataSet.Rd,[],25);
% add noise
inputExtended = repmat(input,10,1);
inputNoise = inputExtended.*normrnd(1,0.01,size(inputExtended));

net = train(emptyNet,inputNoise',repmat(paramDataSet.Params(:,2),10,1)');
genFunction(net,'HHbNet');
