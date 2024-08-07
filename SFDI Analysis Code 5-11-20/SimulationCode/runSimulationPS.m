%% RUN_SIMULATION_PS
% Run a single point-source Monte Carlo Simulation on the GPU. Accepts 
% optical properties as arguments and returns photon data as the output.
%
%--- Inputs ---
% MUA: list of absorption coefficients (mm^-1)
% MUSP: list of scattering coefficients (mm^-1)
% G: list of scattering anisotropy values for henyey-Greenstein phase 
%   funciton (unitless)
% N: list of real refractive indices (unitless)
% Z: list of layer boundary depths for multi-layer simulations (mm)
% NPHOT: number of photons to simulate (photons)
%
%--- Outputs ---
% R: data structure containing photon data from simulation
%   R.PPATH: nLayer-by-nPhot array of cumulative partial pathlengths for 
%       each layer (mm)
%   R.MOM: nLayer-by-nPhot array of cumulative partial momentum transfer
%       (1-cos(scattering angle)) for each layer (unitless)
%   R.NSCAT: nLayer-by-nPhot array of total scattering events for each
%       layer (scattering events)
%   R.FINALX: nPhot array of final x positions of each photon (mm)
%   R.FINALY: nPhot array of final y positions of each photon (mm)
%   R.FINALZ: nPhot array of final z positions of each photon (mm)
%   R.FINALUX: nPhot array of final x directional cosine for each photon
%   R.FINALUY: nPhot array of final y directional cosine for each photon
%   R.FINALUZ: nPhot array of final z directional cosine for each photon
%   R.FINALW: nPhot array of final photon weight for each photon
%   R.EXIT: nPhot array of photon exit flags
%       (1: photon reflected/transmitted, 0: photon absorbed)
%
%--- EXAMPLES ---
% r = runSimulationPS(0.02,1.5,0.9,1.4,10,10^7);
% This runs a homogeneous medium simulation with air as the surroundings. The
% z parameter represents a medium with a total thickness of 10 mm. 10 million
% photons are simulated. Note that with thinner media, many photons will
% transmit, so extra care is needed to filter out photons that will not be
% detected via diffuse reflection.
%
% r = runSimulationPS([0.01,0.02],[1,1.5],[0.9,0.9],[1.4,1.4],[2,1000],10^7);
% This runs a two-layer medium simulation with air as the surroundings. The
% z parameters represents a medium with a top layer of 2 mm and a bottom layer
% thickness of 1000 - 2 = 998 mm. The total thickness of the slab is 1000
% mm. 10 million photons are simulated.
%
function r = runSimulationPS(mua,musp,g,n,z,nPhot)
% Generate a random integer seed for the simualtion. This helps make each
% consecutive simulation start with different initial conditions
randSeed = randi([1,10000]);

% populate optical properties for each layer. By default, the medium is
% assumed to be in air (mua,musp=0 and g,n=1) and the top boundary is z=0mm
mua = [0 mua];
musp = [0 musp];
g = [1 g];
n = [1 n];
layerZ = [0 z];

% Initialize photon position and directions. For a point-source it is
% convenient to have the initial position to be at the origin. By defualt,
% the initial directional cosines are [ux=0,uy=0,uz=1], meaning that each
% photon is initially propagating purely in z direction
x = 0;
y = 0;
z = 0;
ux = 0;
uy = 0;
uz = 1;

% Prepare the GPU for the simulation by pre-compiling the kernel
k = initGPUKernel(nPhot,'MC_PS');

% Define the maximum number of scattering events allowed. This helps the
% simulation finish instead of waiting for a couple photons that haven't
% been absorbed/reflected/transmitted (useful when absorption is very low)
scatMax = 500000;
scatMax = gpuArray(int32(scatMax));

% Get the number of layers in the medium
nLayers = length(mua);
% Send photon and medium data structures to GPU global memory
% Check to make sure every optical property is defined for each layer
pData = gpuArray(single([x,y,z,ux,uy,uz]));
try
    mData = [mua;musp./(1-g);g;n;layerZ];
catch
    error('Optical properties must be defined for each layer.');
end
mData = gpuArray(single(mData(:)'));

% Allocate memory for returned photon data on the GPU. GPUARRAY simply
% sends data to the GPU to do calculations and returns the data when
% finished
finalx = gpuArray(single(zeros(1,nPhot)));
finaly = gpuArray(single(zeros(1,nPhot)));
finalz = gpuArray(single(zeros(1,nPhot)));
finalux = gpuArray(single(zeros(1,nPhot)));
finaluy = gpuArray(single(zeros(1,nPhot)));
finaluz = gpuArray(single(zeros(1,nPhot)));
finalw = gpuArray(single(zeros(1,nPhot)));
maxZ = gpuArray(single(zeros(1,nPhot)));
exit = gpuArray(int32(zeros(1,nPhot)));
ppath = gpuArray(single(zeros(1,nPhot*(nLayers-1))));
mom = gpuArray(single(zeros(1,nPhot*(nLayers-1))));
nscat = gpuArray(int32(zeros(1,nPhot*(nLayers-1))));
nPhot = gpuArray(int32(nPhot));
randSeed = gpuArray(int32(randSeed));
nLayers = gpuArray(int32(nLayers));

% This function calls the GPU kernel k with all the GPU data previosuly
% defined and retuns the same arrays popualted by the GPU
[ppath,mom,nscat,finalx,finaly,finalz,finalux,finaluy,finaluz,finalw,maxZ,exit] = feval(k,nPhot,...
    nLayers,randSeed,scatMax,pData,mData,ppath,mom,nscat,...
    finalx,finaly,finalz,finalux,finaluy,finaluz,finalw,maxZ,exit);

% The GATHER function copies GPU data back to the CPU for MATLAB to use.
% All the data generated by the GPU is copied back to the CPU and loaded
% into a single data structure r
ppath = gather(ppath);
mom = gather(mom);
nscat = gather(nscat);
r.ppath = reshape(ppath,(nLayers-1),[]);
r.mom = reshape(mom,(nLayers-1),[]);
r.nscat = reshape(nscat,(nLayers-1),[]);
r.finalx = gather(finalx);
r.finaly = gather(finaly);
r.finalz = gather(finalz);
r.finalux = gather(finalux);
r.finaluy = gather(finaluy);
r.finaluz = gather(finaluz);
r.finalw = gather(finalw);
r.maxZ = gather(maxZ);
r.exit = gather(exit);
end