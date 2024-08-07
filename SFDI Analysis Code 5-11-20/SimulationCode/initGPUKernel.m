%% INIT_GPU_KERNEL
% Prepares CUDA GPU code to be run from MATLAB code.
%
%--- Inputs ---
% NPHOT: number of photons to be simulated (number of photons)
% PTXFILENAME: name of pre-compiled .ptx file (string w/o file extention)
%
%--- Outputs ---
% K: MATLAB CUDA kernel handle for use in FEVAL function
%
%--- EXAMPLES ---
% k = initGPUKernel(10^7,'MC_PS');
% This prepares the point-source kernel to simulate 10 million photons
%
% k = initGPUKernel(10^7,'MC_SFD');
% This prepares the spatial frequency domain kernel to simulate 10 million photons
function k = initGPUKernel(nPhot,ptxFileName)
% Load kernel from raw .cu and .ptx files. PTX (parallel thread execution)
% files are pre-compiled CUDA files to be run on the GPU. They contain
% everything but device memory allocation and data return routines. This is
% handeled in MATLAB'S GPU library
k = parallel.gpu.CUDAKernel([ptxFileName '.ptx'],[ptxFileName '.cu']);

% Define GPU execution architecture based on the number of photons
k.ThreadBlockSize = [k.MaxThreadsPerBlock,1,1];
k.GridSize = [ceil(nPhot/k.MaxThreadsPerBlock),1];
end