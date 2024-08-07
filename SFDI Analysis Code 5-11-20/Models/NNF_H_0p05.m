function [y1] = NNF_H_0p05(x1)
%NNF_H_0P05 neural network simulation function.
%
% Generated by Neural Network Toolbox function genFunction, 06-Dec-2019 13:14:45.
% 
% [y1] = NNF_H_0p05(x1) takes these arguments:
%   x = 2xQ matrix, input #1
% and returns:
%   y = 1xQ matrix, output #1
% where Q is the number of samples.

%#ok<*RPMT0>

% ===== NEURAL NETWORK CONSTANTS =====

% Input 1
x1_step1.xoffset = [0.00113107484963635;0.237537409132807];
x1_step1.gain = [10.1668945937144;0.420166092857931];
x1_step1.ymin = -1;

% Layer 1
b1 = [5.9729407689105702417;-3.9765995434281120602;6.9047658680581793789;2.3729344666635787675;-4.7560127425295037895;2.3752707985680729941;2.390445466119087925;1.5444429992994004319;0.45323386414070404582;-0.35892717656719952402;1.0518300667352562527;-0.89173294673830305612;0.19262981495301206847;1.553874100374373235;-3.2485536717498142423;-3.4547395113254846244;3.4817789990299714553;3.1679763722028160267;4.1316814517702509235;-8.4479190358432436625];
IW1_1 = [-2.9787770157065103227 -3.4587019369929508095;3.4736405970500228158 2.2778071339327730094;1.0894710429613396752 5.9824401478395614618;0.50298917024024702993 1.6479176120304248609;-0.0094800716840307396732 -5.912755513775520555;-0.014442357811778380045 3.1768709626281674829;-0.041794122896093778563 3.0100682537809522188;-2.1500924632723315 -1.0370771448597870812;-1.9670983908575885302 -0.4802024293823992096;-0.52492020691212126771 4.6368812627150060024;1.1292528461930948502 0.0095593254230726034532;0.39217076827005459272 -1.6238093321083015841;0.49930272241718348658 -0.86133193552083897604;1.2049581800387307862 -1.2279615901900828856;-0.42471886508312167718 4.8123416911245451288;-4.2136117912013144604 -4.0865550428295396301;2.9861864918568050697 -0.79914272050080426002;2.5931661602350004614 0.36557065906350422191;1.7295539157248609463 1.439672106206147717;-6.393615082458125265 0.35388130326601957565];

% Layer 2
b2 = -0.5426148686257857845;
LW2_1 = [0.0095001603221814903688 -0.0046043147429472193466 0.075285291851387797779 0.47845408968855973608 0.03968935686678673469 -0.70596455434396510054 0.99808680905549562734 0.010364781456906845686 0.011555529186120299875 0.0020726721032850808045 -0.41017855114382295589 -0.15686176597967985136 -0.26429719927160583159 -0.051968550414735033283 -0.0014228290122209562556 -0.00016676775012408758501 -0.082333238719764376135 -0.40990414752416326483 2.0072232331596819854 1.5655766605625351939];

% Output 1
y1_step1.ymin = -1;
y1_step1.gain = 2.5004659320444;
y1_step1.xoffset = 0.0348720662295818;

% ===== SIMULATION ========

% Dimensions
Q = size(x1,2); % samples

% Input 1
xp1 = mapminmax_apply(x1,x1_step1);

% Layer 1
a1 = tansig_apply(repmat(b1,1,Q) + IW1_1*xp1);

% Layer 2
a2 = repmat(b2,1,Q) + LW2_1*a1;

% Output 1
y1 = mapminmax_reverse(a2,y1_step1);
end

% ===== MODULE FUNCTIONS ========

% Map Minimum and Maximum Input Processing Function
function y = mapminmax_apply(x,settings)
  y = bsxfun(@minus,x,settings.xoffset);
  y = bsxfun(@times,y,settings.gain);
  y = bsxfun(@plus,y,settings.ymin);
end

% Sigmoid Symmetric Transfer Function
function a = tansig_apply(n,~)
  a = 2 ./ (1 + exp(-2*n)) - 1;
end

% Map Minimum and Maximum Output Reverse-Processing Function
function x = mapminmax_reverse(y,settings)
  x = bsxfun(@minus,y,settings.ymin);
  x = bsxfun(@rdivide,x,settings.gain);
  x = bsxfun(@plus,x,settings.xoffset);
end
