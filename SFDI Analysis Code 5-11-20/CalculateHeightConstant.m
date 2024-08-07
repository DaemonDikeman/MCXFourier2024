%% Import Height Calibration Data Set
[dataSet,phantomID] = importDataSet([0,0.05,0.1,0.15,0.2],[0,1000],'MI');
dataSet = demodDataSet(dataSet);
%% Reference Phase Map
fxID = 1;
I1 = dataSet{phantomID}.Profile(:,:,1,fxID);
I2 = dataSet{phantomID}.Profile(:,:,2,fxID);
I3 = dataSet{phantomID}.Profile(:,:,3,fxID);
%{
smooth = 1;
I1 = imgaussfilt(I1,smooth);
I2 = imgaussfilt(I2,smooth);
I3 = imgaussfilt(I3,smooth);
%}
wrapppedPhase = atan2((sqrt(3)*(I1-I3)),(2*I2 - I1 - I3));
figure
imagesc(wrapppedPhase)
refPhase = unwrap(wrapppedPhase')';
%% Sample phase map 1
% -1cm below reference phantom
ID = 3;
disp(dataSet{ID}.Name)
fxID = 1;
I1 = dataSet{ID}.Profile(:,:,1,fxID);
I2 = dataSet{ID}.Profile(:,:,2,fxID);
I3 = dataSet{ID}.Profile(:,:,3,fxID);
%{
smooth = 1;
I1 = imgaussfilt(I1,smooth);
I2 = imgaussfilt(I2,smooth);
I3 = imgaussfilt(I3,smooth);
%}
wrapppedPhase = atan2((sqrt(3)*(I1-I3)),(2*I2 - I1 - I3));
phase1 = unwrap(wrapppedPhase')';
x1 = mean2(phase1(100:300,100:300)-refPhase(100:300,100:300));
%% Sample phase map 2
% +1cm
ID = 5;
disp(dataSet{ID}.Name)
fxID = 1;
I1 = dataSet{ID}.Profile(:,:,1,fxID);
I2 = dataSet{ID}.Profile(:,:,2,fxID);
I3 = dataSet{ID}.Profile(:,:,3,fxID);
%{
smooth = 1;
I1 = imgaussfilt(I1,smooth);
I2 = imgaussfilt(I2,smooth);
I3 = imgaussfilt(I3,smooth);
%}
wrapppedPhase = atan2((sqrt(3)*(I1-I3)),(2*I2 - I1 - I3));
phase2 = unwrap(wrapppedPhase')';
x2 = mean2(phase2(100:300,100:300)-refPhase(100:300,100:300));
%% Calculate Proportionality Constant
actual_height = [-1,1];
proportional_height = [x1,x2];
s = (actual_height(2) - actual_height(1))/(proportional_height(2)-proportional_height(1));
disp(['Height Factor = ', num2str(s)]);