%% PDF Estimation
function [x,P] = pdfEstimate(X,n_bins,XBounds,smoothSize)
%n_bins = floor(sqrt(length(X)));
if isempty(X)
    P = 0;
    x = 0;
else
%n_bins = ceil(sqrt(length(X)));
bin_width = (XBounds(2) - XBounds(1))/n_bins;
P = zeros(n_bins);
x = linspace(XBounds(1),XBounds(2),n_bins);
for i = 1:n_bins
    lower = XBounds(1) + (i-1)*bin_width;
    upper = XBounds(1) + (i)*bin_width;
    idx = X >= lower & X < upper;
    P(i) = sum(idx);
end
% Smooth and normalize
P = smoothdata(P, 'movmean', smoothSize);
P = P/trapz(x,P);
P(isnan(P)) = 0;
P = P';
end
end