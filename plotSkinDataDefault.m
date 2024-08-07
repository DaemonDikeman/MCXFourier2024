% Create Figure
f=figure;
subplot(2,3,1);
%image(out.color);
%axis image; 
axis off;
title(out.name);
cidx=[1:5];
% Loop through Chromophores
for c=1:5
    subplot(2,3,c+1);
    imagesc(squeeze(out.MaskMaps).*squeeze(out.ChromMaps(:,:,cidx(c))),[0 out.limits(cidx(c))])
    %axis image; 
    axis off; colorbar;    
    title(out.ctitles{cidx(c)});
end