function [dref] = mcxcwbbdref(detp, cfg, fluence)
%
%  [dref] = mcxcwbbdref(detp, cfg)
%
%  Compute CW diffuse reflectance from MC detected photon profiles
%  while using (only) bounding box detectors
% 
%  author: Shijie Yan (yan.shiji <at> northeastern.edu)
%  edited by: Daemon Dikeman (daemon.dikeman <at> maine.edu)
%
%  input:
%    detp: profiles of detected photons
%    cfg:  a struct, or struct array. Each element of cfg defines 
%         the parameters associated with a MC simulation.
%
%  output:
%    dref: CW diffuse reflectance at detectors
%
%    this file is part of Monte Carlo eXtreme (MCX)
%    License: GPLv3, see http://mcx.sf.net for details
%    see Yao2018
%
    detweight = mcxdetweight(detp, cfg.prop, cfg.unitinmm);
    detnum = length(unique(detp.detid));
    detweightsum = zeros(detnum, 1);
    for i = 1 : length(detp.detid)
        detweightsum = detweightsum + detweight(i);
    end
    boxbounds = cfg.bc(7:12);
    bnum = double(uint16(bin2dec(boxbounds)));
    bnum2 = bnum;
    
    [d1, d2, d3] = size(cfg.vol);
    d1 = d1*cfg.unitinmm;
    d2 = d2*cfg.unitinmm;
    d3 = d3*cfg.unitinmm;
    
    area = 0;
    if (bnum2-32) >= 0
        bnum2 = bnum2 - 32;
        area = area + d2*d3;
    end
    if (bnum2-16)>= 0
        bnum2 = bnum2 - 16;
        area = area + d1*d3;
    end
    if (bnum2-8) >= 0
        bnum2 = bnum2 - 8;
        area = area + d1*d2;
    end
    if (bnum2-4) >= 0
        bnum2 = bnum2 - 4;
        area = area + d2*d3;
    end
    if (bnum2-2) >= 0
        bnum2 = bnum2 - 2;
        area = area + d1*d3;
    end
    if (bnum2-1) >= 0
        area = area + d1*d2;
    end

    dref = detweightsum ./ (fluence.stat.energytot * area); % 
end