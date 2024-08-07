function [Rd, Rdlog, pix] = mcxdemodRd(fluence, cfg, icx)
    % Setting up to grab the three phases
    idx1 = icx*3 -2; % Phase 1
    idx2 = icx*3 -1; % Phase 2
    idx3 = icx*3;    % Phase 3

    X = size(fluence(idx1).dref,1); % X dimension in voxels
    Y = size(fluence(idx1).dref,2); % Y dimension in voxels

    scale = cfg.unitinmm;
    illumscale = (X - 2*cfg.margin/scale) * (Y - 2*cfg.margin/scale); % Determining the size of the illuminated area in mm^2
    %illumscale = (X*scale - 2*cfg.margin) * (Y*scale - 2*cfg.margin); % Determining the size of the illuminated area in mm^2
    areascale = 1/illumscale;
    scalefact = areascale*cfg.nphoton; % photon density factor (in photons/mm^2)
    
    p1 = fluence(idx1).dref(:,:,1)./scalefact; % Adjusting phase-pixel intensities
    p2 = fluence(idx2).dref(:,:,1)./scalefact; % based on overall photon density
    p3 = fluence(idx3).dref(:,:,1)./scalefact; % of light projection
    dm = zeros(X,Y);
    weightsum = 0;
    normavg = (fluence(idx1).stat.normalizer + fluence(idx2).stat.normalizer + fluence(idx3).stat.normalizer)/3;

    %%
    if icx == 1
        for i = 1+cfg.margin/scale:(X-cfg.margin/scale)
            for j = 1+cfg.margin/scale:(Y-cfg.margin/scale)
                demval = (sqrt(2)/3) *    ...
                    (((p1(i,j)-p2(i,j))^2 ...
                    + (p2(i,j)-p3(i,j))^2 ...
                    + (p3(i,j)-p1(i,j))^2) ^ 0.5);
                % demval = p1(i,j);
                scaleval = demval * scalefact;
                dm(i,j) = scaleval;
                weightsum = weightsum + scaleval;
            end
        end
    else
        for i = 1+cfg.margin/scale:(X-cfg.margin/scale)
            for j = 1+cfg.margin/scale:(Y-cfg.margin/scale)
                demval =(sqrt(2)/3) *     ...
                    (((p1(i,j)-p2(i,j))^2 ...
                    + (p2(i,j)-p3(i,j))^2 ...
                    + (p3(i,j)-p1(i,j))^2) ^ 0.5);
                scaleval = demval * scalefact;
                dm(i,j) = scaleval;
                weightsum = weightsum + scaleval;
            end
        end
    end
    %%
    zedm = (size(dm,1)*size(dm,2)) - nnz(dm);
    zrat = zedm / (size(dm,1)*size(dm,2));
    
    if icx == 1
        energytot = (fluence(idx1).stat.energytot + fluence(idx2).stat.energytot + fluence(idx3).stat.energytot)/3;
        energynorm = (fluence(idx1).stat.normalizer + fluence(idx2).stat.normalizer + fluence(idx3).stat.normalizer)/3;
    else
        energytot = (fluence(idx1).stat.energytot + fluence(idx2).stat.energytot + fluence(idx3).stat.energytot)/3;
        energynorm = (fluence(idx1).stat.normalizer + fluence(idx2).stat.normalizer + fluence(idx3).stat.normalizer)/3;
    end
    
    
    %%
    Rd = weightsum * 1 / energytot;

    Rdlog = log10(Rd);
    
    Rdtest = weightsum / energytot / energynorm;
    
    Rdtestlog = log10(Rdtest);
    
    pix = dm;
end