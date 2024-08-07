function [coeff,fit] = spectrumFit2(spectrum1,spectrum2,mua_measured)
    fun = @(params)specFit(params,spectrum1,spectrum2,mua_measured);
    options = optimoptions('fmincon','display','none');
    params0 = [0,0];
    [coeff,fit] = fmincon(fun,params0,[],[],[],[],[0,0],[],[],options);
    function error = specFit(params,spec1,spec2,mua_measured)
        error = mean((mua_measured - (params(1)*spec1 + params(2)*spec2)).^2);
    end
end