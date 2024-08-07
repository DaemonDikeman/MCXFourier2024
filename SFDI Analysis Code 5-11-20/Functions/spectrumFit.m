function [coeff,fit] = spectrumFit(spectrum,mua_measured)
    fun = @(params)specFit(params,spectrum,mua_measured);
    options = optimoptions('fmincon','display','none');
    params0 = 0;
    [coeff,fit] = fmincon(fun,params0,[],[],[],[],0,[],[],options);
    function error = specFit(params,spectrum,mua_measured)
        error = mean((mua_measured - params(1)*spectrum).^2);
    end
end