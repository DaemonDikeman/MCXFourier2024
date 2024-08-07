function [a,bMie,fit] = scatFit1(lambdas,musp_measured)

fun = @(params)scatFit(params,lambdas,musp_measured);
options = optimoptions('fmincon','display','none');
params0 = [0,0];
[soln,fit] = fmincon(fun,params0,[],[],[],[],[0,0],[10,10],[],options);
a = soln(1);
bMie = soln(2);

function error = scatFit(params,lambdas,musp_measured)
    af = params(1);
    bMief = params(2);
    musp = af*((lambdas/500).^(-bMief));
    error = mean((musp_measured - musp).^2);
end

end