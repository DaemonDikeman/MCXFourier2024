function musp = getScattering2(lambdas,a,bMie,fRay)
    musp = a*(fRay*(lambdas/500).^(-4) + (1-fRay)*(lambdas/500).^(-bMie));
end