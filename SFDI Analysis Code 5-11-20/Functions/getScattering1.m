function musp = getScattering1(lambdas,a,bMie)
    musp = a*((lambdas/500).^(-bMie));
end