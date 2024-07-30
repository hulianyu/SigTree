function p = simplifiedChi2cdf(x,v)
% Simplified GAMCDF Gamma cumulative distribution function.
z = x/2;
a = v/2;
p = gammainc(z, a,'upper');
end