function S = resample_digital(t, X, fx, phi)
    nc = length(X);
    delay = phi / nc;
    T = nc / fx;
    t = t + T*(delay+1);
    S = X(mod(floor(t*fx), nc) + 1);
end