function [PRIF, iono] = ionocorr(C1, f1, P2, f2)
    
    PRIF = f1^2/(f1^2 - f2^2) * C1 - f2^2/(f1^2 - f2^2) * P2;
    iono = f2^2/(f1^2 - f2^2) * (P2 - C1);

end