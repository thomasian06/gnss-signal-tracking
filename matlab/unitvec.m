function [r,u] = unitvec(x,dim)
% function [r,u] = unitvec(x,dim)
% Input is array, x
% dim is the dimension that has 3 elements to be normalized.
% i.e. if it is an nx3 you would have (x,2)
% Outputs are length of the vector, and the unit vector

r = sqrt(sum(x.^2,dim));
if dim == 1
    u = x./repmat(r,3,1);
else 
    u = x./repmat(r,1,3);
end
