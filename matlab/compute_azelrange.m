function [AZ, EL, RANGE] = compute_azelrange(user, sat, R, e2, tol)
% Function compute_azelrange
%
% Takes a user and a satellite position in ECEF and computes the azimuth,
% elevation, and range to that satellite

    LOS_ENU = compute_LOS_ENU(user, sat, R, e2, tol);
    RANGE = sqrt(sum((sat-user).^2 ,1));
    AZ = atan2d(LOS_ENU(1, :), LOS_ENU(2, :));
    EL = asind(LOS_ENU(3, :));  
end


function R = R1(theta)
% Function R1
%
% takes an angle in degrees and generates the R1 rotation matrix for it
    R = [1  0           0
         0  cosd(theta) sind(theta)
         0 -sind(theta) cosd(theta)];
end

function R = R3(theta)
% Function R3
%
% takes an angle in degrees and generates the R3 rotation matrix for it
    R = [cosd(theta) sind(theta) 0
        -sind(theta) cosd(theta) 0
         0           0           1];
end

function [phi_gd, lambda, h] = ECEF2LLA_single(x, y, z, Re, e2, tol)
% Function ECEF2LLA_single
%
% Takes 3x1 vector v which are xyz coordinates in ECEF
%
% Returns geodetic latitude, longitude, and altitude 
    lambda = atan2d(y, x); % calculate longitude
    rho = sqrt(x^2 + y^2); % calculate 2d range
    r = sqrt(x^2 + y^2 + z^2); % calculate range from center
    phi_gd = 0; % initialize geodetic latitude
    phi_gd_old = asind(z / r); % guess geocentric latitude
    C = 0; % initialize C
    while true % iterate over algorithm to converge on geodetic latitude
        C = Re / sqrt(1 - e2 * sind(phi_gd_old)^2); 
        phi_gd_new = atand((z + e2*C*sind(phi_gd_old)) / rho);
        if norm(phi_gd_new - phi_gd_old) < tol
            phi_gd = phi_gd_new;
            break;
        end
        phi_gd_old = phi_gd_new;
    end
    h = rho / cosd(phi_gd) - C; % solve for altitude
end

function lla = ECEF2LLA(v, Re, e2, tol)
% Function ECEF2LLA
%
% Takes 3xM vector v where each column is xyz coordinates in ECEF and calls
% ECEF2LLA_single via arrayfun
%
% Returns geodetic latitude, longitude, and altitude 
    x = v(1, :);
    y = v(2, :);
    z = v(3, :);
    [phi_gd, lambda, h] = ...
        arrayfun(@(x, y, z) ECEF2LLA_single(x, y, z, Re, e2, tol), ...
        x, y, z);
    lla = [phi_gd; lambda; h];
end

function R = ECEF2ENU(lat, lon)
% Function ECEF2ENU
%
% Takes a geodetic latitude and longitude in degrees and returns the
% transformation matrix to rotate from ECEF to ENU
    R = R1(90 - lat)*R3(90+lon);
end

function LOS_ENU = compute_LOS_ENU(user, sat, R, e2, tol)
% Function compute_LOS_ENU
%
% Takes a user and a satellite position in ECEF and computes the LOS unit
% vector to the satellite
%
% Coordinates are expressed in 3xM format (columns are coordinates)
    v = sat - user;
    v = v ./ sqrt(sum(v.^2, 1));
    lla = ECEF2LLA(user, R, e2, tol);
    
    LOS_ENU = ECEF2ENU(lla(1), lla(2)) * v;
    
end