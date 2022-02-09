function RANGE = compute_expected_range(user, ephem_all, prn, t)
%%
    c = 299792458;
    [~, sat] = broadcast_eph2pos(ephem_all, t, prn);
    rg = sqrt(sum((user - sat').^2));
    for i = 1:5
        tt = [t(:, 1), t(:, 2) - rg'/c];
        [~, sat] = broadcast_eph2pos(ephem_all, tt, prn);
        rg = sqrt(sum((user - sat').^2));
        phi = 2*pi/86164.0905 * (rg/c); 
        sat_r = [cos(phi).*sat(:, 1)' + sin(phi).*sat(:, 2)'
            -sin(phi).*sat(:, 1)' + cos(phi).*sat(:, 2)'
            sat(:, 3)'];
        rg = sqrt(sum((user - sat_r).^2));
    end
    RANGE = rg;
end