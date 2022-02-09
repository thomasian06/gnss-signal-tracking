function [wn,tow,wn2]= cal2gps2021(date)
%==========================================================================
% [wn,tow,wn2]= cal2gps2021(date)
%
% Converts a MATLAB datetime array into GPS week number and TOW 
%
%
% Author: P. Axelrad
% Date: September 9, 2021
%
%
% INPUT:               Description                                   Units
%
%  date   - MATLAB datetime array assumed to be GPS time represented as a
%  calendar date
%
%
% OUTPUT:       
%    
%  wn - GPS week number since January 6, 1980.
%  tow - time of week in seconds
%  wn2 - GPS week number modulo 1024 to be used in almanac
%==========================================================================


% GPS start date January 6, 1980
gpsrefdate = datetime(1980,01,06);  
dt = between(gpsrefdate,date,{'Weeks' 'Time'});

[wn,tow]=split(dt,{'Weeks','Time'});
tow = seconds(tow);
wn2 = rem(wn,1024);

end
