function [gps_ephem,gps_ephem_cell] = read_GPSyuma(yumafilename,rollovers)

%==========================================================================
%==========================================================================
% [gps_ephem,gps_ephem_cell] = read_GPSyuma(yuma_filename)
%
% Reads in a GPS YUMA almanac and constructs a matrix of all ephemeris 
%  values.
%
%
% Author: Ben K. Bradley
% Date: 02/22/2011
%
%
% INPUT:               Description                                   Units
%
%  yumafilename   - name of YUMA almanac file to read in             string
%  rollovers  - number of 1024 week rollovers   (integer)
%                if omitted it assumes 0.
%
%
% OUTPUT:       
%    
%  gps_ephem      - matrix of gps satellite orbit parameters         (nx25)
%  
%                  col1: prn, PRN number of satellite
%                  col2: M0, mean anomaly at reference time, rad
%                  col3: blank (zero)
%                  col4: ecc, eccentricity of orbit
%                  col5: sqrt_a, square root of semi-major axis, m^0.5
%                  col6: Loa, longitude of ascending node of orbit plane at weekly epoch, rad
%                  col7: incl, inclination angle at reference time, rad
%                  col8: perigee, argument of perigee, rad
%                  col9: ra_rate, rate of change of right ascension, rad/s
%                 col10: blank (zero)
%                 col11: blank (zero)
%                 col12: blank (zero)
%                 col13: blank (zero)
%                 col14: blank (zero)
%                 col15: blank (zero)
%                 col16: blank (zero)
%                 col17: Toe, reference time ephemeris (seconds into GPS week)
%                 col18: blank (zero)
%                 col19: GPS_week, GPS Week Number (to go with Toe)
%                 col20: blank (zero)
%                 col21: Af0, satellite clock bias (sec)
%                 col22: Af1, satellite clock drift (sec/sec)
%                 col23: blank (zero)
%                 col24: blank (zero)
%                 col25: health, satellite health (0=good and usable)
%
%   P. Axelrad - added cell array for people who prefer that.
%   P. Axelrad 9/15/2019 - added an input for the number of week rollovers
%                                     
% Coupling:
%
%   none
%
% References:
% 
%   [1] Interface Control Document: IS-GPS-200D
%         < http://www.navcen.uscg.gov/gps/geninfo/IS-GPS-200D.pdf >
%
%==========================================================================
%==========================================================================

if nargin < 2
    rollovers = 0;
    disp('Assuming rollovers = 0, WN will be less than 1024')
end

% Open desired YUMA almanac ===============================================
if (exist(yumafilename,'file') == 2)

    fid = fopen(yumafilename,'r');
else
    error(sprintf('Unable to find YUMA almanac: %s',yumafilename), 'ERROR!');
end


gps_ephem_cell={};

j = 1;

% Enter main loop to read the rest of the ephemeris file
%==========================================================================
while 1
    
    % Load next line in ephemeris file
    tline = fgetl(fid); %junk line of text eg. ***Week ....***
    
    
    % If the next line is not a character then the end of the file has been
    %   reached and the while loop is exited
    if ~ischar(tline), break, end
    
    
    tline = fgetl(fid); %first line of real satellite info
    
    
    % If the next line is not a character then the end of the file has been
    %   reached and the while loop is exited
    if ~ischar(tline), break, end
   

        eph={};
    
        %-----------------------------------------------------------------
        % Read in variables of the FIRST line this satellite's ephemeris
        %-----------------------------------------------------------------
        eph.prn = str2num(tline(27:end));   %prn number of satellite
        
        %-----------------------------------------------------------------
        % SECOND LINE
        %-----------------------------------------------------------------
        tline  = fgetl(fid); %read in second line of satellite ephemeris
        
        eph.health = str2num(tline(27:end)); % health of satellite
        
        
        %-----------------------------------------------------------------
        % THIRD LINE
        %-----------------------------------------------------------------
        tline = fgetl(fid); %read in third line of satellite ephemeris
        
        eph.ecc   = str2num(tline(27:end)); % eccentricity
        
        
        %-----------------------------------------------------------------
        % FOURTH LINE
        %-----------------------------------------------------------------
        tline = fgetl(fid); %read in fourth line of satellite ephemeris
        
        eph.Toe   = str2num(tline(27:end)); % Reference Time Ephemeris (sec into GPS week)
        
        
        %-----------------------------------------------------------------                    
        % FIFTH LINE
        %-----------------------------------------------------------------
        tline = fgetl(fid); %read in fifth line of satellite ephemeris
        
        eph.incl  = str2num(tline(27:end)); % Inclination Angle at Toe, rad
        
        
        %-----------------------------------------------------------------
        % SIXTH LINE
        %-----------------------------------------------------------------
        tline   = fgetl(fid); %read in sixth line of satellite ephemeris
        
        eph.ra_rate = str2num(tline(27:end)); % Rate of change of Longitude
                                          %  of ascending node, rad/s
                                                    
        %-----------------------------------------------------------------
        % SEVENTH LINE
        %-----------------------------------------------------------------
        tline   = fgetl(fid); %read in sixth line of satellite ephemeris
        
        eph.sqrt_a  = str2num(tline(27:end)); % square root of Semi-Major axis (m^1/2)                                           
                                                    
                                                    
        
        %-----------------------------------------------------------------
        % SEVENTH LINE 
        %-----------------------------------------------------------------
        tline = fgetl(fid); %read in seventh line of satellite ephemeris
        
        eph.Loa   = str2num(tline(27:end)); % Longitude of ascending node at beginning
                                        %  of GPS week, rad
        
        %-----------------------------------------------------------------
        % EIGHTH LINE 
        %-----------------------------------------------------------------
        tline   = fgetl(fid); %read in eighth line of satellite ephemeris
        
        eph.perigee = str2num(tline(27:end)); % Argument of Perigee, rad
        
        
        %-----------------------------------------------------------------
        % NINTH LINE 
        %-----------------------------------------------------------------
        tline = fgetl(fid); %read in ninth line of satellite ephemeris
        
        eph.M0    = str2num(tline(27:end)); % Mean Anomaly, rad
        
        
        %-----------------------------------------------------------------
        % TENTH LINE 
        %-----------------------------------------------------------------
        tline = fgetl(fid); %read in tenth line of satellite ephemeris
        
        eph.Af0   = str2num(tline(27:end)); % clock bias, seconds
        
        
        %-----------------------------------------------------------------
        % ELEVENTH LINE 
        %-----------------------------------------------------------------
        tline = fgetl(fid); %read in eleventh line of satellite ephemeris
        
        eph.Af1   = str2num(tline(27:end)); % clock drift, s/s
        
        
        %-----------------------------------------------------------------
        % TWELFTH LINE 
        %-----------------------------------------------------------------
        tline = fgetl(fid); %read in twelfth line of satellite ephemeris
        
        GPS_week = str2num(tline(27:end)); % GPS week number
%         if GPS_week < 1024
%             GPS_week = GPS_week+1024; % Penny added to avoid modulo 1024
%         end
        eph.GPS_week = GPS_week+rollovers*1024;
        
        
        % Variables not included in the YUMA almanac (but are in broadcast ephemeris)
        eph.delta_n = 0;    eph.i_rate = 0;
        eph.Cuc     = 0;    eph.Cus    = 0;
        eph.Crc     = 0;    eph.Crs    = 0;
        eph.Cic     = 0;    eph.Cis    = 0;
        eph.IODE    = 0;    eph.Toc    = 0; %Time of clock
        eph.Af2     = 0;
        
        
        
        gps_ephem(j,:) = [eph.prn eph.M0 eph.delta_n eph.ecc eph.sqrt_a ...
            eph.Loa eph.incl eph.perigee eph.ra_rate eph.i_rate eph.Cuc ...
            eph.Cus eph.Crc eph.Crs eph.Cic eph.Cis eph.Toe eph.IODE ...
            eph.GPS_week eph.Toc eph.Af0 eph.Af1 eph.Af2 0 eph.health];
        
        gps_ephem_cell{j}=eph;
        
        
        j = j + 1; 
        
        tline = fgetl(fid);%get blank line inbetween satellite info
        
        % If the next line is not a character then the end of the file has been
        %   reached and the while loop is exited
        if ~ischar(tline), break, end
    
    
end

fclose(fid);





