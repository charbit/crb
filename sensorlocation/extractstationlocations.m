function xsensors_m = extractstationlocations()
data=importdata('Export (3).xls');
[nbraw,nbcol] = size(data);
stations = data(2:nbraw,1);
for numselect=26%38:59
    zone = 31;
    hemisphere = 'N';
    xsensors_m.name        = cell(nbraw-1,1);
    xsensors_m.coordinates = zeros(nbraw-1,3);
    cp=0;
    for ir = 1:nbraw-1
        if str2double(stations{ir}(2:3))==numselect
            cp=cp+1;
            lat  = str2double(data{ir+1,4});
            lon  = str2double(data{ir+1,5});
            elev = str2double(data{ir+1,6});
            [x,y] = ll2utm (lat,lon, zone, hemisphere);
            xsensors_m.coordinates(cp,1:2) = [x,y];
            xsensors_m.coordinates(cp,3) = elev*1000;
            xsensors_m.name{cp} = stations{ir};
        end
    end
    if not(cp==0)
        xsensors_m.coordinates = xsensors_m.coordinates(1:cp,:);
        xsensors_m.coordinates = xsensors_m.coordinates-ones(cp,1)*xsensors_m.coordinates(1,:);
        xsensors_m.name = xsensors_m.name(1:cp);
        eval(sprintf('save %s xsensors_m',xsensors_m.name{1}(1:3)))
    end
end


%============================================================
function [x,y] = ll2utm (lat,lon, zone, hemisphere)
%============================================================
% LL2UTM
%        Conversion (latitude,longitude)
%          en coordonnees UTM planaires
% SYNOPSIS
%	[x,y] = ll2utm(lat, lon, zone, hemisphere)
%	lat : latitude en degr?
%   lon : longitude en degr?
%   zone: numero de zone (France en zone 30,31 et 32)
%   (voir par exemple http://www.dmap.co.uk/utmworld.htm)
%   hemisphere: 'N' pour Nord et 'S' pour Sud
%	x,y : coordonnees planaires en UTM
%   UTM : Universal Transverse Mercator
%=============================================================
RADIUS = 6378137;            % rayon de la terre a l'equateur
FLAT = 1 / 298.257223563;    % applatissement WGS-84
M0 = 0;                      % en UTM la latitude a l'origine est toujours nulle
K_0 = 0.9996;                % Facteur d'echelle du meridien central (Central Meridian)
LARGEUR_ZONE=6;              % intervalle de 6 degres E/O, chacun de 3 degres d'E/O
FE = 500000;                 % faux Est
FN = 10000000;               % faux Nord
DEG2RADS = pi/180;
%=============================================================
% qqs verifications
if (max(abs(lat)) > 90)
    error('la latitude ne doit pas exceder 90 degre');
    return;
end
if ((zone < 1) | (zone > 60))
    error ('les numeros de zones utm vont uniquement entre 1 a 60');
    return;
end
%=====================================
% qq valeurs intermediaires
e2  = 2*FLAT - FLAT*FLAT;
e4  = e2 * e2;
e6  = e4 * e2;
ep2 = e2/(1-e2);
%=====================================
lat = lat * DEG2RADS;
lon = lon * DEG2RADS;
sinLat = sin(lat);
tanLat = tan(lat);
cosLat = cos(lat);
sin2Lat = sin(2*lat);
sin4Lat = sin(4*lat);
sin6Lat = sin(6*lat);
T = tanLat.*tanLat;
C = ep2 * (cosLat.*cosLat);
% les As
lambda_0  = (-180 + zone*LARGEUR_ZONE - LARGEUR_ZONE/2).*DEG2RADS;
A = (lon - lambda_0).*cosLat;
A2 = A.*A;
A3 = A2.*A;
A4 = A3.*A;
A5 = A4.*A;
A6 = A5.*A;
sinLat2 = sinLat.*sinLat;
%=====================================
% rayon au point de latitude considere
radius_lat = RADIUS ./ (sqrt (1-e2*sinLat2));
% les Ms
CM1 = 1 - e2/4 - 3*e4/64 - 5*e6/256;
CM2 = 3*e2/8 + 3*e4/32 + 45*e6/1024;
CM3 = 15*e4/256 + 45*e6/1024;
CM4 = 35*e6/3072;
M = RADIUS.*(CM1.*lat - CM2.*sin2Lat + CM3.*sin4Lat - CM4.*sin6Lat);
% x
X1 = A3/6;
CX1 = 1 - T + C;
X2 = A5/120;
CX2 = 5 - 18*T + T.*T + 72*C - 58*ep2;
x = FE + K_0.*radius_lat.*(A + CX1.*X1 + CX2.*X2);
% y
Y1=A2/2;
CY2 = 5 - T + 9*C + 4*C.*C;
Y2 = A4/24;
CY3 = 61 - 58*T + T.*T + 600*C - 330*ep2;
Y3 = A6/720;
y = (hemisphere=='S')*FN + ...
    K_0.*(M - M0 + radius_lat.*tanLat.*(Y1 + CY2.*Y2 + CY3.*Y3));
%======================== END ======================================