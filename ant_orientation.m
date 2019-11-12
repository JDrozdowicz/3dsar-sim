% Copyright 2019 Jedrzej Drozdowicz
% 
% Licensed under the Apache License, Version 2.0 (the "License");
% you may not use this file except in compliance with the License.
% You may obtain a copy of the License at
% 
%     http://www.apache.org/licenses/LICENSE-2.0
% 
% Unless required by applicable law or agreed to in writing, software
% distributed under the License is distributed on an "AS IS" BASIS,
% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing permissions and
% limitations under the License.

function [az, el] = ant_orientation(ant_pos, p)
%ANT_ORIENTATION Calculate azimuth and elevation angles of point p with
%respect to the antenna

ix = [1,0,0].';
iy = [0,1,0].';
iz = [0,0,1].';

dr = ant_pos(4:6).';
tr = ant_pos(7);

% zaczynam obliczenia

v = cross(iz, dr);
c = dot(iz, dr);

vx = [0, -v(3), v(2); v(3), 0, -v(1); -v(2), v(1), 0];

R = eye(3) + vx + vx^2 .*(1/(1+c));

if  abs(dot(dr,-iz) - 1) < eps  %isequal(dr,-iz)
    ixp = -ix;
    iyp = iy;
else
    ixp = R*ix;
    iyp = R*iy;
end
izp = dr;

ixpp = ixp*cos(tr) + iyp*sin(tr);
iypp = -ixp*sin(tr) + iyp*cos(tr);
izpp = izp;

ppp = p - ant_pos(1:3).';

pxpp = dot(ppp,ixpp);
pypp = dot(ppp,iypp);
pzpp = dot(ppp,izpp);

az = atan2(pypp,pzpp);
el = atan2(pxpp,pzpp);

