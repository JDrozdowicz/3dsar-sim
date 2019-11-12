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


function [raw_data] = sim_resp(fc, B, cell_size, start_range, end_range, points, tx_pos, rx_pos, tx_pattern, rx_pattern)
%SIM_RESP Calculate radar response.
%   Calculates FMCW radar response for given parameters, scene and antenna
%   pairs positions. For monostatic radar use tx_pos = rx_pos.
%   size(tx_pos) must be equal to size(rx_pos)

% Constants
% speed of light [m/s]
c = 3e8;

% range axis
r = start_range:cell_size:end_range;

% bistatic range for each point
r_d = (sqrt(sum((tx_pos(1:3) - points(:,1:3)).^2,2)).' + sqrt(sum((rx_pos(1:3) - points(:,1:3)).^2,2)).')/2;

ant_gain = ones(size(points,1),1);

% Include TX antenna pattern
for iter_ant = 1:size(points,1)
    [az, el] = ant_orientation(tx_pos,points(iter_ant,1:3).');
    ant_gain(iter_ant) = ant_gain(iter_ant).*tx_pattern(az,el);
end

% Include RX antenna pattern
for iter_ant = 1:size(points,1)
    [az, el] = ant_orientation(tx_pos,points(iter_ant,1:3).');
    ant_gain(iter_ant) = ant_gain(iter_ant).*rx_pattern(az,el);
end

% raw data generation
raw_data = sum(points(:,4).*ant_gain.*sinc((r-r_d.')./(c/(2*B))).*exp(1j.*2.*(r-r_d.')./(c/fc).*2.*pi+points(:,5)),1);


end

