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


function [raw_data] = sim_resp(fc, B, T, cell_size, start_range, end_range, points, tx_pos, rx_pos, fd)
%SIM_RESP Calculate radar response.
%   Calculates FMCW radar response for given parameters, scene and antenna
%   pairs positions. For monostatic radar use tx_pos = rx_pos.
%   size(tx_pos) must be equal to size(rx_pos)

% Check the Doppler Frequency
if ~exist('fd', 'var')
    fd = 0;
end

% Modulation speed
alpha = B/T;

% Constants
% speed of light [m/s]
c = 3e8;

% range axis
r = start_range:cell_size:end_range;
% time axis
t = 2*r/c;

% bistatic range for each point
r_d = (sqrt(sum((tx_pos(1:3) - points(:,1:3)).^2,2)).' + sqrt(sum((rx_pos(1:3) - points(:,1:3)).^2,2)).')/2;
% bistatic time for each point
t_d = 2*r_d/c;

% raw data generation
raw_data = sum(points(:,4).*(1-abs((t-t_d.')/T)).*sinc((fd + alpha.*(t-t_d.')).*T.*(1-abs((t-t_d.')/T))).*exp(-1j.*pi.*fd.*(t-t_d.')).*exp(1j.*2.*pi.*(t-t_d.').*fc).*exp(1j.*points(:,5)),1);

end

