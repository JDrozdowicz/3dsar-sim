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


function [fd] = doppler(fc, pos, tx_pos, rx_pos, v)
%DOPPLER Calculate the Doppler shift.
%   Calculates the doppler shift for given parameters. Positions and
%   velocity must be supplied as vectors.

% Constants
% speed of light [m/s]
c = 3e8;

% For clarity: v = dr/dt, which means the target moving AWAY FROM the
% antenna has positive velocity, and target moving TOWARDS the antenna has
% negative velocity

% Calculate the relative positions:
pos_rel_tx = pos(1:3) - tx_pos(1:3);
pos_rel_rx = pos(1:3) - rx_pos(1:3);

% Normalize
pos_rel_tx_norm = pos_rel_tx/norm(pos_rel_tx);
pos_rel_rx_norm = pos_rel_rx/norm(pos_rel_rx);

% Relative velocity
v_rel_tx = dot(v,pos_rel_tx_norm);
v_rel_rx = dot(v,pos_rel_rx_norm);

% Doppler frequency
fd = -(v_rel_tx+v_rel_rx)/c*fc;

end

