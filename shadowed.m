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

function [S] = shadowed(tx_pos,rx_pos,vertices,transparency,points)
%SHADOWED Calculate the magnitude of shadowed point

% Calculate the shadowing factor (0 - no shadowing, 1 - single pass, 2 -
% double pass)
V = tx_pos;
T = vertices.' - V.';

lambda = T\(points-V).';

shadowing_factor = and(all(lambda >= 0), sum(lambda)>=1).';

V = rx_pos;
T = vertices.' - V.';

lambda = T\(points-V).';

shadowing_factor = shadowing_factor + and(all(lambda >= 0), sum(lambda)>=1).';

S = transparency.^shadowing_factor;

end

