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

function [S] = reflected(tx_pos,rx_pos,vertices,roughness)
%REFLECTED Calculate the magnitude of reflected energy

% based on https://www.scratchapixel.com/lessons/3d-basic-rendering/phong-shader-BRDF

% calculate the roughness factor
roughness_factor = (tan(pi/2-min(max(roughness,0),1)*pi/2));

% Calculate the normal vectors of the face
normalv = cross(diff(vertices(1:2,:)),diff(vertices(2:3,:)));
N = normalv/norm(normalv);
% Calculate the central point of the face
centerp = mean(vertices);


%calculate the directional vectors relative to the antenna
L = centerp - tx_pos;
Ln = L/norm(L);

V = centerp - rx_pos;
Vn = V/norm(V);

% Phong
R = -Ln + 2*dot(N, Ln)*N;
Rn = R/norm(R);

S = max(0,dot(Rn,Vn)).^roughness_factor;

end

