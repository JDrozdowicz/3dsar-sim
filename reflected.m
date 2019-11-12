function [S] = reflected(tx_pos,rx_pos,vertices,roughness)
%REFLECTED Calculate the magnitude of reflected energy

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

