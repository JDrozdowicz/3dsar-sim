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

function [points] = face2points(vertices,d,magnitude,phase)
%FACE2POINTS Convert 3D face to 3D points

% Calculate the normal vector of the face
normalv = cross(diff(vertices(1:2,:)),diff(vertices(2:3,:)));
N = normalv/norm(normalv);

% Calculate the central point of the face
centerp = mean(vertices);

% Move face to the origin
verticesc = vertices - centerp;

%Calculate rotation matrix R*N = [0;0;1]
% https://math.stackexchange.com/questions/180418/calculate-rotation-matrix-to-align-vector-a-to-vector-b-in-3d/476311#476311
v = cross(N.',[0;0;1]);
s = norm(v);
c = dot(N.',[0;0;1]);

vx = [0, -v(3), v(2); v(3), 0, -v(1); -v(2), v(1), 0];

R = eye(3) + vx + vx^2 .*((1-c)/s^2);

%Rotate vertices
verticesr = (R*verticesc.').';

verticesp = verticesr(:,1:2);


vpx = verticesp(:,1);
vpy = verticesp(:,2);

xmax = max(vpx);
xmin = min(vpx);
ymax = max(vpy);
ymin = min(vpy);


area = 0.5*abs(det([verticesp.';[1,1,1]]));

npoints = ceil(area/d.^2);

points = zeros(npoints,3);

iter_points = 0;

while iter_points < npoints
    pointx = xmin + rand().*(xmax-xmin);
    pointy = ymin + rand().*(ymax-ymin);
    
    T = (verticesp(1:2,:) - verticesp(3,:)).';
    
    lambda = T\(([pointx,pointy]-verticesp(3,:)).');
    
    if all(lambda > 0) && (sum(lambda) < 1)
        iter_points = iter_points + 1;
        points(iter_points,:) = transpose(R)*[pointx;pointy;0]+centerp.';
    end


    
end

points = [points, magnitude.*ones(size(points,1),1), phase.*ones(size(points,1),1)];

end

