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

clear all;
close all;

% Radar parameters
% carrier frequency [Hz]
fc = 24e9; 
% bandwidth [Hz]
B = 1000e6; 
% pulse duration [s]
T = 1e-3;

% Constants
% speed of light [m/s]
c = 3e8;

% Imaging parameters
% start range [m]
start_range = 20;
% end range [m]
end_range = 120;
% range cell size [m]
cell_size = c/(2*B)/5;
% face sampling density [m]
sampling_density = cell_size;

% Scene definition
% Scene can be defined by hand (as in this example), but also imported from
% an external tool.
% Scene consists of points and faces defined by vertices
% Reflecting points, not belonging to any face (x,y,z,magnitude,phase)
point = [0,0,0,0,0];
% Load vertices  and reflecting faces
load('complex_scene');

% Antenna trajectory definition
% Antenna positions [x,y,z,dir_x,dir_y,dir_z,rotation)
[ant_rot(1), ant_rot(2), ant_rot(3)] = sph2cart(deg2rad(0),deg2rad(-30),1);
tx_pos = [-50,0,50,ant_rot,0];
rx_pos = tx_pos;


% Antenna beam pattern definition (azimuth 30 deg BW, elevation 60 deg BW)
ant_pat = @(az,el) abs(sinc(2*az/deg2rad(30)).^0.5.*sinc(2*el/deg2rad(60)).^0.5);

% END OF INPUT PARAMETERS SECTION

for iter_faces = 1:numel(face)
    % Convert faces to points
    points{iter_faces} = face2points(vertex(face(iter_faces).v,:),sampling_density, face(iter_faces).magnitude, face(iter_faces).phase);

end


% MAIN LOOP

raw_data = zeros(size(tx_pos,1),ceil((end_range-start_range)/cell_size));

for iter_pos = 1:size(tx_pos,1)
    
    % those points will be modified
    curr_points = points;
    curr_point = point;
    
    for iter_faces = 1:numel(face)
        % Calculate reflectivity
        curr_points{iter_faces}(:,4) = curr_points{iter_faces}(:,4).*reflected(tx_pos(iter_pos,1:3),rx_pos(iter_pos,1:3),vertex(face(iter_faces).v,:),face(iter_faces).roughness);
        % Calculate shadowing for points belonging to faces
        for iter_faces2 = 1:numel(face)
            if(iter_faces ~= iter_faces2)
                curr_points{iter_faces}(:,4) = curr_points{iter_faces}(:,4).*shadowed(tx_pos(iter_pos,1:3),rx_pos(iter_pos,1:3),vertex(face(iter_faces2).v,:),face(iter_faces2).transparency,curr_points{iter_faces}(:,1:3));
            end
        end
        % Calculate shadowing for points not belonging to faces
        curr_point(:,4) = curr_point(:,4).*shadowed(tx_pos(iter_pos,1:3),rx_pos(iter_pos,1:3),vertex(face(iter_faces).v,:),face(iter_faces).transparency,curr_point(:,1:3));
    
    end
    
    for iter_faces = 1:numel(face)
        curr_point = [curr_point; curr_points{iter_faces}];
    end
    
    % Calculate antenna gain factor for points
    ant_gain = ones(size(curr_point,1),1);
    % Include TX antenna pattern
    for iter_ant = 1:size(curr_point,1)
        [az, el] = ant_orientation(tx_pos(iter_pos,:),curr_point(iter_ant,1:3).');
        ant_gain(iter_ant) = ant_gain(iter_ant).*ant_pat(az,el);
    end
    % Include RX antenna pattern
    for iter_ant = 1:size(curr_point,1)
        [az, el] = ant_orientation(rx_pos(iter_pos,:),curr_point(iter_ant,1:3).');
        ant_gain(iter_ant) = ant_gain(iter_ant).*ant_pat(az,el);
    end
    curr_point(:,4) = curr_point(:,4).*ant_gain;
    
    % Calculate response for points
    raw_data(iter_pos,:) = sim_resp(fc, B, T, cell_size, start_range, end_range, curr_point, tx_pos(iter_pos,:), rx_pos(iter_pos,:));
    
end

% Plot scene
figure;
% Plot faces
hold on;
for iter_faces = 1:numel(face)
    fill3(vertex(face(iter_faces).v,1),vertex(face(iter_faces).v,2),vertex(face(iter_faces).v,3), face(iter_faces).magnitude, 'FaceAlpha', 1 - face(iter_faces).transparency, 'EdgeColor', 'k', 'EdgeAlpha', 0.2);
end
% %Plot points derived from faces:
% for iter_faces = 1:numel(face)
%     plot3(points{iter_faces}(:,1),points{iter_faces}(:,2),points{iter_faces}(:,3),'.k');
% end
% Plot trajectory
plot3(tx_pos(:,1),tx_pos(:,2),tx_pos(:,3),'r-x');
plot3(rx_pos(:,1),rx_pos(:,2),rx_pos(:,3),'b-o');
% Plot visible points
visible = curr_point(:,4) > eps;
vispoint = curr_point(visible,:);
vstep = 100;
% scatter3(curr_point(visible,1),curr_point(visible,2),curr_point(visible,3),1+50*curr_point(visible,4),curr_point(visible,4),'filled','MarkerEdgeColor','none');
scatter3(vispoint(1:vstep:end,1),vispoint(1:vstep:end,2),vispoint(1:vstep:end,3),1+20*vispoint(1:vstep:end,4),'k','filled','MarkerEdgeColor','none','MarkerFaceAlpha',0.1);

hold off;

axis equal;
xlabel('x [m]');
ylabel('y [m]');
zlabel('z [m]');
view(3);

figure;
plot(start_range:cell_size:end_range,abs(raw_data));
xlabel('x [m]');
ylabel('y [m]');

