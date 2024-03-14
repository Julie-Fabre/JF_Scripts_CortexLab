function cl_convert_allen_to_paxinos()


% Initial AP, ML, DV values (example)
AP = 1; % Example value
ML = 1; % Example value
DV = 1; % Example value

% Scaling factors
scale = [0.952, -1.031, 0.885] ./ 100; % [ML, AP, DV]

% Scale the coordinates
AP_scaled = AP * scale(2);
ML_scaled = ML * scale(1);
DV_scaled = DV * scale(3);

% Tilt angle
theta_degrees = 15; % Tilt by 15 degrees
theta_radians = deg2rad(theta_degrees);

% Rotation matrix for tilting around ML axis
R = [cos(theta_radians) 0 sin(theta_radians); 0 1 0; -sin(theta_radians) 0 cos(theta_radians)];

% Apply rotation to the scaled AP and DV (ignoring ML in rotation)
coordinates_scaled = [AP_scaled; ML_scaled; DV_scaled];
coordinates_rotated = R * coordinates_scaled;

% Output the rotated coordinates
disp('Rotated Coordinates:');
disp(coordinates_rotated);

end