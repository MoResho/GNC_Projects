% Plot Earth
hold off;
[Xe, Ye, Ze] = sphere(100);
earth_radius = 6378;
% surf(Xe * earth_radius, Ye * earth_radius, Ze * earth_radius, 'FaceColor', [0.28, 0.56, 1]);
surf(Xe * earth_radius, Ye * earth_radius, -Ze * earth_radius, 'FaceColor', 'texturemap', 'CData', imread("EARTH.jpg"), 'EdgeColor', 'none');
hold on;
axis equal;
grid on;

% Define Orbital Elements
a = 22000;      % Semi-major axis
e = 0;          % Eccentricity
incl = 0;       % Inclination
RA = 0;         % Right Ascension of the Ascending Node (RAAN)
w = 0;          % Argument of Perigee
TA = 0;         % True Anomaly

% Convert Orbital Elements to Position and Velocity Vectors
[R0, V0] = COE2RV(a, e, TA, RA, incl, w);
X0 = [R0; V0];

% Calculate Orbital Period
mu = 398600;
T_orbit = 2 * pi * sqrt(a^3 / mu);

% Integrate the Equations of Motion
options = odeset('RelTol', 1e-8, 'AbsTol', 1e-10);
t_span = 0:10:T_orbit * 2;
[t, state] = ode45(@drfforbit, t_span, X0, options);

% Extract Position Components
X = state(:, 1);
Y = state(:, 2);
Z = state(:, 3);

% Plot the Orbital Trajectory
plot3(X, Y, Z, '.r');
axis off;
surf(500*Xe+X(end),500*Ye+Y(end),500*Ze+Z(end));

% Function to Calculate Derivatives of Position and Velocity
function ydot = drfforbit(t, y)
    r = y(1:3);      % Position vector [km]
    v = y(4:6);      % Velocity vector [km/s]
    mu = 398600;     % Earth's gravitational parameter [km^3/s^2]
    R = norm(r);     % Distance from Earth's center [km]

    % Acceleration
    a = -mu * r / R^3;
    
    % Derivative of state vector
    ydot = [v; a];
end

% Function to Convert Orbital Elements to Position and Velocity Vectors
function [R0, V0] = COE2RV(a, e, TA, RA, incl, w)
    % Convert degrees to radians for trigonometric functions
    RA_rad = deg2rad(RA);
    incl_rad = deg2rad(incl);
    w_rad = deg2rad(w);
    TA_rad = deg2rad(TA);
    
    % Gravitational parameter
    mu = 398600; % [km^3/s^2]
    
    % Semi-latus rectum
    p = a * (1 - e^2);
    
    % Specific angular momentum
    h = sqrt(mu * p);
    
    % Position in Perifocal Coordinate System
    rp = (p / (1 + e * cos(TA_rad)) ) * [cos(TA_rad); sin(TA_rad); 0];
    
    % Velocity in Perifocal Coordinate System
    vp = (mu / h) * [-sin(TA_rad); e + cos(TA_rad); 0];
    
    % Rotation Matrices
    % Rotation about Z-axis by RAAN
    R3_RAAN = [cos(RA_rad), -sin(RA_rad), 0;
               sin(RA_rad),  cos(RA_rad), 0;
               0,            0,           1];
           
    % Rotation about X-axis by Inclination
    R1_incl = [1, 0,             0;
              0, cos(incl_rad), -sin(incl_rad);
              0, sin(incl_rad),  cos(incl_rad)];
          
    % Rotation about Z-axis by Argument of Periapsis
    R3_arg = [cos(w_rad), -sin(w_rad), 0;
              sin(w_rad),  cos(w_rad), 0;
              0,           0,          1];
          
    % Combined Rotation Matrix: R = R3(RAAN) * R1(incl) * R3(argument of periapsis)
    Q = R3_RAAN * R1_incl * R3_arg;
    
    % Transform to ECI Frame
    R0 = Q * rp;
    V0 = Q * vp;
end