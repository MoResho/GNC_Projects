%% -------------------- Initial Conditions --------------------
V = 30;                         % Velocity of the aircraft (m/s)
theta0 = 2.1471;                % Initial pitch angle (deg.)
u0 = V * cosd(theta0);         % Initial forward velocity (m/s)
w0 = V * sind(theta0);         % Initial vertical velocity (m/s)
theta0_rad = deg2rad(theta0);  % Initial pitch angle (rad)
q0 = 0;                         % Initial pitch rate (rad/s)
x0 = 0;                         % Initial x-position (m)
z0 = 0;                         % Initial z-position (m)

% State vector [u, w, theta, q, x, z]
initial_conditions = [u0 w0 theta0_rad q0 x0 z0]';

%% -------------------- Time Span --------------------
tspan = [0 500];                % Simulate for 500 seconds

%% -------------------- ODE Solver --------------------
[t, states] = ode45(@EOM, tspan, initial_conditions);

%% -------------------- Extract States --------------------
u     = states(:,1);
w     = states(:,2);
theta = states(:,3);
q     = states(:,4);
x     = states(:,5);
z     = states(:,6);
alpha = atan2(w, u);

%% -------------------- Plot Velocities --------------------
figure(1);
subplot(2,1,1);
plot(t, u, 'b', 'linewidth', 1.5); grid on;
title('Forward velocity (u)');
xlabel('Time (s)');
ylabel('u (m/s)');

subplot(2,1,2);
plot(t, w, 'r', 'linewidth', 1.5); grid on;
title('Vertical velocity (w)');
xlabel('Time (s)');
ylabel('w (m/s)');

%% -------------------- Plot Angles --------------------
figure(2);
subplot(2,1,1);
plot(t, theta * (180/pi), 'm', 'linewidth', 1.5); grid on;
title('Pitch angle (\theta)');
xlabel('Time (s)');
ylabel('\theta (deg)');

subplot(2,1,2);
plot(t, alpha * (180/pi), 'c', 'linewidth', 1.5); grid on;
title('AoA (\alpha)');
xlabel('Time (s)');
ylabel('\alpha (deg)');

%% -------------------- Plot Pitch Rate & Altitude --------------------
figure(3);
subplot(2,1,1);
plot(t, q * (180/pi), 'g', 'linewidth', 1.5); grid on;
title('Pitch rate (q)');
xlabel('Time (s)');
ylabel('q (deg/s)');

subplot(2,1,2);
plot(t, -z, 'k', 'linewidth', 1.5); grid on;
title('Altitude (x-z)');
xlabel('Time (s)');
ylabel('z (m)');

%% -------------------- Equations of Motion --------------------
function dydt = EOM(t, y)
    % Constants
    m = 13.5; Iyy = 1.135; g = 9.81;
    rho = 1.225; S = 0.55; c = 0.19;
    T_max = 2 * g; delta_t = 0.5;

    % Aerodynamic coefficients
    CL0 = 0.28; CL_alpha = 3.45; CL_delta_e = -0.36;
    CD0 = 0.03; CD_alpha = 0.3;
    CM0 = -0.024; CM_alpha = -0.38; CM_delta_e = -0.5;

    % State extraction
    u = y(1); w = y(2); theta = y(3); q = y(4); x = y(5); z = y(6);
    V = sqrt(u^2 + w^2);
    alpha = atan2(w, u);
    delta_e = -4.3791 * (pi/180);

    % Aerodynamic forces and moments
    CL = CL0 + CL_alpha * alpha + CL_delta_e * delta_e;
    CD = CD0 + CD_alpha * alpha;
    CM = CM0 + CM_alpha * alpha + CM_delta_e * delta_e;

    L = 0.5 * rho * V^2 * S * CL;
    D = 0.5 * rho * V^2 * S * CD;
    My = 0.5 * rho * V^2 * S * c * CM;
    T = T_max * delta_t;

    % Forces in body axes
    Fx = -D * cos(alpha) + L * sin(alpha) + T;
    Fz = -L * cos(alpha) - D * sin(alpha);

    % Accelerations
    A_xe = Fx/m - g * sin(theta);
    A_ze = Fz/m + g * cos(theta);

    % State derivatives
    du     = A_xe - q * w;
    dw     = A_ze + q * u;
    dtheta = q;
    dq     = My / Iyy;
    dx     = u * cos(theta) + w * sin(theta);
    dz     = -u * sin(theta) + w * cos(theta);

    dydt = [du dw dtheta dq dx dz]';
end
