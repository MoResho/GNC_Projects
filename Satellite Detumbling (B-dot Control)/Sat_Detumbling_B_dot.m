clc;
clear;

%% --- Simulation Parameters ---
dt = 1;                    % Time step [s]
T = 100000;                 % Total simulation time [s]
N = T/dt;                  % Number of steps
time = 0:dt:T;             % Time vector

%% --- Satellite Parameters ---
I = diag([0.01, 0.01, 0.02]);     % Inertia matrix [kg·m²]
invI = inv(I);                   
k_bdot = -1e3;                    % B-dot gain
m_max = 1;                     % Max dipole [A·m²]

%% --- Initial Conditions ---
w = deg2rad([1; -1; 0.5]);        % Angular velocity [rad/s]
q = [1; 0; 0; 0];                 % Quaternion [w x y z] (scalar-first)

%% --- Orbit Parameters (example location) ---
lat = 0;                         % Latitude [deg]
lon = 0;                         % Longitude [deg]
alt = 500;                       % Altitude [km]

%% --- Preallocation ---
w_hist = zeros(3, N+1);
q_hist = zeros(4, N+1);
B_body_hist = zeros(3, N+1);

w_hist(:,1) = w;
q_hist(:,1) = q;

startDate = datetime(2024, 6, 7, 12, 0, 0);  % UTC start time

%% --- Main Simulation Loop ---
for i = 1:N
    % Get current magnetic field in ECI
    t_now = startDate + seconds((i-1)*dt);
    decimalYear = year(t_now) + (day(t_now, 'dayofyear') - 1)/365;
    [Bx, By, Bz] = igrfmagm(alt, lat, lon, decimalYear);
    B_eci = 1e-9 * [Bx(1); By(1); Bz(1)];  % Convert nT to T

    % Rotate B field into body frame
    C = quat2dcm(q');                   
    B_body = C * B_eci;
    B_body_hist(:,i) = B_body;

    % Estimate dB/dt
    if i == 1
        dB = zeros(3,1);
    else
        dB = (B_body - B_body_hist(:, i-1)) / dt;
    end

    % B-dot control law with saturation
    m = k_bdot * dB;
    if norm(m) > m_max
        m = m_max * m / norm(m);
    end

    % Dynamics update
    Tm = cross(m, B_body);                     % Magnetic torque
    dw = invI * (Tm - cross(w, I*w));          % Angular acceleration
    w = w + dw * dt;                           % Update angular velocity

    % Quaternion update
    Omega = [  0,   -w(1), -w(2), -w(3);
              w(1),   0,   w(3), -w(2);
              w(2), -w(3),  0,    w(1);
              w(3),  w(2), -w(1),  0  ];
    q = q + 0.5 * Omega * q * dt;
    q = q / norm(q);                           % Normalize

    % Store history
    w_hist(:,i+1) = w;
    q_hist(:,i+1) = q;
end

%% --- Plot Angular Velocity ---
figure;
plot(time, w_hist(1,:), 'b', 'DisplayName', '\omega_x'); hold on;
plot(time, w_hist(2,:), 'r', 'DisplayName', '\omega_y');
plot(time, w_hist(3,:), 'g', 'DisplayName', '\omega_z');
xlabel('Time [s]');
ylabel('Angular Velocity [rad/s]');
title('Angular Velocity vs Time (B-dot Control)');
legend; grid on;

%% --- Visualize Attitude with poseplot ---
figure;
axis equal;
for i = 1:100:length(q_hist)
    poseplot(quaternion(q_hist(:,i)'));
    title(sprintf('Time = %d s', time(i)));
    drawnow;
%     exportgraphics(gcf,'Satellite_Detumbling.gif','Append',true);
end
