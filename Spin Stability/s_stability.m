% A case Study for Spin Stabilizaton of a CubeSat around z-axis

I1 = 0.1;  % Moment of inertia around x-axis
I2 = 0.2;  % Moment of inertia around y-axis
I3 = 1;  % Moment of inertia around z-axis

% Initial angular velocity
omega0 = [0.01; 0.01; 5];  % Initial spin about z-axis (major axis)

% Simulation time
tspan = [0 10];  

% Solve Euler equations
[t, omega] = ode45(@(t, w) eulerEquations(t, w, I1, I2, I3), tspan, omega0);

% Convert to Euler angles
euler_angles = zeros(length(t), 3);
euler_angles(1, :) = [0, 0, 0];

for i = 2:length(t)
    dt = t(i) - t(i-1);
    omega_body = omega(i, :)';
    
    % Convert angular velocity to Euler angle rates
    phi_dot = omega_body(1) + omega_body(2) * sin(euler_angles(i-1,1)) * tan(euler_angles(i-1,2)) + omega_body(3) * cos(euler_angles(i-1,1)) * tan(euler_angles(i-1,2));
    theta_dot = omega_body(2) * cos(euler_angles(i-1,1)) - omega_body(3) * sin(euler_angles(i-1,1));
    psi_dot = omega_body(2) * sin(euler_angles(i-1,1)) / cos(euler_angles(i-1,2)) + omega_body(3) * cos(euler_angles(i-1,1)) / cos(euler_angles(i-1,2));
    
    % Integrate Euler angle rates to get Euler angles
    euler_angles(i,1) = euler_angles(i-1,1) + phi_dot * dt;
    euler_angles(i,2) = euler_angles(i-1,2) + theta_dot * dt;
    euler_angles(i,3) = euler_angles(i-1,3) + psi_dot * dt;

    % Plotting the CubeSat Spinning around z_axis
    q = angle2quat(euler_angles(i,1), euler_angles(i,2), euler_angles(i,3), 'XYZ');
    poseplot(quaternion(q));
%     pause(0.1);
    drawnow;
%      exportgraphics(gcf,'CubeSat_Spin_major.gif','Append',true);
end

% Results
figure;
subplot(2,1,1);
plot(t, omega);
xlabel('Time (s)'); ylabel('\omega (rad/s)');
title('Angular Velocity'); legend('\omega_1', '\omega_2', '\omega_3');

subplot(2,1,2);
plot(t, euler_angles);
xlabel('Time (s)'); ylabel('Euler Angles (rad)');
title('Euler Angles'); legend('\phi (Roll)', '\theta (Pitch)', '\psi (Yaw)');

function dwdt = eulerEquations(~, w, I1, I2, I3)
    % angular velocities
    w1 = w(1);
    w2 = w(2);
    w3 = w(3);
    
    % Euler's equations
    dw1dt = ((I2 - I3) / I1) * w2 * w3;
    dw2dt = ((I3 - I1) / I2) * w3 * w1;
    dw3dt = ((I1 - I2) / I3) * w1 * w2;
    
    % Return w_derivatives
    dwdt = [dw1dt; dw2dt; dw3dt];
end
