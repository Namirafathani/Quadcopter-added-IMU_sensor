%% Prepare workspace
clear; close all; clc;

% select waypoint type

% select trajectory type
% 1 - Circular 2 - spiral 3 - eight
% 4 - rectangular
param.traj = 2;

% Time change
% 1 - 1 * T, 2 - 2 * T, 3 - 3 * T
param.time_change = 3.1;

%% Quadcopter Variables
param.g = 9.81; 
param.m = 0.468;
param.l = 0.225;
param.Jr = 3.357e-5;
param.Kt = 2.98e-6;
param.Kb = 1.140e-7;
Kdx = 0.25; Kdy = 0.25; Kdz = 0.25;
param.Kd = diag([Kdx; Kdy; Kdz]);
Ixx = 4.856e-3; Iyy = 4.856e-3; Izz = 8.801e-3;
param.I = diag([Ixx; Iyy; Izz]);

%% Platform Parameters
platform_length = 0.5;      % Length of the platform (x-direction)
platform_width = 0.5;       % Width of the platform (y-direction)
platform_amplitude = 3.0;   % Oscillation amplitude (meters)
platform_frequency = 0.05;   % Oscillation frequency (Hz)

% For Nonlinear PD with helix trajectory
if param.traj == 2
    PI_var.Kp = 6; PI_var.Ki = 1;
    PI_var.Kp_z = 0.5; PI_var.Ki_z = 0.07;
    PI_var.Kp_xy = 0.09; PI_var.Ki_xy = 0.01;
    
    PI_var.Kp_x = 0.11; PI_var.Ki_x = 0.02;
    PI_var.Kp_y = 0.11; PI_var.Ki_y = 0.01;

elseif param.traj == 1
    PI_var.Kp = 9.2; PI_var.Ki = 1.5;
    PI_var.Kp_z = 0.44; PI_var.Ki_z = 1;
    PI_var.Kp_xy = 0.102; PI_var.Ki_xy = 0.01;
end

if param.traj == 3
    PI_var.Kp = 6.0; PI_var.Ki = 1.5;
    PI_var.Kp_z = 5; PI_var.Ki_z = 1;
    PI_var.Kp_xy = 0.20; PI_var.Ki_xy = 0.01;

    PI_var.Kp_x = 0.14; PI_var.Ki_x = 0.02;
    PI_var.Kp_y = 0.19; PI_var.Ki_y = 0.01;

elseif param.traj == 4
    PI_var.Kp = 10.0; PI_var.Ki = 15;
    PI_var.Kp_z = 4; PI_var.Ki_z = 1;
    PI_var.Kp_xy = 0.04; PI_var.Ki_xy = 0.01;
end

%% Initialize Simulation
start_time = 0;
T = 15; % in seconds
end_time = param.time_change * T;
dt = 0.001;
sim_time = (start_time:dt:end_time)';
n_sim = length(sim_time);

xi =  zeros(n_sim, 3);
xi_desired =  zeros(n_sim, 3);
nu_I = zeros(n_sim, 3);
eta = zeros(n_sim, 3);
omega_I = zeros(n_sim, 3);
control_values = zeros(n_sim, 4);

% Preallocate IMU data storage
imu_gyro_data = zeros(n_sim, 3);    % Angular velocity (body frame)
imu_accel_data = zeros(n_sim, 3);   % Linear acceleration (body frame)

%% Simulation
for idx_t = 1:n_sim

    % extract current time step data
    xi_current = xi(idx_t, :)';
    nu_I_current = nu_I(idx_t, :)';
    eta_current = eta(idx_t, :)';
    omega_I_current = omega_I(idx_t, :)';

    % desired position is a circular or spiral trajectory
    if param.traj == 1 
       [xi_des, xi_des_d, xi_des_dd] = circular_traj(sim_time(idx_t), T);
    elseif param.traj == 2
       [xi_des, xi_des_d, xi_des_dd] = spiral_traj(sim_time(idx_t), T);
    end

    if param.traj == 3
        [xi_des, xi_des_d, xi_des_dd] = eight_traj(sim_time(idx_t), T);
    elseif param.traj == 4
        [xi_des, xi_des_d, xi_des_dd] = rectangular_traj(idx_t, T);
    end

    % convert from inertia frame to body frame
    omega_B_current = inertia2body(omega_I_current, eta_current);
    
    control_input = PI(omega_I_current, eta_current, xi_current, nu_I_current, param, PI_var, xi_des, xi_des_d, xi_des_dd); % omega square

    nu_I_dot = linear_acceleration(control_input, eta_current, nu_I_current, param);
    
    % Compute IMU Measurement in Simulation Loop
    % Compute IMU measurements (current timestep)
    R = rotation(eta_current);  % Rotation matrix (body to inertial)
    gravity_inertial = [0; 0; -param.g];  % Gravity vector in inertial frame

    % Non-gravitational acceleration in inertial frame
    non_grav_accel_inertial = nu_I_dot - gravity_inertial;

    % Convert to body frame (IMU measurements)
    imu_accel_body = R' * non_grav_accel_inertial;  % Specific force
    imu_gyro_body = omega_B_current;               % Angular rates

    % Store IMU data
    imu_gyro_data(idx_t, :) = imu_gyro_body';
    imu_accel_data(idx_t, :) = imu_accel_body';
   
    omega_B_dot = rotational_acceleration(control_input, omega_B_current, param);
    
    % linear motion evolution using discretized Euler
    nu_I_next = nu_I_current + dt * nu_I_dot;
    xi_next = xi_current + dt * nu_I_next;

    % rotational motion evolution using discretized Euler
    omega_B_next = omega_B_current + dt * omega_B_dot;
    % convert from bodyframe to inertia frame
    omega_I_next = body2inertia(omega_B_next, eta_current);
    eta_next = eta_current + dt * omega_I_next;

    % fill in next time step data
    xi(idx_t + 1, :) = xi_next';
    nu_I(idx_t + 1, :) = nu_I_next';
    eta(idx_t + 1, :) = eta_next';
    omega_I(idx_t + 1, :) = omega_I_next';
    control_values(idx_t, :) = sqrt(control_input)';
    
    xi_desired(idx_t, :) = xi_des';
end

%% 3d Animation Plots
figure;
ax = gca;
hold(ax, 'on');
grid(ax, 'on');
view(ax, 3);  % 3D view
xlabel(ax, 'X Position');
ylabel(ax, 'Y Position');
zlabel(ax, 'Z Position');
title(ax, 'Quadcopter Simulation - X Configuration');

all_data = [xi; xi_desired];
max_vals = max(all_data);
min_vals = min(all_data);
margin = 0.1 * (max_vals - min_vals);
axis(ax, [min_vals(1)-margin(1), max_vals(1)+margin(1), ...
          min_vals(2)-margin(2), max_vals(2)+margin(2), ...
          min_vals(3)-margin(3), max_vals(3)+margin(3)]);
daspect(ax, [1 1 1]);  % Equal aspect ratio

% Create moving platform (oscillating in x)
[X_plat, Y_plat] = meshgrid(...
    linspace(-platform_length/2, platform_length/2, 10), ...
    linspace(-platform_width/2, platform_width/2, 10));
Z_plat = zeros(size(X_plat));
platform = surf(X_plat, Y_plat, Z_plat, ...
    'FaceColor', [1.0 0.71 0.76], ...
    'EdgeColor', 'none', ...
    'FaceAlpha', 0.5, ...
    'DisplayName','Landing Platform');

% Initialize plot objects
actual = plot3(ax, nan, nan, nan, 'b-', 'LineWidth', 1.5, 'DisplayName','Actual Path');
desired = plot3(ax, nan, nan, nan, 'r--', 'LineWidth', 1.5,'DisplayName','Desired Path');
actual_point = plot3(ax, nan, nan, nan, 'bo', 'MarkerSize', 8, 'MarkerFaceColor', 'b','HandleVisibility','off');
desired_point = plot3(ax, nan, nan, nan, 'rs', 'MarkerSize', 8, 'MarkerFaceColor', 'none','HandleVisibility','off');
legend('show', 'Location', 'best');

% Optional: Bring platform to bottom layer
uistack(platform, 'bottom')

% Animation loop
for i = 1:n_sim
    % --- Platform Motion (Sinusoidal in x) ---
    x_platform = platform_amplitude * sin(2*pi*platform_frequency*sim_time(i));
    
    % Update platform position
    set(platform, 'XData', X_plat + x_platform);



    if mod(i,400) == 1 || i == n_sim
    % Update actual trajectory (up to current point)
    set(actual, 'XData', xi(1:i,1), 'YData', xi(1:i,2), 'ZData', xi(1:i,3));
    
    % Update desired trajectory (up to current point)
    set(desired, 'XData', xi_desired(1:i,1), 'YData', xi_desired(1:i,2), 'ZData', xi_desired(1:i,3));
    
    % Update current position markers
    set(actual_point, 'XData', xi(i,1), 'YData', xi(i,2), 'ZData', xi(i,3));
    set(desired_point, 'XData', xi_desired(i,1), 'YData', xi_desired(i,2), 'ZData', xi_desired(i,3));

    drawnow;
    end
end


%% Plots
figure
hold on
plot(sim_time, control_values(:, 1), 'LineWidth', 2)
plot(sim_time, control_values(:, 2), 'LineWidth', 2)
plot(sim_time, control_values(:, 3), 'LineWidth', 2)
plot(sim_time, control_values(:, 4), 'LineWidth', 2)
hold off
xlabel('Time [s]')
ylabel('control input \omega[rad/s]')
legend('\omega_1', '\omega_2', '\omega_3', '\omega_4', 'Location', 'Best')
grid on

figure
hold on
plot(sim_time, xi(1:end-1, 1), 'LineWidth', 2)
plot(sim_time, xi_desired(:, 1), 'LineWidth', 2)
hold off
xlabel('Time [s]')
ylabel('position - x [m]')
legend('x actual', 'x desired', 'Location', 'Best')
grid on

figure
hold on
plot(sim_time, xi(1:end-1, 2), 'LineWidth', 2)
plot(sim_time, xi_desired(:, 2), 'LineWidth', 2)
hold off
xlabel('Time [s]')
ylabel('position - y [m]')
legend('y actual', 'y desired', 'Location', 'Best')
grid on

figure
hold on
plot(sim_time, xi(1:end-1, 3), 'LineWidth', 2)
plot(sim_time, xi_desired(:, 3), 'LineWidth', 2)
hold off
xlabel('Time [s]')
ylabel('position - z [m]')
legend('z actual', 'z desired', 'Location', 'Best')
grid on

figure
hold on
plot3(xi(1:end-1, 1), xi(1:end-1, 2), xi(1:end-1, 3), 'LineWidth', 2)
plot3(xi_desired(1:end-1, 1), xi_desired(1:end-1, 2), xi_desired(1:end-1, 3),'--', 'LineWidth', 2)
plot3(xi_desired(1,1), xi_desired(1,2), xi_desired(1,3),'o', 'LineWidth', 6,'MarkerSize', 6)
plot3(xi(1,1), xi(1,2), xi(1,3),'*', 'LineWidth', 6,'MarkerSize', 6)
    
xlabel('x [m]')
ylabel('y [m]')
zlabel('z [m]')
grid on
grid minor
legend('Actual Trajectory','Desired Trajectory','Start of desired trajectory','Start of actual trajectory')

%% IMU Sensor Plots
% Angular Velocity (Gyroscope)
figure;
subplot(3,1,1);
plot(sim_time, imu_gyro_data(:,1), 'LineWidth', 1.5);
ylabel('p [rad/s]');
title('IMU Gyroscope (Body Frame)');
grid on;

subplot(3,1,2);
plot(sim_time, imu_gyro_data(:,2), 'LineWidth', 1.5);
ylabel('q [rad/s]');
grid on;

subplot(3,1,3);
plot(sim_time, imu_gyro_data(:,3), 'LineWidth', 1.5);
xlabel('Time [s]');
ylabel('r [rad/s]');
grid on;

% Linear Acceleration (Accelerometer)
figure;
subplot(3,1,1);
plot(sim_time, imu_accel_data(:,1), 'LineWidth', 1.5);
ylabel('a_x [m/s^2]');
title('IMU Accelerometer (Body Frame)');
grid on;

subplot(3,1,2);
plot(sim_time, imu_accel_data(:,2), 'LineWidth', 1.5);
ylabel('a_y [m/s^2]');
grid on;

subplot(3,1,3);
plot(sim_time, imu_accel_data(:,3), 'LineWidth', 1.5);
xlabel('Time [s]');
ylabel('a_z [m/s^2]');
grid on;

%% Functions
function nu_I_dot = linear_acceleration(control_input, eta, nu_I, param)
    % Compute translational/linear motion using Newtonian method
    % gravity
    nu_I_dot_1 = -param.g * [0, 0, 1]';
    R = rotation(eta);
    % total thrust
    T_B = [0, 0, param.Kt * sum(control_input)]';
    nu_I_dot_2 =  (1 / param.m) * R * T_B;
    % drag force
    nu_I_dot_3 = -(1 / param.m) * param.Kd * nu_I;
    
    nu_I_dot = nu_I_dot_1 + nu_I_dot_2 + nu_I_dot_3;
end

function omega_B_dot = rotational_acceleration(control_input, omega_B, param)
    % Compute rotational motion using Euler method
    tau = torques(control_input, param);
    gamma = computeGamma(control_input, omega_B, param);
    % using direct representation for ease of computation
    omega_B_dot = inv(param.I) * (tau - gamma + cross(omega_B, param.I * omega_B));
end

function gamma = computeGamma(control_input, omega_B, param)
    % compute the gyroscopic forces causes by combined rotation of motors
    omega = sqrt(control_input);
    omega_gamma = -omega(1) + omega(2) - omega(3) + omega(4);
    gamma = param.Jr * cross(omega_B, [0; 0; 1]) * omega_gamma;
end

% Compute torques, given current inputs, length, drag coefficient, and thrust coefficient.
function tau = torques(control_input, param)
    % Inputs are values for ${\omega_i}^2$
    tau = [
        param.l * param.Kt * (-control_input(2) + control_input(4))
        param.l * param.Kt * (-control_input(1) + control_input(3))
        param.Kb * (-control_input(1) + control_input(2) - control_input(3) + control_input(4))
    ];
end

function R = rotation(eta)
    % Compute rotational matrix
    phi = eta(1); theta = eta(2); psi = eta(3);
    
    R11 = cos(psi) * cos(theta);
    R12 = cos(psi) * sin(theta) * sin(phi) - sin(psi) * cos(phi);
    R13 = cos(psi) * sin(theta) * cos(phi) + sin(psi) * sin(phi);
    
    R21 = sin(psi) * cos(theta);
    R22 = sin(psi) * sin(theta) * sin(phi) + cos(psi) * cos(phi);
    R23 = sin(psi) * sin(theta) * cos(phi) - cos(psi) * sin(phi);
    
    R31 = -sin(theta);
    R32 = cos(theta) * sin(phi);
    R33 = cos(theta) * cos(phi);
    
    R = [R11 R12 R13;...
         R21 R22 R23;...
         R31 R32 R33];
end

function I2B = inertia2body(omega_I, eta)
    phi = eta(1); theta = eta(2); psi = eta(3);
    W = [1   0        -sin(theta);...
         0   cos(phi)  cos(theta) * sin(phi);...
         0  -sin(phi)  cos(theta) * cos(phi)];
    I2B = W * omega_I;
end

function B2I = body2inertia(omega_B, eta)
    phi = eta(1); theta = eta(2); psi = eta(3);
    W_inv = [1   sin(phi) * tan(theta)    cos(phi) * tan(theta);...
             0   cos(phi)                -sin(phi);...
             0   sin(phi) / cos(theta)    cos(phi) / cos(theta)];
    B2I = W_inv * omega_B;  
end

%% Controls
function ctrl = PI(omega_I, eta, xi, nu_I, param, PI_var, xi_des, xi_des_d, xi_des_dd)
    ctrl = zeros(4, 1);
    int_err = zeros(4,1);
    % Unload parameters
    phi = eta(1); theta = eta(2); psi = eta(3);
    phi_d = omega_I(1); theta_d = omega_I(2); psi_d = omega_I(3);
    z = xi(3); z_d = nu_I(3);

    int_z = int_err(1);
    int_phi = int_err(2);
    int_theta = int_err(3);
    int_psi = int_err(4);

    % desired position computation - x and y desired translated to phi and theta
    z_des = xi_des(3); z_des_d = xi_des_d(3);
    % desired angle
    eta_des = pos2ang(xi, xi_des, xi_des_d, xi_des_dd, nu_I, param, PI_var);
    % desired angle derivative
    eta_des_d = zeros(size(eta));

    % unload values
    phi_des = eta_des(1); phi_des_d = eta_des_d(1); 
    theta_des = eta_des(2); theta_des_d = eta_des_d(2);
    psi_des = eta_des(3); psi_des_d = eta_des_d(3);

    %T = (param.g + PD_var.Kd_z * (z_des_d - z_d) + PD_var.Kp_z * (z_des - z)) * ( param.m / (cos(phi) * cos(theta)) );
    %tau_phi = (PD_var.Kd * (phi_des_d - phi_d) + PD_var.Kp * (phi_des - phi)) * param.I(1,1);
    %tau_theta = (PD_var.Kd * (theta_des_d - theta_d) + PD_var.Kp * (theta_des - theta)) * param.I(2,2);
    %tau_psi = (PD_var.Kd * (psi_des_d - psi_d) + PD_var.Kp * (psi_des - psi)) * param.I(3,3);
    
    int_z = int_z + (z_des_d - z_d);
    int_phi = int_phi + (phi_des_d - phi_d);
    int_theta = int_theta + (theta_des_d - theta_d);
    int_psi = int_psi + (psi_des_d - psi_d);

    T = (param.g + PI_var.Kp_z * (z_des - z) + PI_var.Ki_z * (int_z)) * (param.m / (cos(phi) * cos(theta)));
    tau_phi = (PI_var.Kp * (phi_des - phi) + PI_var.Ki * (int_phi)) * param.I(1,1);
    tau_theta = (PI_var.Kp * (theta_des - theta) + PI_var.Ki * (int_theta)) * param.I(2,2);
    tau_psi = (PI_var.Kp * (psi_des - psi) + PI_var.Ki * (int_psi)) * param.I(3,3);

    ctrl(1) = ( T / (4 * param.Kt) ) - ( tau_theta / (2 * param.Kt * param.l) ) - ( tau_psi / (4 * param.Kb) );
    ctrl(2) = ( T / (4 * param.Kt) ) - ( tau_phi / (2 * param.Kt * param.l) ) + ( tau_psi / (4 * param.Kb) );
    ctrl(3) = ( T / (4 * param.Kt) ) + ( tau_theta / (2 * param.Kt * param.l) ) - ( tau_psi / (4 * param.Kb) );
    ctrl(4) = ( T / (4 * param.Kt) ) + ( tau_phi / (2 * param.Kt * param.l) ) + ( tau_psi / (4 * param.Kb) );
end

%% Desired Trajectory
% position trajectory
function [xi_des, xi_des_d, xi_des_dd]  = circular_traj(t, T)
    radius = 3;
    ang = 2 * pi * t/T;

    xi_des = zeros(3, 1);
    xi_des_d = zeros(3, 1);
    xi_des_dd = zeros(3, 1);

    % define circular trajectory for the x and y direction
    xi_des(1) = radius * cos(ang);
    xi_des(2) = radius * sin(ang);
    xi_des(3) = t/T; %  equation of line

    % define derivative
    xi_des_d(1) = 0;%-2 * pi * radius * sin(ang)/T;
    xi_des_d(2) = 0;%2 * pi * radius * cos(ang)/T;
    xi_des_d(3) = 0;%1/T;

    % define double derivative
    xi_des_dd(1) = 0;%-4 * radius * pi * pi * cos(ang) / (T * T);
    xi_des_dd(2) = 0;%-4 * radius * pi * pi * sin(ang) / (T * T);
    xi_des_dd(3) = 0;
end

function [xi_des, xi_des_d, xi_des_dd]  = spiral_traj(t, T)
    radius = 3;
    ang = 2* pi * t/T;
    h = 5;
    
    xi_des = zeros(3, 1);
    xi_des_d = zeros(3, 1);
    xi_des_dd = zeros(3, 1);
    % define circular trajectory for the x and y direction
    xi_des(1) = radius * cos(ang);
    xi_des(2) = h * t/T;
    xi_des(3) = t/T; %  equation of line
    % define derivative
    xi_des_d(1) = 0; %-2 * pi * radius * sin(ang)/T;
    xi_des_d(2) = 0; %h / T;
    xi_des_d(3) = 0; %1/T;
    % define double derivative
    xi_des_dd(1) = 0; %-4 * radius * pi * pi * cos(ang) / (T * T);
    xi_des_dd(2) = 0;
    xi_des_dd(3) = 0;
end

function [xi_des, xi_des_d, xi_des_dd] = eight_traj(t,T)
    % Parameters
    hover_time = 29;     % Initial hover duration (seconds)
    amplitude = 6;      % Size of the figure-eight (meters)
    height = 2;         % Flight altitude (meters)
    
    % Initial hover phase (first 5 seconds)
    if t < hover_time
        xi_des = [0; 0; height];
        xi_des_d = [0; 0; 0];
        xi_des_dd = [0; 0; 0];
        return;
    end
    
    % Lemniscate phase (after hover)
    t_lem = t - hover_time;                     % Time since hover ended
    period_lem = T - hover_time;                % Period for one lemniscate cycle
    theta = 2 * pi * t_lem / period_lem;        % Phase angle
    
    % Position calculations
    x = amplitude * sin(theta);
    y = amplitude * sin(theta) .* cos(theta);
    z = height;
    
    % Velocity calculations
    w = 2 * pi / period_lem;  % Angular velocity
    x_dot = 0;%amplitude * w * cos(theta);
    y_dot = 0; %amplitude * w * (cos(theta).^2 - sin(theta).^2);
    z_dot = 0;
    
    % Acceleration calculations
    x_ddot = 0;%-amplitude * w^2 * sin(theta);
    y_ddot = 0;%-4 * amplitude * w^2 * sin(theta) .* cos(theta);
    z_ddot = 0;
    
    % Outputs
    xi_des = [x; y; z];
    xi_des_d = [x_dot; y_dot; z_dot];
    xi_des_dd = [x_ddot; y_ddot; z_ddot];
    
    %a = 2;      % Scaling factor for 8 size
    %h = 1;      % Total height gain over period T
    %om = 2 * pi / T;

    % Calculate angular position
    %ang = om * t;
    %s = sin(ang);
    %c = cos(ang);
    %denom = 1 + s.^2; % Denominator term

    % Position (x, y, z)
    %xi_des = zeros(3, 1);
    %xi_des(1) = a * c / denom;          % x: Lemniscate component
    %xi_des(2) = a * s .* c / denom;     % y: Lemniscate component
    %xi_des(3) = h * t / T;              % z: Linear rise

    % First derivatives (velocity)
    %xi_des_d = zeros(3, 1);
    %dxdang = a * (-s.*denom - c.*(2*s.*c)) / (denom.^2);
    %dydang = a * ((c.^2 - s.^2).*denom - (s.*c).*(2*s.*c)) / (denom.^2);

    %xi_des_d(1) = dxdang * om;  % x-velocity
    %xi_des_d(2) = dydang * om;  % y-velocity
    %xi_des_d(3) = h / T;           % z-velocity (constant)

    % Second derivatives (acceleration)
    %xi_des_dd = zeros(3, 1);
    %d2xdang2 = a * ( ...
    %    -c.*denom - 4*s.*(c.*s) - (denom + 2*s.^2).*(-s) + ...
    %    4*s.*(c.^2).*denom - 2*( -s.*denom - 2*c.^2.*s ).*(2*s.*c) ...
    %    ) / (denom.^3);
    
    %d2ydang2 = a * ( ...
    %    -4*s.*c.*denom - 2*(c.^2 - s.^2).*(2*s.*c) - ...
    %    (2*c.^2 - 2*s.^2).*(2*s.*c) - (s.*c).*(-4*s.^2 + 4*c.^2) ...
    %    ) / (denom.^3);
    
    %xi_des_dd(1) = d2xdang2 * om^2; % x-acceleration
    %xi_des_dd(2) = d2ydang2 * om^2; % y-acceleration
    %xi_des_dd(3) = 0;                  % z-acceleration (zero)
end

function [xi_des, xi_des_d, xi_des_dd] = rectangular_traj(t, T)
    % Corner points [x, y] for the rectangle
    corners = [-1, -0.5;  % Bottom-left
                1, -0.5;  % Bottom-right
                1,  0.5;  % Top-right
               -1,  0.5]; % Top-left
    z_height = 1;        % Fixed flight height
    
    % Segment durations (equal time for each edge)
    seg_time = T / 4;
    
    % Normalized time within current lap
    t_norm = mod(t, T);
    
    % Determine current segment and progress ratio
    seg = floor(t_norm / seg_time) + 1;
    s = mod(t_norm, seg_time) / seg_time;  % 0-1 progress in segment
    
    % Wrap around for the last segment
    if seg > 4
        seg = 1;
    end
    
    % Get start and end points for current segment
    start_pt = corners(seg, :);
    
    % Determine next point (wrap around for last segment)
    if seg == 4
        end_pt = corners(1, :);
    else
        end_pt = corners(seg+1, :);
    end
    
    % Calculate position (linear interpolation)
    x = start_pt(1) + s * (end_pt(1) - start_pt(1));
    y = start_pt(2) + s * (end_pt(2) - start_pt(2));
    
    % Calculate velocity (constant between points)
    vx = (end_pt(1) - start_pt(1)) / seg_time;
    vy = (end_pt(2) - start_pt(2)) / seg_time;
    
    % Outputs
    xi_des = [x; y; z_height];
    xi_des_d = [vx; vy; 0];
    xi_des_dd = [0; 0; 0];  % Acceleration = 0 (simplified)
end

% convert position trajectory to angle
function eta_des = pos2ang(xi, xi_des, xi_des_d, xi_des_dd, nu_I, param, PI_var)
    x_des = xi_des(1); x_des_d = xi_des_d(1); x_des_dd = xi_des_dd(1);
    y_des = xi_des(2); y_des_d = xi_des_d(2); y_des_dd = xi_des_dd(2);
    z_des = xi_des(3); z_des_d = xi_des_d(3); z_des_dd = xi_des_dd(3);
    
    x = xi(1); x_d = nu_I(1);
    y = xi(2); y_d = nu_I(2);
    z = xi(3); z_d = nu_I(3);
    
    % Psi_des will be fixed to zero
    eta_des = zeros(3, 1);

    % define the integral error
    int_err_pos = zeros(2,1);

    int_x = int_err_pos(1);
    int_y = int_err_pos(2);

    int_x = int_x + (x_des_d - x_d);
    int_y = int_y + (y_des_d - y_d);

    T_des = param.m * ( z_des_dd + param.g) + param.Kd(3, 3) * z_des_d; 
    eta_des(1) = -(param.m * y_des_dd + param.Kd(2, 2) * y_des_d + PI_var.Kp_xy * (y_des - y) + PI_var.Ki_xy * (int_y)) / T_des;
    eta_des(2) = (param.m * x_des_dd + param.Kd(1, 1) *  x_des_d  + PI_var.Kp_xy * (x_des - x) + PI_var.Ki_xy * (int_x)) / T_des;

    % only for Kp and Ki with x dan y define
    %eta_des(1) = -(param.m * y_des_dd + param.Kd(2, 2) * y_des_d + PI_var.Kp_y * (y_des - y) + PI_var.Ki_y * (int_y)) / T_des;
    %eta_des(2) = (param.m * x_des_dd + param.Kd(1, 1) *  x_des_d  + PI_var.Kp_x * (x_des - x) + PI_var.Ki_x * (int_x)) / T_des;

end