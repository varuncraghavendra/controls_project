% quadrotor_PID_withEstimator.m
clear; close all; clc;

%% Simulation parameters
dt      = 0.02;             % time step [s]
T_total = 40;               % total simulation time [s]
time    = 0:dt:T_total;

%% Trajectory selection
trajType = 'circle';        % 'circle' or 'square'
R        = 3;               % radius or half‐side length [m]
v_ref    = 0.8;             % reference speed [m/s]

%% State initialization (true plant)
x   = 0;  y   = 0;
vx  = 0;  vy  = 0;

%% Initial estimate states
x_hat  = 0; y_hat  = 0;
vx_hat = 0; vy_hat = 0;

%% PID gains (for controller)
Kp_xy = 2.0;
Ki_xy = 0.05;
Kd_xy = 1.0;

int_err_x  = 0; prev_err_x = 0;
int_err_y  = 0; prev_err_y = 0;

%% Estimator gains (tune these)
L_pos = 2.0;   % gain on position error into position estimate
L_vel = 1.0;   % gain on position error into velocity estimate

%% Storage
X    = zeros(size(time));  Y    = zeros(size(time));
Xd   = zeros(size(time));  Yd   = zeros(size(time));
Xhat = zeros(size(time));  Yhat = zeros(size(time));

for k = 1:length(time)
    t = time(k);
    % Desired trajectory
    switch trajType
        case 'circle'
            Xd(k) = R * cos((v_ref/R)*t);
            Yd(k) = R * sin((v_ref/R)*t);
        case 'square'
            period = 4*(2*R)/v_ref;
            tau    = mod(t,period);
            side_t = (2*R)/v_ref;
            if tau < side_t
                Xd(k) = -R + v_ref*tau;  Yd(k) = -R;
            elseif tau < 2*side_t
                Xd(k) =  R;              Yd(k) = -R + v_ref*(tau-side_t);
            elseif tau < 3*side_t
                Xd(k) =  R - v_ref*(tau-2*side_t); Yd(k) = R;
            else
                Xd(k) = -R;              Yd(k) = R - v_ref*(tau-3*side_t);
            end
        otherwise
            error('Unknown trajectory type');
    end
    
    % Controller error (based on estimated state)
    err_x = Xd(k) - x_hat;
    err_y = Yd(k) - y_hat;
    int_err_x = int_err_x + err_x * dt;
    int_err_y = int_err_y + err_y * dt;
    der_err_x = (err_x - prev_err_x)/dt;
    der_err_y = (err_y - prev_err_y)/dt;
    prev_err_x = err_x; prev_err_y = err_y;
    
    ux = Kp_xy*err_x + Ki_xy*int_err_x + Kd_xy*der_err_x;
    uy = Kp_xy*err_y + Ki_xy*int_err_y + Kd_xy*der_err_y;
    
    % True plant update
    x   = x   + vx  * dt;
    y   = y   + vy  * dt;
    vx  = vx  + ux  * dt;
    vy  = vy  + uy  * dt;
    
    % Measurement (we measure position only)
    meas_x = x + 0.0*randn;
    meas_y = y + 0.0*randn;
    
    % Estimator update: continuous time approx with Euler step
    % dx_hat/dt   = vx_hat + L_pos*(meas_x - x_hat)
    % dvx_hat/dt  = ux       + L_vel*(meas_x - x_hat)
    x_hat_dot  = vx_hat + L_pos*(meas_x - x_hat);
    vx_hat_dot = ux      + L_vel*(meas_x - x_hat);
    y_hat_dot  = vy_hat + L_pos*(meas_y - y_hat);
    vy_hat_dot = uy      + L_vel*(meas_y - y_hat);
    
    x_hat  = x_hat  + x_hat_dot  * dt;
    vx_hat = vx_hat + vx_hat_dot * dt;
    y_hat  = y_hat  + y_hat_dot  * dt;
    vy_hat = vy_hat + vy_hat_dot * dt;
    
    % Store
    X(k)    = x;   Y(k)    = y;
    Xhat(k) = x_hat; Yhat(k) = y_hat;
end

%% Metrics
err_x   = Xd - X;
err_y   = Yd - Y;
final_err_x = err_x(end); final_err_y = err_y(end);
MAE_x = mean(abs(err_x));  MAE_y = mean(abs(err_y));
RMSE_x = sqrt(mean(err_x.^2)); RMSE_y = sqrt(mean(err_y.^2));

fprintf('\n*** Trajectory Tracking + Estimator Metrics ***\n');
fprintf('Trajectory type: %s\n', trajType);
fprintf('Controller PID gains: Kp=%.3f, Ki=%.3f, Kd=%.3f\n', Kp_xy, Ki_xy, Kd_xy);
fprintf('Estimator gains: L_pos=%.3f, L_vel=%.3f\n', L_pos, L_vel);
fprintf('Final error: X = %.4f m, Y = %.4f m\n', final_err_x, final_err_y);
fprintf('MAE: X = %.4f m, Y = %.4f m\n', MAE_x, MAE_y);
fprintf('RMSE: X = %.4f m, Y = %.4f m\n', RMSE_x, RMSE_y);

%% Plot true vs estimated trajectory
figure;
plot(Xd, Yd, 'r--','LineWidth',1.5); hold on;
plot(X,  Y,  'b-','LineWidth',1.5);
plot(Xhat, Yhat, 'g:','LineWidth',1.5);
legend('Desired','True','Estimated');
xlabel('X (m)'); ylabel('Y (m)');
title(sprintf('Tracking with State Estimator — %s', trajType));
axis equal; grid on;
