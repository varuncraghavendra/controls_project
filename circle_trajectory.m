% quadrotor_PID_metrics.m
% Simulate planar x-y tracking of square or circle with PID
% and compute performance metrics instead of error plots.

clear; close all; clc;

%% Simulation parameters
dt      = 0.02;            % time step (s)
T_total = 40;              % total simulation time (s)
time    = 0:dt:T_total;

%% Trajectory selection
trajType = 'circle';       % 'circle' or 'square'
R        = 3;              % radius or half‐side length (m)
v_ref    = 0.8;            % reference speed (m/s)

%% State initialization
if strcmp(trajType, 'square')
    x  = -R; y  = -R;      % start at one corner
else
    x  =  R; y  =  0;      % start point for circle (just for reference)
end
vx = 0; vy = 0;

%% PID gains — modify these
Kp_xy = 1.8;
Ki_xy = 1.05;
Kd_xy = 1.0;

int_err_x  = 0; int_err_y  = 0;
prev_err_x = 0; prev_err_y = 0;

%% Storage for tracking
X   = zeros(size(time));
Y   = zeros(size(time));
Xd  = zeros(size(time));
Yd  = zeros(size(time));

for k = 1:length(time)
    t = time(k);
    % Desired trajectory
    switch trajType
        case 'circle'
            Xd(k) = R * cos( (v_ref / R) * t );
            Yd(k) = R * sin( (v_ref / R) * t );
        case 'square'
            period = 4 * (2*R) / v_ref;
            tau    = mod(t, period);
            side_t = (2*R)/v_ref;
            if tau < side_t
                Xd(k) = -R + v_ref * tau;  Yd(k) = -R;
            elseif tau < 2*side_t
                Xd(k) =  R;                Yd(k) = -R + v_ref*(tau-side_t);
            elseif tau < 3*side_t
                Xd(k) =  R - v_ref*(tau-2*side_t);  Yd(k) = R;
            else
                Xd(k) = -R;                Yd(k) = R - v_ref*(tau-3*side_t);
            end
        otherwise
            error('Unknown trajectory type');
    end
    
    % Compute errors
    err_x     = Xd(k) - x;
    err_y     = Yd(k) - y;
    int_err_x = int_err_x + err_x*dt;
    int_err_y = int_err_y + err_y*dt;
    der_err_x = (err_x - prev_err_x) / dt;
    der_err_y = (err_y - prev_err_y) / dt;
    prev_err_x = err_x;
    prev_err_y = err_y;
    
    % PID control
    ux = Kp_xy*err_x + Ki_xy*int_err_x + Kd_xy*der_err_x;
    uy = Kp_xy*err_y + Ki_xy*int_err_y + Kd_xy*der_err_y;
    
    % Dynamics (simple double‐integrator model)
    x  = x + vx*dt;
    y  = y + vy*dt;
    vx = vx + ux*dt;
    vy = vy + uy*dt;
    
    % Store
    X(k) = x;
    Y(k) = y;
end

%% Compute performance metrics
errors_x = Xd - X;
errors_y = Yd - Y;

final_err_x = errors_x(end);
final_err_y = errors_y(end);

MAE_x = mean(abs(errors_x));
MAE_y = mean(abs(errors_y));

RMSE_x = sqrt(mean(errors_x.^2));
RMSE_y = sqrt(mean(errors_y.^2));

fprintf('\n*** PID Trajectory Tracking Metrics ***\n');
fprintf('Trajectory type: %s\n', trajType);
fprintf('PID gains: Kp=%.3f, Ki=%.3f, Kd=%.3f\n', Kp_xy, Ki_xy, Kd_xy);
fprintf('Final error: X = %.4f m, Y = %.4f m\n', final_err_x, final_err_y);
fprintf('Mean Absolute Error (MAE):   X = %.4f m, Y = %.4f m\n', MAE_x, MAE_y);
fprintf('Root Mean Square Error (RMSE): X = %.4f m, Y = %.4f m\n', RMSE_x, RMSE_y);

%% Plot trajectories
figure;
plot(Xd, Yd, 'r--','LineWidth',1.5); hold on;
plot(X,  Y,  'b-','LineWidth',1.5);
legend('Desired','Actual');
xlabel('X (m)'); ylabel('Y (m)');
title(sprintf('Trajectory Tracking — %s', trajType));
axis equal;
grid on;
