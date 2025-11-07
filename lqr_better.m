% quadrotor_PID_vs_LQR_Tuned.m
% Compare PID vs improved LQR controller for planar quadrotor trajectory tracking.

clear; close all; clc;

%% Simulation parameters
dt      = 0.02;             % time step [s]
T_total = 40;               % total simulation time [s]
time    = 0:dt:T_total;

%% Trajectory selection
trajType = 'circle';        % 'circle' or 'square'
R        = 3;               % radius / half-side length [m]
v_ref    = 0.8;             % reference speed [m/s]

%% Plant (true) initial states for both
x_PID   = -R; y_PID   = -R;
vx_PID  = 0;  vy_PID  = 0;

x_LQR   = -R; y_LQR   = -R;
vx_LQR  = 0;  vy_LQR  = 0;

%% Controller: PID gains
Kp_xy = 2.0; Ki_xy = 0.05; Kd_xy = 1.0;
int_err_x_PID  = 0; prev_err_x_PID = 0;
int_err_y_PID  = 0; prev_err_y_PID = 0;

%% Controller: Improved LQR formulation for plant (x,y,vx,vy)
% Add integral error state to improve steady‐state tracking
% State vector: [ex; ey; vx; vy; int_ex; int_ey]
A_aug = [ 0 0 1 0 0 0;
          0 0 0 1 0 0;
          0 0 0 0 0 0;
          0 0 0 0 0 0;
          1 0 0 0 0 0;
          0 1 0 0 0 0 ];
B_aug = [ 0 0;
          0 0;
          1 0;
          0 1;
          0 0;
          0 0 ];
      
% Tuning Q and R
Q_aug = diag([50,50,5,5,20,20]);     % heavier on position error and integral error
R_aug = diag([1,1]);                 % moderate control effort penalty

K_LQR_aug = lqr(A_aug, B_aug, Q_aug, R_aug);

%% Storage for plotting
Xd = zeros(size(time)); Yd = zeros(size(time));
X_PID = zeros(size(time));  Y_PID = zeros(size(time));
X_LQR = zeros(size(time));  Y_LQR = zeros(size(time));

%% Initial error states for LQR
int_ex = 0; int_ey = 0;

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
    
    %% PID controller step ---
    err_x_PID = Xd(k) - x_PID;
    err_y_PID = Yd(k) - y_PID;
    int_err_x_PID = int_err_x_PID + err_x_PID*dt;
    int_err_y_PID = int_err_y_PID + err_y_PID*dt;
    der_err_x_PID = (err_x_PID - prev_err_x_PID)/dt;
    der_err_y_PID = (err_y_PID - prev_err_y_PID)/dt;
    prev_err_x_PID = err_x_PID;
    prev_err_y_PID = err_y_PID;
    ux_PID = Kp_xy*err_x_PID + Ki_xy*int_err_x_PID + Kd_xy*der_err_x_PID;
    uy_PID = Kp_xy*err_y_PID + Ki_xy*int_err_y_PID + Kd_xy*der_err_y_PID;
    x_PID  = x_PID + vx_PID*dt;
    y_PID  = y_PID + vy_PID*dt;
    vx_PID = vx_PID + ux_PID*dt;
    vy_PID = vy_PID + uy_PID*dt;
    X_PID(k) = x_PID;
    Y_PID(k) = y_PID;
    
    %% LQR controller step ---
    % Compute error states for LQR
    ex = x_LQR - Xd(k);
    ey = y_LQR - Yd(k);
    int_ex = int_ex + ex*dt;
    int_ey = int_ey + ey*dt;
    state_LQR = [ ex; ey; vx_LQR; vy_LQR; int_ex; int_ey ];
    
    u_LQR = -K_LQR_aug * state_LQR;
    ux_LQR = u_LQR(1);
    uy_LQR = u_LQR(2);
    
    x_LQR  = x_LQR + vx_LQR*dt;
    y_LQR  = y_LQR + vy_LQR*dt;
    vx_LQR = vx_LQR + ux_LQR*dt;
    vy_LQR = vy_LQR + uy_LQR*dt;
    
    X_LQR(k) = x_LQR;
    Y_LQR(k) = y_LQR;
end

%% Performance metrics
err_x_PID = Xd - X_PID;   err_y_PID = Yd - Y_PID;
MAE_x_PID = mean(abs(err_x_PID));  MAE_y_PID = mean(abs(err_y_PID));
RMSE_x_PID = sqrt(mean(err_x_PID.^2)); RMSE_y_PID = sqrt(mean(err_y_PID.^2));

err_x_LQR = Xd - X_LQR;   err_y_LQR = Yd - Y_LQR;
MAE_x_LQR = mean(abs(err_x_LQR));  MAE_y_LQR = mean(abs(err_y_LQR));
RMSE_x_LQR = sqrt(mean(err_x_LQR.^2)); RMSE_y_LQR = sqrt(mean(err_y_LQR.^2));

fprintf('\n--- PID Controller Metrics ---\n');
fprintf('MAE (X,Y): %.4f, %.4f [m]\n', MAE_x_PID, MAE_y_PID);
fprintf('RMSE (X,Y): %.4f, %.4f [m]\n', RMSE_x_PID, RMSE_y_PID);

fprintf('\n--- LQR Controller Metrics (improved) ---\n');
fprintf('MAE (X,Y): %.4f, %.4f [m]\n', MAE_x_LQR, MAE_y_LQR);
fprintf('RMSE (X,Y): %.4f, %.4f [m]\n', RMSE_x_LQR, RMSE_y_LQR);

%% 3D Visualization
figure;
plot3(Xd, Yd, zeros(size(time)), 'k--','LineWidth',1.5); hold on;
plot3(X_PID, Y_PID, zeros(size(time))+0.2, 'b-','LineWidth',1.2);
plot3(X_LQR, Y_LQR, zeros(size(time))-0.2, 'r-','LineWidth',1.2);
xlabel('X (m)'); ylabel('Y (m)'); zlabel('Z (m)');
legend('Desired trajectory','PID actual','LQR actual');
title('3D Visualization of PID vs LQR (tuned LQR)');
grid on; axis equal;
view(45,30);
