% quadrotor_PID_vs_LQR_3Laps_Animated_ImprovedPID.m
% Planar quadrotor trajectory tracking: PID vs Integral-Augmented LQR
% - PID is tuned to be "good" (not terrible)
% - LQR is tuned to be better overall
% - Simulation runs for 3 full laps of the reference trajectory
% - Metrics (MAE, RMSE) and a conclusion are printed after each run
% - 3 animated subplots:
%     1. Desired Trajectory
%     2. PID-controlled drone
%     3. LQR-controlled drone

clear; close all; clc;

%% Trajectory & simulation setup
trajType = 'circle';    % 'circle' or 'square'
R        = 3;           % radius (circle) / half-side (square) [m]
v_ref    = 1.8;         % reference speed [m/s]
numLaps  = 3;           % "iterations" = number of full loops

dt       = 0.02;        % time step [s]

% Compute period and total time based on trajectory type and laps
switch trajType
    case 'circle'
        period = 2*pi*R / v_ref;        % time for 1 lap
    case 'square'
        period = 4*(2*R)/v_ref;         % time for 1 lap
    otherwise
        error('Unknown trajectory type');
end

T_total = numLaps * period;            % total simulation time [s]
time    = 0:dt:T_total;
N       = numel(time);

%% Initial states (planar position + velocity) for both controllers
x_PID   = -R;  y_PID   = -R;
vx_PID  = 0;   vy_PID  = 0;

x_LQR   = -R;  y_LQR   = -R;
vx_LQR  = 0;   vy_LQR  = 0;

%% ===== PID controller gains (improved, with anti-windup) =====
% These gains give decent tracking; LQR is still tuned to outperform.
Kp_xy = 2.2;
Ki_xy = 0.08;
Kd_xy = 0.6;

int_err_x_PID  = 0; prev_err_x_PID = 0;
int_err_y_PID  = 0; prev_err_y_PID = 0;

% Simple integral anti-windup limits
int_limit = 5.0;   % [m*s]

%% ===== LQR: Integral-augmented error-state design =====
% State vector:
%   x_aug = [ex; ey; vx; vy; int_ex; int_ey]
% Error dynamics (continuous-time):
%   ex_dot     = vx
%   ey_dot     = vy
%   vx_dot     = ux
%   vy_dot     = uy
%   int_ex_dot = ex
%   int_ey_dot = ey

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

% LQR tuning: strong position + integral weighting,
% moderate penalty on control inputs.
Q_aug = diag([130, 130, 12, 12, 70, 70]);   % strong ex,ey,int_ex,int_ey
R_aug = diag([0.8, 0.8]);                   % moderate input cost

K_LQR_aug = lqr(A_aug, B_aug, Q_aug, R_aug);

%% Storage for trajectories
Xd      = zeros(1, N);  Yd      = zeros(1, N);      % reference
X_PID   = zeros(1, N);  Y_PID   = zeros(1, N);      % PID position
X_LQR   = zeros(1, N);  Y_LQR   = zeros(1, N);      % LQR position

%% Initial integral error states for LQR
int_ex = 0;
int_ey = 0;

%% ======= Main simulation loop (dynamics only; no plotting) =======
for k = 1:N
    t = time(k);

    % -------------------------------
    % Desired trajectory (planar)
    % -------------------------------
    switch trajType
        case 'circle'
            omega = v_ref / R;
            Xd(k) = R * cos(omega * t);
            Yd(k) = R * sin(omega * t);

        case 'square'
            % Axis-aligned square with side 2R, starting from (-R,-R)
            period_sq = 4 * (2 * R) / v_ref;
            tau    = mod(t, period_sq);
            side_t = (2 * R) / v_ref;

            if tau < side_t
                % Bottom edge: (-R,-R) -> ( R,-R)
                Xd(k) = -R + v_ref * tau;
                Yd(k) = -R;
            elseif tau < 2 * side_t
                % Right edge: ( R,-R) -> ( R, R)
                Xd(k) =  R;
                Yd(k) = -R + v_ref * (tau - side_t);
            elseif tau < 3 * side_t
                % Top edge: ( R, R) -> (-R, R)
                Xd(k) =  R - v_ref * (tau - 2 * side_t);
                Yd(k) =  R;
            else
                % Left edge: (-R, R) -> (-R,-R)
                Xd(k) = -R;
                Yd(k) =  R - v_ref * (tau - 3 * side_t);
            end

        otherwise
            error('Unknown trajectory type');
    end

    % ==========================================================
    % Improved PID controller (position PID -> accelerations)
    % ==========================================================
    err_x_PID = Xd(k) - x_PID;
    err_y_PID = Yd(k) - y_PID;

    int_err_x_PID = int_err_x_PID + err_x_PID * dt;
    int_err_y_PID = int_err_y_PID + err_y_PID * dt;

    % Anti-windup clamping
    int_err_x_PID = max(min(int_err_x_PID, int_limit), -int_limit);
    int_err_y_PID = max(min(int_err_y_PID, int_limit), -int_limit);

    der_err_x_PID = (err_x_PID - prev_err_x_PID) / dt;
    der_err_y_PID = (err_y_PID - prev_err_y_PID) / dt;

    prev_err_x_PID = err_x_PID;
    prev_err_y_PID = err_y_PID;

    ux_PID = Kp_xy * err_x_PID + Ki_xy * int_err_x_PID + Kd_xy * der_err_x_PID;
    uy_PID = Kp_xy * err_y_PID + Ki_xy * int_err_y_PID + Kd_xy * der_err_y_PID;

    % Euler integration: PID plant
    x_PID  = x_PID + vx_PID * dt;
    y_PID  = y_PID + vy_PID * dt;
    vx_PID = vx_PID + ux_PID * dt;
    vy_PID = vy_PID + uy_PID * dt;

    X_PID(k) = x_PID;
    Y_PID(k) = y_PID;

    % ==========================================================
    % LQR controller on augmented error state
    % ==========================================================
    ex = x_LQR - Xd(k);
    ey = y_LQR - Yd(k);

    int_ex = int_ex + ex * dt;
    int_ey = int_ey + ey * dt;

    xLQR_aug = [ex; ey; vx_LQR; vy_LQR; int_ex; int_ey];

    u_LQR_vec = -K_LQR_aug * xLQR_aug;  % [ux_LQR; uy_LQR]
    ux_LQR    = u_LQR_vec(1);
    uy_LQR    = u_LQR_vec(2);

    % Euler integration: LQR plant
    x_LQR  = x_LQR + vx_LQR * dt;
    y_LQR  = y_LQR + vy_LQR * dt;
    vx_LQR = vx_LQR + ux_LQR * dt;
    vy_LQR = vy_LQR + uy_LQR * dt;

    X_LQR(k) = x_LQR;
    Y_LQR(k) = y_LQR;
end

%% ======= Performance metrics after 3 laps =======
err_x_PID = Xd - X_PID;   err_y_PID = Yd - Y_PID;
MAE_x_PID = mean(abs(err_x_PID));    MAE_y_PID  = mean(abs(err_y_PID));
RMSE_x_PID = sqrt(mean(err_x_PID.^2)); RMSE_y_PID = sqrt(mean(err_y_PID.^2));

err_x_LQR = Xd - X_LQR;   err_y_LQR = Yd - Y_LQR;
MAE_x_LQR = mean(abs(err_x_LQR));    MAE_y_LQR  = mean(abs(err_y_LQR));
RMSE_x_LQR = sqrt(mean(err_x_LQR.^2)); RMSE_y_LQR = sqrt(mean(err_y_LQR.^2));

fprintf('\n=== Controller Metrics over %d laps ===\n', numLaps);
fprintf('\n--- PID Controller (Improved) ---\n');
fprintf('MAE (X,Y):   %.4f, %.4f [m]\n', MAE_x_PID,  MAE_y_PID);
fprintf('RMSE (X,Y):  %.4f, %.4f [m]\n', RMSE_x_PID, RMSE_y_PID);

fprintf('\n--- LQR Controller (Integral-Augmented, Tuned) ---\n');
fprintf('MAE (X,Y):   %.4f, %.4f [m]\n', MAE_x_LQR,  MAE_y_LQR);
fprintf('RMSE (X,Y):  %.4f, %.4f [m]\n', RMSE_x_LQR, RMSE_y_LQR);

% Summary scalars
RMSE_PID_avg = mean([RMSE_x_PID,  RMSE_y_PID]);
RMSE_LQR_avg = mean([RMSE_x_LQR,  RMSE_y_LQR]);
MAE_PID_avg  = mean([MAE_x_PID,   MAE_y_PID]);
MAE_LQR_avg  = mean([MAE_x_LQR,   MAE_y_LQR]);

fprintf('\n=== Conclusion for this run ===\n');
fprintf('Average PID  MAE = %.4f m, RMSE = %.4f m\n', MAE_PID_avg,  RMSE_PID_avg);
fprintf('Average LQR  MAE = %.4f m, RMSE = %.4f m\n', MAE_LQR_avg, RMSE_LQR_avg);

if (RMSE_LQR_avg < RMSE_PID_avg) && (MAE_LQR_avg < MAE_PID_avg)
    fprintf(['LQR outperforms the improved PID controller:\n' ...
             '  - Lower average MAE and RMSE over all 3 laps.\n' ...
             '  - Strong weights on position and integral error in Q lead to\n' ...
             '    tighter steady-state tracking and less drift along the path.\n']);
elseif (RMSE_LQR_avg <= RMSE_PID_avg) || (MAE_LQR_avg <= MAE_PID_avg)
    fprintf(['LQR is comparable or slightly better than the improved PID.\n' ...
             'You can further separate performance by increasing Q (ex,ey,int)\n' ...
             'or slightly reducing the PID gains.\n']);
else
    fprintf(['Improved PID still outperforms LQR with this tuning.\n' ...
             'Tweak Q_aug and R_aug to make LQR more aggressive, or soften PID.\n']);
end

%% ======= Animated plotting with 3 subplots (2 drones) =======
anim_pause = 0.03;   % increase to slow down animation further

% Axis limits
allX = [Xd, X_PID, X_LQR];
allY = [Yd, Y_PID, Y_LQR];
margin = 0.5;
x_min = min(allX) - margin; x_max = max(allX) + margin;
y_min = min(allY) - margin; y_max = max(allY) + margin;

figure('Name','PID vs LQR Trajectory Animation (3 laps)');

tiledlayout(1,3);

% ----- Subplot 1: Desired trajectory -----
ax1 = nexttile(1);
hold(ax1, 'on'); grid(ax1, 'on'); axis(ax1, 'equal');
xlim(ax1, [x_min, x_max]);
ylim(ax1, [y_min, y_max]);
title(ax1, 'Desired Trajectory');
xlabel(ax1,'X (m)'); ylabel(ax1,'Y (m)');

hTrailDes = animatedline(ax1, 'Color','k', 'LineStyle','--', 'LineWidth',1.5);
hDroneDes = plot(ax1, Xd(1), Yd(1), 'ko', 'MarkerFaceColor','y', 'MarkerSize',7);
legend(ax1, {'Desired path','Desired point'}, 'Location','best');

% ----- Subplot 2: PID drone -----
ax2 = nexttile(2);
hold(ax2, 'on'); grid(ax2, 'on'); axis(ax2, 'equal');
xlim(ax2, [x_min, x_max]);
ylim(ax2, [y_min, y_max]);
title(ax2, 'PID-Controlled Quadrotor (Improved)');
xlabel(ax2,'X (m)'); ylabel(ax2,'Y (m)');

% Desired as faint grey reference
plot(ax2, Xd, Yd, 'Color',[0.7 0.7 0.7], 'LineStyle',':');
hTrailPID = animatedline(ax2, 'Color',[0 0.4470 0.7410], 'LineWidth',1.5);     % blue path
hDronePID = plot(ax2, X_PID(1), Y_PID(1), 'ro', ...
    'MarkerFaceColor','r', 'MarkerSize',7);                                     % red drone
legend(ax2, {'Desired path','PID path','PID drone'}, 'Location','best');

% ----- Subplot 3: LQR drone -----
ax3 = nexttile(3);
hold(ax3, 'on'); grid(ax3, 'on'); axis(ax3, 'equal');
xlim(ax3, [x_min, x_max]);
ylim(ax3, [y_min, y_max]);
title(ax3, 'LQR-Controlled Quadrotor');
xlabel(ax3,'X (m)'); ylabel(ax3,'Y (m)');

% Desired as faint grey reference
plot(ax3, Xd, Yd, 'Color',[0.7 0.7 0.7], 'LineStyle',':');
hTrailLQR = animatedline(ax3, 'Color',[0.8500 0.3250 0.0980], 'LineWidth',1.5); % orange path
hDroneLQR = plot(ax3, X_LQR(1), Y_LQR(1), 'mo', ...
    'MarkerFaceColor','m', 'MarkerSize',7);                                     % magenta drone
legend(ax3, {'Desired path','LQR path','LQR drone'}, 'Location','best');

% ----- Animation loop -----
for k = 1:N
    % Desired plot
    addpoints(hTrailDes, Xd(k), Yd(k));
    set(hDroneDes, 'XData', Xd(k), 'YData', Yd(k));

    % PID plot
    addpoints(hTrailPID, X_PID(k), Y_PID(k));
    set(hDronePID, 'XData', X_PID(k), 'YData', Y_PID(k));

    % LQR plot
    addpoints(hTrailLQR, X_LQR(k), Y_LQR(k));
    set(hDroneLQR, 'XData', X_LQR(k), 'YData', Y_LQR(k));

    drawnow;
    pause(anim_pause);
end
