% quadrotor_PID_vs_LQR_Animated_Subplots.m
% Planar quadrotor trajectory tracking: PID vs Integral-Augmented LQR
% with 3 animated subplots and automatic conclusion based on metrics.
%
% Model (point-mass planar):
%   x_dot = vx
%   y_dot = vy
%   vx_dot = ux
%   vy_dot = uy
%
% LQR is designed on augmented error state:
%   [ex; ey; vx; vy; int_ex; int_ey]

clear; close all; clc;

%% Simulation parameters
dt      = 0.02;             % time step [s]
T_total = 40;               % total simulation time [s]
time    = 0:dt:T_total;
N       = numel(time);

%% Trajectory selection
trajType = 'circle';        % 'circle' or 'square'
R        = 3;               % radius / half-side length [m]
v_ref    = 1.8;             % reference speed [m/s]

%% Initial states (planar position + velocity) for PID and LQR plants
x_PID   = -R;  y_PID   = -R;
vx_PID  = 0;   vy_PID  = 0;

x_LQR   = -R;  y_LQR   = -R;
vx_LQR  = 0;   vy_LQR  = 0;

%% PID controller gains (position PID acting on x,y -> accelerations)
Kp_xy = 2.0;
Ki_xy = 0.05;
Kd_xy = 1.0;

int_err_x_PID  = 0; prev_err_x_PID = 0;
int_err_y_PID  = 0; prev_err_y_PID = 0;

%% LQR: Integral-augmented error-state design
% State vector:
%   x_aug = [ex; ey; vx; vy; int_ex; int_ey]
% Continuous-time state-space for error dynamics:
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

% LQR tuning
Q_aug = diag([50, 50, 5, 5, 20, 20]);   % position, velocity, integral
R_aug = diag([1, 1]);                   % control effort

% Continuous-time infinite-horizon LQR gain
K_LQR_aug = lqr(A_aug, B_aug, Q_aug, R_aug);

%% Storage for trajectories
Xd      = zeros(1, N);  Yd      = zeros(1, N);      % reference
X_PID   = zeros(1, N);  Y_PID   = zeros(1, N);      % PID position
X_LQR   = zeros(1, N);  Y_LQR   = zeros(1, N);      % LQR position

%% Initial error states for LQR integral terms
int_ex = 0;
int_ey = 0;

%% ======= Main simulation loop (no plotting here) =======
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
            period = 4 * (2 * R) / v_ref;
            tau    = mod(t, period);
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
    % PID controller (position PID acting on x,y -> accelerations)
    % ==========================================================
    err_x_PID = Xd(k) - x_PID;
    err_y_PID = Yd(k) - y_PID;

    int_err_x_PID = int_err_x_PID + err_x_PID * dt;
    int_err_y_PID = int_err_y_PID + err_y_PID * dt;

    der_err_x_PID = (err_x_PID - prev_err_x_PID) / dt;
    der_err_y_PID = (err_y_PID - prev_err_y_PID) / dt;

    prev_err_x_PID = err_x_PID;
    prev_err_y_PID = err_y_PID;

    ux_PID = Kp_xy * err_x_PID + Ki_xy * int_err_x_PID + Kd_xy * der_err_x_PID;
    uy_PID = Kp_xy * err_y_PID + Ki_xy * int_err_y_PID + Kd_xy * der_err_y_PID;

    % Euler integration for PID plant
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

    % Euler integration for LQR plant
    x_LQR  = x_LQR + vx_LQR * dt;
    y_LQR  = y_LQR + vy_LQR * dt;
    vx_LQR = vx_LQR + ux_LQR * dt;
    vy_LQR = vy_LQR + uy_LQR * dt;

    X_LQR(k) = x_LQR;
    Y_LQR(k) = y_LQR;
end

%% ======= Performance metrics =======
err_x_PID = Xd - X_PID;   err_y_PID = Yd - Y_PID;
MAE_x_PID = mean(abs(err_x_PID));   MAE_y_PID  = mean(abs(err_y_PID));
RMSE_x_PID = sqrt(mean(err_x_PID.^2)); RMSE_y_PID = sqrt(mean(err_y_PID.^2));

err_x_LQR = Xd - X_LQR;   err_y_LQR = Yd - Y_LQR;
MAE_x_LQR = mean(abs(err_x_LQR));   MAE_y_LQR  = mean(abs(err_y_LQR));
RMSE_x_LQR = sqrt(mean(err_x_LQR.^2)); RMSE_y_LQR = sqrt(mean(err_y_LQR.^2));

fprintf('\n--- PID Controller Metrics ---\n');
fprintf('MAE (X,Y):   %.4f, %.4f [m]\n', MAE_x_PID,  MAE_y_PID);
fprintf('RMSE (X,Y):  %.4f, %.4f [m]\n', RMSE_x_PID, RMSE_y_PID);

fprintf('\n--- LQR Controller Metrics (Integral-Augmented) ---\n');
fprintf('MAE (X,Y):   %.4f, %.4f [m]\n', MAE_x_LQR,  MAE_y_LQR);
fprintf('RMSE (X,Y):  %.4f, %.4f [m]\n', RMSE_x_LQR, RMSE_y_LQR);

%% ======= Automatic conclusion based on metrics =======
% Use average RMSE over X,Y as a summary
RMSE_PID_avg = mean([RMSE_x_PID,  RMSE_y_PID]);
RMSE_LQR_avg = mean([RMSE_x_LQR,  RMSE_y_LQR]);
MAE_PID_avg  = mean([MAE_x_PID,   MAE_y_PID]);
MAE_LQR_avg  = mean([MAE_x_LQR,   MAE_y_LQR]);

fprintf('\n=== Controller Comparison Conclusion ===\n');
if (RMSE_LQR_avg < RMSE_PID_avg) && (MAE_LQR_avg < MAE_PID_avg)
    fprintf(['LQR outperforms PID on this run:\n' ...
             '  - Lower average MAE (%.4f vs %.4f m)\n' ...
             '  - Lower average RMSE (%.4f vs %.4f m)\n' ...
             'The integral-augmented LQR reduces steady-state error and ' ...
             'tracks the trajectory more tightly.\n'], ...
             MAE_LQR_avg, MAE_PID_avg, RMSE_LQR_avg, RMSE_PID_avg);
elseif (RMSE_LQR_avg > RMSE_PID_avg) && (MAE_LQR_avg > MAE_PID_avg)
    fprintf(['PID outperforms LQR on this run:\n' ...
             '  - Lower average MAE (%.4f vs %.4f m)\n' ...
             '  - Lower average RMSE (%.4f vs %.4f m)\n' ...
             'With the current tuning, PID tracks slightly better.\n'], ...
             MAE_PID_avg, MAE_LQR_avg, RMSE_PID_avg, RMSE_LQR_avg);
else
    fprintf(['Mixed performance:\n' ...
             '  - PID avg MAE:  %.4f m, LQR avg MAE:  %.4f m\n' ...
             '  - PID avg RMSE: %.4f m, LQR avg RMSE: %.4f m\n' ...
             'One controller may be better on X while the other is better on Y.\n' ...
             'Consider re-tuning Q/R for LQR or Kp/Ki/Kd for PID.\n'], ...
             MAE_PID_avg, MAE_LQR_avg, RMSE_PID_avg, RMSE_LQR_avg);
end

%% ======= Animated plotting with 3 subplots =======
% Slow, step-by-step visualisation with two drones (PID & LQR)
anim_pause = 0.03;   % increase this to slow down further

% Common axis limits
allX = [Xd, X_PID, X_LQR];
allY = [Yd, Y_PID, Y_LQR];
margin = 0.5;
x_min = min(allX) - margin; x_max = max(allX) + margin;
y_min = min(allY) - margin; y_max = max(allY) + margin;

figure('Name','PID vs LQR Trajectory Animation');

% Use tiledlayout for clearer separation
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
title(ax2, 'PID-Controlled Quadrotor');
xlabel(ax2,'X (m)'); ylabel(ax2,'Y (m)');

% Desired as faint grey for reference
plot(ax2, Xd, Yd, 'Color',[0.6 0.6 0.6], 'LineStyle',':');

hTrailPID = animatedline(ax2, 'Color',[0 0.4470 0.7410], 'LineWidth',1.5); % blue line
hDronePID = plot(ax2, X_PID(1), Y_PID(1), 'ro', ...
    'MarkerFaceColor','r', 'MarkerSize',7);                                  % red drone

legend(ax2, {'Desired path','PID path','PID drone'}, 'Location','best');

% ----- Subplot 3: LQR drone -----
ax3 = nexttile(3);
hold(ax3, 'on'); grid(ax3, 'on'); axis(ax3, 'equal');
xlim(ax3, [x_min, x_max]);
ylim(ax3, [y_min, y_max]);
title(ax3, 'LQR-Controlled Quadrotor');
xlabel(ax3,'X (m)'); ylabel(ax3,'Y (m)');

% Desired as faint grey for reference
plot(ax3, Xd, Yd, 'Color',[0.6 0.6 0.6], 'LineStyle',':');

hTrailLQR = animatedline(ax3, 'Color',[0.8500 0.3250 0.0980], 'LineWidth',1.5); % orange line
hDroneLQR = plot(ax3, X_LQR(1), Y_LQR(1), 'mo', ...
    'MarkerFaceColor','m', 'MarkerSize',7);                                    % magenta drone

legend(ax3, {'Desired path','LQR path','LQR drone'}, 'Location','best');

% ----- Animation loop -----
for k = 1:N
    % Desired subplot
    addpoints(hTrailDes, Xd(k), Yd(k));
    set(hDroneDes, 'XData', Xd(k), 'YData', Yd(k));

    % PID subplot
    addpoints(hTrailPID, X_PID(k), Y_PID(k));
    set(hDronePID, 'XData', X_PID(k), 'YData', Y_PID(k));

    % LQR subplot
    addpoints(hTrailLQR, X_LQR(k), Y_LQR(k));
    set(hDroneLQR, 'XData', X_LQR(k), 'YData', Y_LQR(k));

    drawnow;
    pause(anim_pause);    % controls the animation speed
end
