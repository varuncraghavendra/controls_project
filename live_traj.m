% quadrotor_PID_vs_LQR_Tuned_Animated.m
% Planar quadrotor trajectory tracking: PID vs Integral-Augmented LQR
% with 3D animation of the LQR quadrotor and moving target.
%
% Model: point-mass planar dynamics
%   x_dot = vx
%   y_dot = vy
%   vx_dot = ux
%   vy_dot = uy
%
% LQR is designed on an augmented error state:
%   [ex; ey; vx; vy; int_ex; int_ey]
%
% Animation style inspired by Simulink 3D animation helpers
% (NewFigure, AnimEulerTar, matrixB2I, etc.)

clear; close all; clc;

%% Simulation parameters
dt      = 0.02;             % time step [s]
T_total = 40;               % total simulation time [s]
time    = 0:dt:T_total;
N       = numel(time);

g       = 9.81;             % gravity (for attitude approximation in anim)

%% Trajectory selection
trajType = 'circle';        % 'circle' or 'square'
R        = 3;               % radius / half-side length [m]
v_ref    = 1.8;             % reference speed [m/s]

%% Initial states (planar position + velocity) for PID and LQR plants
x_PID   = -R;  y_PID   = -R;
vx_PID  = 0;   vy_PID  = 0;

x_LQR   = -R;  y_LQR   = -R;
vx_LQR  = 0;   vy_LQR  = 0;

%% PID controller gains (position level on x,y)
Kp_xy = 2.0;
Ki_xy = 0.05;
Kd_xy = 1.0;

int_err_x_PID  = 0; prev_err_x_PID = 0;
int_err_y_PID  = 0; prev_err_y_PID = 0;

%% LQR: Integral-augmented error-state design
% State vector:
%   x_aug = [ex; ey; vx; vy; int_ex; int_ey]
% Plant model (double integrator in x and y) for error dynamics:
%   ex_dot = vx
%   ey_dot = vy
%   vx_dot = ux
%   vy_dot = uy
%   int_ex_dot = ex
%   int_ey_dot = ey
%
% Continuous-time state-space:
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

% LQR tuning (similar spirit to quadrotor paper: heavier on position + integrals)
Q_aug = diag([50, 50, 5, 5, 20, 20]);   % position, velocity, integral
R_aug = diag([1, 1]);                   % control effort penalty

% Continuous-time infinite-horizon LQR gain
K_LQR_aug = lqr(A_aug, B_aug, Q_aug, R_aug);

%% Storage for plotting / animation
Xd      = zeros(1, N);  Yd      = zeros(1, N);      % reference
X_PID   = zeros(1, N);  Y_PID   = zeros(1, N);      % PID position
X_LQR   = zeros(1, N);  Y_LQR   = zeros(1, N);      % LQR position
VX_LQR  = zeros(1, N);  VY_LQR  = zeros(1, N);      % LQR velocities
UX_LQR  = zeros(1, N);  UY_LQR  = zeros(1, N);      % LQR accelerations

%% Initial error states for LQR integral terms
int_ex = 0;
int_ey = 0;

%% Main simulation loop
for k = 1:N
    t = time(k);

    % -------------------------------
    % Desired trajectory (planar)
    % -------------------------------
    switch trajType
        case 'circle'
            % Constant speed circular motion: x = R cos(omega t), y = R sin(omega t)
            omega = v_ref / R;
            Xd(k) = R * cos(omega * t);
            Yd(k) = R * sin(omega * t);

        case 'square'
            % Axis-aligned square with side 2R, starting from (-R,-R),
            % moving counter-clockwise at constant speed v_ref.
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
    % Error defined as plant - reference (sign absorbed into K)
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
    VX_LQR(k) = vx_LQR;
    VY_LQR(k) = vy_LQR;
    UX_LQR(k) = ux_LQR;
    UY_LQR(k) = uy_LQR;
end

%% Performance metrics
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

%% Static 3D Visualization (as before)
figure;
plot3(Xd,    Yd,    zeros(size(time)), 'k--','LineWidth',1.5); hold on;
plot3(X_PID, Y_PID, 0.2*ones(size(time)), 'b-','LineWidth',1.2);
plot3(X_LQR, Y_LQR,-0.2*ones(size(time)), 'r-','LineWidth',1.2);
xlabel('X (m)'); ylabel('Y (m)'); zlabel('Z (m)');
legend('Desired trajectory','PID actual','LQR actual','Location','best');
title('PID vs LQR Trajectory Tracking (Planar Model)');
grid on; axis equal;
view(45,30);

%% Build data for 3D animation (LQR vs moving target)
% Use planar position for x,y, fix altitude z=0.5 m
XYZ_LQR = [X_LQR.'  , Y_LQR.'  , 0.5*ones(N,1)];
VXYZ_LQR = [VX_LQR.', VY_LQR.', zeros(N,1)];
TAR      = [Xd.'      , Yd.'      , 0.5*ones(N,1)];

% Approximate attitude from accelerations:
% For small angles, ux ≈ g * theta , uy ≈ -g * phi
phi   = -atan2(UY_LQR, g);   % roll
theta =  atan2(UX_LQR, g);   % pitch
psi   = zeros(1, N);         % yaw (keep zero for planar tracking)
EulerAngles = [phi.' theta.' psi.'];

%% 3D Animation (LQR + target)
figure;
xlimit = [min(X_LQR)-1, max(X_LQR)+1];
ylimit = [min(Y_LQR)-1, max(Y_LQR)+1];
zlimit = [0, 2];
width  = 800;
height = 650;
NewFigure(xlimit, ylimit, zlimit, -43, 25, width, height);
pause(1);
AnimEulerTar(time.', XYZ_LQR, EulerAngles, VXYZ_LQR, TAR);

%% ========= Local Functions (animation utilities) =========

function NewFigure(xlimVals, ylimVals, zlimVals, viewx, viewy, w, h)
    ax = gca;
    set(ax, 'XLim', xlimVals, 'YLim', ylimVals, 'ZLim', zlimVals);
    view(viewx, viewy);
    x0 = 50; y0 = 50;
    set(gcf, 'Position', [x0, y0, w, h]);
    hold on;
    grid on;
    xlabel('X (m)');
    ylabel('Y (m)');
    zlabel('Z (m)');
    title('LQR Quadrotor 3D Animation with Moving Target');
end

function AnimEulerTar(t_plot, XYZs, EulerAngles, VXYZs, Tars)
    t_section = 0;
    curve    = animatedline('LineWidth', 1);
    curveTR  = animatedline('LineWidth', 1, 'LineStyle', ':');

    for i = 1:length(t_plot)
        if abs(t_plot(i) - t_section) < 1e-4
            Euler = EulerAngles(i,:);   % [phi, theta, psi]
            XYZ   = XYZs(i,:);          % [x,y,z]
            VXYZ  = VXYZs(i,:);         % [vx,vy,vz]
            TR    = Tars(i,:);          % target position

            O = eye(3);
            T_BtoI = matrixB2I(Euler(1), Euler(2), Euler(3));
            O_I    = T_BtoI * O;        % body axes in inertial frame

            addpoints(curve,   XYZ(1), XYZ(2), XYZ(3));
            addpoints(curveTR, TR(1),  TR(2),  TR(3));
            head = scatter3(TR(1), TR(2), TR(3), 'filled', ...
                            'MarkerFaceColor','k', 'MarkerEdgeColor','k');

            line1 = drawline(XYZ, O_I(:,1), 'b--', 1.5);  % body x-axis
            line2 = drawline(XYZ, O_I(:,2), 'g--', 1.5);  % body y-axis
            line3 = drawline(XYZ, O_I(:,3), 'r--', 1.5);  % body z-axis
            line5 = extendline(XYZ, O_I(:,1), 'b:');      % heading line

            frame1 = drawline(XYZ, 0.5*O_I(:,1)+0.5*O_I(:,2), 'k', 2.0);
            frame2 = drawline(XYZ, 0.5*O_I(:,1)-0.5*O_I(:,2), 'k', 2.0);
            frame3 = drawline(XYZ,-0.5*O_I(:,1)+0.5*O_I(:,2), 'k', 2.0);
            frame4 = drawline(XYZ,-0.5*O_I(:,1)-0.5*O_I(:,2), 'k', 2.0);

            drawnow;
            pause(0.01);

            % Logging / title text
            xlabel(sprintf('t = %.1f s', t_plot(i)));
            vstr = sprintf('Velocity [%.1f, %.1f, %.1f] m/s', ...
                          VXYZ(1), VXYZ(2), VXYZ(3));
            eulStr = sprintf('Euler [%.1f, %.1f, %.1f] deg', ...
                             rad2deg(Euler(1)), rad2deg(Euler(2)), rad2deg(Euler(3)));
            posStr = sprintf('XYZ [%.1f, %.1f, %.1f] m', XYZ(1), XYZ(2), XYZ(3));
            title(sprintf('%s   %s   %s', eulStr, posStr, vstr));

            t_section = t_section + 0.4;   % animation sampling interval

            if i ~= length(EulerAngles)
                delete(line1); delete(line2); delete(line3);
                delete(line5);
                delete(frame1); delete(frame2); delete(frame3); delete(frame4);
                delete(head);
            end
        end
    end
end

function m = matrixB2I(phi, theta, psi)
    % Body-to-Inertial rotation (Z-Y-X convention with negative angles,
    % matching the original helper)
    T_BtoV2   = [ 1      0           0;
                  0  cos(-phi)  sin(-phi);
                  0 -sin(-phi)  cos(-phi)];
    T_V2toV1  = [ cos(-theta)  0 -sin(-theta);
                  0           1  0;
                  sin(-theta)  0  cos(-theta)];
    T_V1toI   = [ cos(-psi)  sin(-psi) 0;
                 -sin(-psi)  cos(-psi) 0;
                  0          0        1];
    m = T_V1toI * T_V2toV1 * T_BtoV2;
end

function line = drawline(p1, p2, color, width)
    if nargin < 4
        width = 1.0;
    end
    pt1 = p1;
    pt2 = pt1 + p2.';
    pts = [pt1; pt2];
    line = plot3(pts(:,1), pts(:,2), pts(:,3), color, 'LineWidth', width);
end

function line = extendline(p1, p2, color)
    pt1 = p1;
    pt2 = pt1 + 20 * p2.';  % extended ray
    pts = [pt1; pt2];
    line = plot3(pts(:,1), pts(:,2), pts(:,3), color, 'LineWidth', 0.5);
end
