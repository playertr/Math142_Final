% horiz_filter.m
% This file processes data from a rocket flight using a multiplicative EKF.
% Assuming Gaussian process and measurement noise, the filter fuses
% magnetometer, gyroscope, and accelerometer data with a video-identified
% horizon to estimate pose.
% Tim Player, November 2019

%% Add Quaternions to Path 
addpath('./quaternions');

%% Load Data
load('data_for_horiz_filter'); 
camera_time = 1/30 * [1:762];
% These data have multiple sources:
% orientation_filter.m has some of Prof. Spjut's original data processing
% code.
% image_orientation.m has the horizon -> pitch, yaw code.
% process_vids.py has the OpenCV horizon identification.

%% Pose estimation using an EKF with orientation deviation states
% Algorithm 4, pg.41.

% state and covariance
q_nb = zeros(4, n);
P = zeros(9,9,n);
P(:,:,1) = [eye(6,6), zeros(6,3);
            zeros(3,6), 0.0001*eye(3,3)];

% gyroscope covariance matrix (process noise)
sigma_omega = 10^-2 * eye(3);

% accelerometer covariance matrix (process noise)
sigma_ac = [accelXvar, 0, 0; 0, accelYvar, 0; 0, 0, accelZvar];

Q = [5 * eye(3), zeros(3,6);
    zeros(3,3), 5*eye(3), zeros(3,3);
    zeros(3,6), 0.01 * eye(3)];
    
% magnetometer covariance matrix (measurement noise)
R =   eye(3,3);

% measurements 
y_m = [magXrs; magYrs; magZrs]; %magnetometer

y_omega = [gyroXrs; gyroYrs; gyroZrs]; %gyroscope

y_a = [accelXrs; accelYrs; accelZrs];

g_v = [0, 0, -localg]';

m_v = magEarthVec / norm(magEarthVec);

% initialize state (pose with orientation deviations)
x = zeros(9, n); % x, y, z, vx, vy, vz, eta_1, eta_2, eta_3

% Initiate quaternion. This is the point to apply the TRIAD method to
% determine the initial alignment of the rocket on the pad with ENU.
[q_nb(1:4,1), RM0] = TRIAD([0;0;localg], magEarthVec, accelPreVec, magVec0);

for i = 2:n
    
    %%% Time update (prediction)
    
    % Update position (eq. 3.74a)
    x(1:3,i) = x(1:3,i-1) + sp * x(4:6,i-1) + ...
        sp^2/2 * (quat2rotm(q_nb(:,i-1)) * y_a(:,i-1) + g_v); 
    
    % Update velocity (eq. 3.74a)
    x(4:6,i) = x(4:6,i-1) + sp * (quat2rotm(q_nb(:,i-1)) * y_a(:,i-1) + g_v);
    
    % Orientation deviation remains zero
    % x(7:9,i) = [0; 0; 0];
    
    % Update orientation linearity point based on gyro
    q_pred = qmult(q_nb(:,i-1), rotvec2quat( sp/2 * y_omega(:,i-1)));
    
    % Linearized Kinematic matrix per appendix B.4
    F = [eye(3), sp*eye(3), -sp^2/2 * quat2rotm(q_pred) * crossmatrix(y_a(:,i));
        zeros(3,3), eye(3), -sp*quat2rotm(q_pred) * crossmatrix(y_a(:,i));
        zeros(3,3), zeros(3,3), eye(3)];
    
    % Matrix relating effect of process noise on state uncertainty
    G = [eye(6), zeros(6,3);...
        zeros(3,6), sp * quat2rotm(q_pred)];
    
    % Uncertainty in prediction
    P_pred = F * P(:,:,i-1)* F.' + G * Q * G.';
    
    %%% Measurement Update (correction)
    
    % Predict magnetometer measurement
    y_m_pred = quat2rotm(q_pred).' * m_v; % the transpose is the inverse, R_bn
    
    % Measurement dependency matrix for magnetometer
    H = [zeros(3,6), quat2rotm(q_pred).' * crossmatrix(m_v)]; % 3 x 9
    
    S = H * P_pred * H.' + R; % 3 x 3
    
    % Kalman gain
    K = P_pred * H.' * inv(S); 
    
    x(:,i) = x(:, i) + K * (y_m(:,i) - y_m_pred);

    P(:,:,i) = P_pred - K*S*K.';
    
    %%% Relinearize (update orientation linearization point)

    q_nb(:,i) = qmult(rotvec2quat(x(7:9,i)/2), q_pred);
end
%% Verify no roll offset to camera
% cost_of_theta(0, q_nb, pitch, yaw, camera_time, timeArr)
% min_theta = fminsearch(@(x)cost_of_theta(x, q_nb, pitch, yaw, camera_time, timeArr), 0.5) 

%% Plot to determine alignment of camera frame
R_nb = quat2rotm(q_nb);

pred_pitch = reshape(asin(R_nb(3,2,:)), 1, length(R_nb));
pred_yaw = reshape(asin(R_nb(3,1,:)), 1, length(R_nb));

% TODO: make sure the cameras are synchronized.
% Make individual cam1time, cam2time series.

figure(2);clf; hold on;
subplot(2,1,1); hold on;
plot(timeArr, pred_yaw);
plot(camera_time, yaw * 2.2);
xlabel('Time');
ylabel('Yaw (Rad)');
title('Yaw Comparison');
legend('IMU Filtered', 'Camera');
xlim([ 0 22])

subplot(2,1,2); hold on;
plot(timeArr, -1 * pred_pitch);
plot(camera_time, pitch);
xlabel('Time');
ylabel('Pitch (Rad)');
title('Pitch Comparison');
legend('IMU Filtered', 'Camera');
xlim([ 0 22])


%% Resample pitch, yaw.
yawrs = interp1(camera_time, yaw * 2.2, timeArr,'pchip');
pitchrs = interp1(camera_time, pitch * -1, timeArr,'pchip');

%% Pose estimation including camera
% Algorithm 4, pg.41.

% state and covariance
q_nb = zeros(4, n);
P = zeros(9,9,n);
P(:,:,1) =  0.0001 * eye(9);

% gyroscope covariance matrix (process noise)
sigma_omega = 10^-2 * eye(3);

% accelerometer covariance matrix (process noise)
sigma_ac = [accelXvar, 0, 0; 0, accelYvar, 0; 0, 0, accelZvar];

sigma_m = [gyroXvar, 0, 0;0, gyroYvar,0 ; 0, 0, gyroZvar];

Q = [ eye(3), zeros(3,6);
    zeros(3,3), eye(3), zeros(3,3);
    %zeros(3,6), 0.01 * eye(3)];
    zeros(3,6),sigma_omega];
    
% measurement noise)
R =   [ eye(3,3), zeros(3,2); % mag
        zeros(1,3), 0.1, 0; %pitch 
        zeros(1,4),0.1]; %yaw

% measurements 
y_m = [magXrs; magYrs; magZrs]; %magnetometer

y_omega = [gyroXrs; gyroYrs; gyroZrs]; %gyroscope

y_a = [accelXrs; accelYrs; accelZrs];

y_p = pitchrs; % camera pitch

y_y = yawrs; % camera yaw

g_v = [0, 0, -localg]';

m_v = magEarthVec / norm(magEarthVec);

% initialize state (pose with orientation deviations)
x = zeros(9, n); % x, y, z, vx, vy, vz, eta_1, eta_2, eta_3

% Initiate quaternion. This is the point to apply the TRIAD method to
% determine the initial alignment of the rocket on the pad with ENU.
[q_nb(1:4,1), RM0] = TRIAD([0;0;localg], magEarthVec, accelPreVec, magVec0);

for i = 2:23000 %apogee occurs at 23 seconds, or 230000 samples
    
    %%% Time update (prediction)
    
    % Update position (eq. 3.74a)
    x(1:3,i) = x(1:3,i-1) + sp * x(4:6,i-1) + ...
        sp^2/2 * (quat2rotm(q_nb(:,i-1)) * y_a(:,i-1) + g_v); 
    
    % Update velocity (eq. 3.74a)
    x(4:6,i) = x(4:6,i-1) + sp * (quat2rotm(q_nb(:,i-1)) * y_a(:,i-1) + g_v);
    
    % Orientation deviation remains zero
    % x(7:9,i) = [0; 0; 0];
    
    % Update orientation linearity point based on gyro
    q_pred = qmult(q_nb(:,i-1), rotvec2quat( sp/2 * y_omega(:,i-1)));
    
    % Linearized Kinematic matrix per appendix B.4
    F = [eye(3), sp*eye(3), -sp^2/2 * quat2rotm(q_pred) * crossmatrix(y_a(:,i));
        zeros(3,3), eye(3), -sp*quat2rotm(q_pred) * crossmatrix(y_a(:,i));
        zeros(3,3), zeros(3,3), eye(3)];
    
    % Matrix relating effect of process noise on state uncertainty
    G = [sp * eye(6), zeros(6,3);...
        zeros(3,6), sp * quat2rotm(q_pred)];
    
    % Uncertainty in prediction
    P_pred = F * P(:,:,i-1)* F.' + G * Q * G.';
    
    %%% Measurement Update (correction)
    
    % Predict magnetometer measurement
    y_m_pred = quat2rotm(q_pred).' * m_v; % the transpose is the inverse, R_bn
    
    R_nb = quat2rotm(q_pred);
    
    y_p_pred = asin(R_nb(3,2));
    y_y_pred = asin(R_nb(3,1));
 
    % Measurement dependency matrix 
    H = [zeros(3,6), R_nb' * crossmatrix(m_v);
        zeros(1,6), 1/sqrt(1 - R_nb(3,2)^2) * [-R_nb(2,2), -R_nb(1,2), 0]; % pitch
        zeros(1,6), 1/sqrt(1 - R_nb(3,1)^2) * [-R_nb(2,1), -R_nb(1,1), 0]];% yaw
    
    S = H * P_pred * H.' + R; % 3 x 3
    
    % Kalman gain
    K = P_pred * H.' * inv(S); 
    
    y = [y_m(:,i); y_p(i); y_y(i)];
    y_pred = [y_m_pred; y_p_pred; y_y_pred];
    
    x(:,i) = x(:, i) + K * (y - y_pred);

    P(:,:,i) = P_pred - K*S*K.';
    
    %%% Relinearize (update orientation linearization point)

    q_nb(:,i) = qmult(rotvec2quat(x(7:9,i)/2), q_pred);
end

%% more plotting
figure()
clf
%plot3(posGlobVec(1,1:n), posGlobVec(2,1:n), posGlobVec(3,1:n));
axis([-500 1250 -1250 500 0 3000])
hold on
plot3(x(1,:), x(2,:), x(3,:), '.');
title('Trajectory of Rocket');
xlabel('East')
ylabel('North')
zlabel('Up')
scatter3(longDel,latDel, GPSAGL);
legend({'Predicted Ascent Trajectory', 'GPS Descent Trajectory'}, 'Location', 'North')

%% 2D plotting
figure(10); clf;

hold on
subplot(2,1,1); hold on;
plot(x(1,:), x(3,:), '.');
title('Trajectory of Rocket');
xlabel('East')
ylabel('Up')
axis([-250 750 0 3000])
scatter(longDel, GPSAGL);


subplot(2,1,2);hold on;
plot(x(2,:), x(3,:), '.');
xlabel('North')
ylabel('Up')
axis([-750 250 0 3000])
scatter(latDel, GPSAGL);
legend({'Predicted Ascent Trajectory', 'GPS Descent Trajectory'}, 'Location', 'Northeast')

%% more plotting (error)
figure()
clf
%plot3(posGlobVec(1,1:n), posGlobVec(2,1:n), posGlobVec(3,1:n));
axis([-500 1250 -1250 500 0 3000])
hold on
plot3(x(1,:), x(2,:), x(3,:), '.');

for i = 1: 3000 : n
    %a, b, c] = ellipsoid(x(1,i), x(2,i), x(3,i),P(1, 1,i), P(2,2,i), P(3,3,i));
    %surf(a, b, c)
    plot_ellipse(P(1:3,1:3,i), [x(1,i), x(2,i), x(3,i)]);
end
title('Trajectory with Uncertainty Ellipses');
xlabel('East')
ylabel('North')
zlabel('Up')
scatter3(longDel,latDel, GPSAGL);

%% Plot global orientation on sphere
nn = 20;                             % Facets On Sphere
[Xs,Ys,Zs] = sphere(nn);
figure(9); clf;
sc = 0.5;
surf(sc * Xs,sc * Ys, sc * Zs)
hold on

% am = quiver3(0,0,0,... %actual magnetic north
%         m_v(1), m_v(2), m_v(3), 1.1,'m');
    
title('Trace of axes of R^{nb} for t \in [0, 1]')
    
xlabel('East');
ylabel('North');
zlabel('Up');
for i = 400: 10 : 1000
    i
    R_nb = quat2rotm(q_nb(:,i));
    
    r1 = quiver3(0,0,0,...
        R_nb(1,1), R_nb(2,1), R_nb(3,1), 1.1,'b'); % yaw
    r2=quiver3(0,0,0,...
        R_nb(1,2), R_nb(2,2), R_nb(3,2), 1.1,'b'); % ptich
    r3=quiver3(0,0,0,...
        R_nb(1,3), R_nb(2,3), R_nb(3,3), 1.1,'r'); % roll
    
    o1 = cross(y_omega(:,i), [R_nb(1,1), R_nb(2,1), R_nb(3,1)]);
    quiver3( R_nb(1,1), R_nb(2,1), R_nb(3,1),...
        o1(1), o1(2), o1(3), 0.1,'k');
    
    o2 = cross(y_omega(:,i), [R_nb(1,2), R_nb(2,2), R_nb(3,2)]);
    quiver3( R_nb(1,2), R_nb(2,2), R_nb(3,2),...
        o2(1), o2(2), o2(3), 0.1,'k');
    
    o3 = cross(y_omega(:,i), [R_nb(1,3), R_nb(2,3), R_nb(3,3)]);
    quiver3( R_nb(1,3), R_nb(2,3), R_nb(3,3),...
        o3(1), o3(2), o3(3), 0.1,'k');
    
%     pred_mn = R_nb * y_m(:,i); %predicted magnetic north
%     quiver3(0,0,0,...
%         pred_mn(1), pred_mn(2), pred_mn(3), 1.1,'k');
    drawnow;
end
hold off
grid on
axis equal

%% Generate rotation axis plot
figure(10); clf; hold on;
axis([-0.5 1 -0.5 1 -0.5 1])

R_nb = eye(3);

r1 = quiver3(0,0,0,...
    R_nb(1,1), R_nb(2,1), R_nb(3,1), 'linewidth', 5, 'Color', 'r');
r2=quiver3(0,0,0,...
    R_nb(1,2), R_nb(2,2), R_nb(3,2), 'linewidth', 5, 'Color', 'b');
r2=quiver3(0,0,0,...
    R_nb(1,2), R_nb(2,2), R_nb(3,2), 'linewidth', 5, 'Color', 'r');
r3=quiver3(0,0,0,...
    R_nb(1,3), R_nb(2,3), R_nb(3,3), 'linewidth', 5, 'Color', 'r');

R_nb = quat2rotm(qmult([1 0 0 0]', rotvec2quat([0.06;0.06;0.06])));

r1 = quiver3(0,0,0,...
    R_nb(1,1), R_nb(2,1), R_nb(3,1), 'linewidth', 5, 'Color', 'b');
r2=quiver3(0,0,0,...
    R_nb(1,2), R_nb(2,2), R_nb(3,2), 'linewidth', 5, 'Color', 'b');
r3=quiver3(0,0,0,...
    R_nb(1,3), R_nb(2,3), R_nb(3,3), 'linewidth', 5, 'Color', 'b');

title('Navigation and Body Frames');
legend('Navigation', 'Body');
xlabel('X');
ylabel('Y');
zlabel('Z');



