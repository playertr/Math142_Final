function q = eang2quat(eta, alpha)
% Euler angle to unit quaternion
% eta: axis of rotation (unit vector [x y z])
% alpha: rotation in radians
% Kok et al. 2017 equation 3.30

q = zeros(1,4);
q(1) = cos(alpha/2);
q(2:4) = - sin(alpha/2) * eta;
end