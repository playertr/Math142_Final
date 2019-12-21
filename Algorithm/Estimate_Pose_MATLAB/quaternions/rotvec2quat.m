function q = rotvec2quat(v)
    v = reshape(v, 3, 1);
    %converts a rotation vector to a quaternion as in Eq. 3.36a
    theta = norm(v);
    q = zeros(4,1);
    q(1) = cos(theta);
    q(2:4) = sin(theta) * v / norm(v);
end