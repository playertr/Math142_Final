function v = quat2rotvec(q)
    % converts unit quaternion to rotation vector
    v = zeros(3, 1);
    theta = 2 * acos(q(1));
    
    v = theta / sin(theta/2) * q(2:4);
    
end