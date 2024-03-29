function R = quat2rotm( q )
%QUAT2ROTM Convert quaternion to rotation matrix
%   R = QUAT2ROTM(Q) converts a unit quaternion, Q, into an orthonormal
%   rotation matrix, R. The input, Q, is an N-by-4 matrix containing N quaternions. 
%   Each quaternion represents a 3D rotation and is of the form q = [w x y z], 
%   with a scalar number as the first value. Each element of Q must be a real number.
%   The output, R, is an 3-by-3-by-N matrix containing N rotation matrices.
%
%   Example:
%      % Convert a quaternion to rotation matrix
%      q = [0.7071 0.7071 0 0];
%      R = quat2rotm(q)
%
%   See also rotm2quat

% Reshape the quaternions in the depth dimension
q = reshape(q,[4 1 size(q,2)]);

s = q(1,1,:);
x = q(2,1,:);
y = q(3,1,:);
z = q(4,1,:);

R = [   1-2*(y.^2+z.^2)   2*(x.*y-s.*z)   2*(x.*z+s.*y)
    2*(x.*y+s.*z) 1-2*(x.^2+z.^2)   2*(y.*z-s.*x)
    2*(x.*z-s.*y)   2*(y.*z+s.*x) 1-2*(x.^2+y.^2)   ];

end
