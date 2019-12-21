function q = qmult(qa, qb)
% Quaternion multiplication
% q1, q2 both in form [q0 q1 q2 q3]
% Kok et al. 2017, eq. 3.27

q = zeros(4,1);
q(1) = qa(1) * qb(1) - dot(qa(2:4), qb(2:4));

q(2:4) = qa(1) * qb(2:4) + ...
    qa(2:4) * qb(1) + ...
    cross(qa(2:4), qb(2:4));
end