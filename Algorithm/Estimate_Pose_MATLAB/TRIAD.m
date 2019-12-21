function [q,RM] = TRIAD(GG, MVG, GL, MVL)
%TRIAD Calc Initial Alignment on pad
%   This function uses the TRIAD method
%   (https://en.wikipedia.org/wiki/Triad_method) to calculate the initial
%   rotation quaternion, q, and the initial rotation matrix, RM, given the
%   gravity column vector in both global, GG, and local, GL, frames, and the
%   magnetic column vector in both global, MVG, and local, ML, frames.

S = GG/norm(GG);
s = GL/norm(GL);
M = cross(GG,MVG)/norm(cross(GG,MVG));
m = cross(GL,MVL)/norm(cross(GL,MVL));
RM = [S M cross(S,M)]*transpose([s m cross(s, m)]);
%q = rotToQuat(RM);
q = rotm2quat(RM);
end

