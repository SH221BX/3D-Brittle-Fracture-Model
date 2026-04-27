function [theta, phi, SIF] = finding_gb_plane_angles(vertices, Gamma)
v1 = vertices(2,:) - vertices(1,:);
v2 = vertices(3,:) - vertices(1,:);
n  = cross(v1, v2);
nn = norm(n);
if nn == 0
    theta = NaN; phi = NaN; SIF = NaN; return
end
n  = n / nn;

theta = atan2d(n(2), n(1));            % azimuth of normal
phi   = acosd(abs(n(3)));              % inclination to XY plane in [0,90]

c = max(1e-12, abs(cosd(phi)));        % guard against divide-by-zero
SIF = Gamma / c;
end
