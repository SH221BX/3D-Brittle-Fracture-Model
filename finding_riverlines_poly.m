function segs3 = finding_riverlines_poly(polyV, phi1, Phi, phi2, stepSpacing, ridgeHeight, faceColor)

if nargin < 5 || isempty(stepSpacing), stepSpacing = 0.35; end
if nargin < 6 || isempty(ridgeHeight), ridgeHeight = 0.25; end
if nargin < 7 || isempty(faceColor),   faceColor   = [0.93 0.69 0.13]; end

col1 = [0.62 0.64 0.67];
col2 = [1 0 0];


lightDir = [0.45; -0.35; 0.82];
ambient  = 0.45;
diffuse  = 0.75;

segs3 = {};

planes = [1 0 0; 0 1 0; 0 0 1];

Rz1 = [cosd(phi1) -sind(phi1) 0; sind(phi1) cosd(phi1) 0; 0 0 1];
Rx  = [1 0 0; 0 cosd(Phi) -sind(Phi); 0 sind(Phi) cosd(Phi)];
Rz2 = [cosd(phi2) -sind(phi2) 0; sind(phi2) cosd(phi2) 0; 0 0 1];
R = Rz1*Rx*Rz2;

rotated_planes = (R*planes')';
n1 = rotated_planes(1,:).'; n1 = n1/norm(n1);

P0 = polyV;
if size(P0,1) < 3, return; end
[P0, ~] = finding_face_center(P0);

pb = mean(P0,1).';

nb = polygon_normal(P0);
if norm(nb) < 1e-14, return; end
nb = nb(:)/norm(nb);
if nb(3) < 0, nb = -nb; end

e1 = (P0(2,:).' - P0(1,:).');
e1 = e1 - dot(e1,nb)*nb;
if norm(e1) < 1e-14
    e1 = cross(nb,[1;0;0]);
    if norm(e1) < 1e-14
        e1 = cross(nb,[0;1;0]);
    end
end
e1 = e1/norm(e1);
e2 = cross(nb,e1);
e2 = e2/norm(e2);

P2 = zeros(size(P0,1),2);
for i = 1:size(P0,1)
    r = P0(i,:).' - pb;
    P2(i,1) = dot(r,e1);
    P2(i,2) = dot(r,e2);
end

d1 = cross(nb,n1);
if norm(d1) < 1e-14, d1 = cross(nb,[1;0;0]); end
u1 = [dot(d1,e1); dot(d1,e2)];
if norm(u1) < 1e-14, u1 = [1;0]; end
u1 = u1/norm(u1);

v = [-u1(2); u1(1)];
v = v/norm(v);

sV = P2*v;
smin = min(sV);
smax = max(sV);

if smax - smin < 1e-12, return; end

nBands = max(2, round((smax - smin)/stepSpacing));
sEdges = linspace(smin, smax, nBands+1);
hVals  = linspace(0, ridgeHeight, nBands);

for k = 1:nBands
    sL = sEdges(k);
    sR = sEdges(k+1);

    B = clip_strip_convex_poly(P2, v, sL, sR);
    if isempty(B) || size(B,1) < 3
        continue
    end

    h = hVals(k);
    B3 = lift_band(B, pb, e1, e2, nb, h);

    if mod(k,2)==1
        baseCol = col1;
    else
        baseCol = col2;
    end

    fcTop = shadeColor(baseCol, nb, lightDir, ambient, diffuse);

    patch('Vertices',B3,'Faces',1:size(B3,1), ...
        'FaceColor',fcTop,'FaceAlpha',0.92, ...
        'EdgeColor','none'); hold on

    if k > 1
        Q = line_segment_in_convex_poly(P2, v, sL);
        if ~isempty(Q)
            h0 = hVals(k-1);
            h1 = hVals(k);

            A0 = (pb + Q(1,1)*e1 + Q(1,2)*e2 + h0*nb).';
            B0 = (pb + Q(2,1)*e1 + Q(2,2)*e2 + h0*nb).';
            A1 = (pb + Q(1,1)*e1 + Q(1,2)*e2 + h1*nb).';
            B1 = (pb + Q(2,1)*e1 + Q(2,2)*e2 + h1*nb).';

            riserVerts = [A0; B0; B1; A1];
            nr = cross(riserVerts(2,:) - riserVerts(1,:), riserVerts(4,:) - riserVerts(1,:));
            if norm(nr) > 1e-14
                nr = nr(:) / norm(nr);
            else
                nr = nb;
            end

            fcSide = shadeColor(darkenColor(baseCol, 0.72), nr, lightDir, ambient, diffuse);

            patch('Vertices',riserVerts,'Faces',[1 2 3 4], ...
                'FaceColor',fcSide,'FaceAlpha',0.95, ...
                'EdgeColor','none'); hold on

            plot3([A1(1) B1(1)], [A1(2) B1(2)], [A1(3) B1(3)], ...
                'Color', darkenColor(baseCol,0.25), 'LineWidth', 0.9);
        end
    end
end

end

function n = polygon_normal(V)
n = [0 0 0];
m = size(V,1);
for i = 1:m
    p = V(i,:);
    q = V(1+mod(i,m),:);
    n = n + cross(p,q);
end
if norm(n) < 1e-14
    n = cross(V(2,:)-V(1,:), V(3,:)-V(1,:));
end
end

function B3 = lift_band(B, pb, e1, e2, nb, h)
B3 = zeros(size(B,1),3);
for i = 1:size(B,1)
    X = pb + B(i,1)*e1 + B(i,2)*e2 + h*nb;
    B3(i,:) = X.';
end
end

function Pout = clip_strip_convex_poly(P, v, sL, sR)
P1 = clip_half_plane(P,  v, sL, true);
Pout = clip_half_plane(P1, v, sR, false);
end

function P2 = clip_half_plane(P, v, s0, keepGreater)
if isempty(P), P2 = []; return; end

tol = 1e-12;
Q = [];
m = size(P,1);

for i = 1:m
    A = P(i,:);
    B = P(1+mod(i,m),:);

    fA = A*v - s0;
    fB = B*v - s0;

    if keepGreater
        inA = fA >= -tol;
        inB = fB >= -tol;
    else
        inA = fA <= tol;
        inB = fB <= tol;
    end

    if inA && inB
        Q = [Q; B];
    elseif inA && ~inB
        t = fA/(fA - fB);
        X = A + t*(B - A);
        Q = [Q; X];
    elseif ~inA && inB
        t = fA/(fA - fB);
        X = A + t*(B - A);
        Q = [Q; X; B];
    end
end

if isempty(Q)
    P2 = [];
    return
end

P2 = unique(round(Q,12),'rows','stable');
if size(P2,1) >= 3
    P2 = order_polygon_2d(P2);
end
end

function Q = line_segment_in_convex_poly(P, v, s0)
m = size(P,1);
v = v(:).';
I = [];

for i = 1:m
    A = P(i,:);
    B = P(1+mod(i,m),:);
    fA = A*v.' - s0;
    fB = B*v.' - s0;

    if abs(fA) < 1e-12
        I = [I; A];
    end

    if fA*fB < 0
        t = fA/(fA - fB);
        X = A + t*(B-A);
        I = [I; X];
    end
end

if isempty(I)
    Q = [];
    return
end

I = unique(round(I,12),'rows','stable');
if size(I,1) < 2
    Q = [];
    return
end

u = null(v);
u = u(:,1).';
proj = I*u.';
[~,ord] = sort(proj);
I = I(ord,:);

Q = [I(1,:); I(end,:)];
end

function P = order_polygon_2d(P)
ctr = mean(P,1);
ang = atan2(P(:,2)-ctr(2), P(:,1)-ctr(1));
[~,ord] = sort(ang);
P = P(ord,:);
end

function c = shadeColor(baseColor, n, lightDir, ambient, diffuse)
n = n(:) / norm(n);
lightDir = lightDir(:) / norm(lightDir);
I = ambient + diffuse * max(dot(n, lightDir), 0);
I = max(0, min(1.2, I));
c = baseColor * I;
c = max(0, min(1, c));
end

function c = darkenColor(baseColor, factor)
c = max(0, min(1, factor * baseColor));
end