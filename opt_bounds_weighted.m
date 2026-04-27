function [Theta1_max, Theta2_max, k1, k2] = opt_bounds_weighted(planes1, planes2, R1, R2, nStarts)

if nargin < 5, nStarts = 20000; end

P  = [planes1 planes2];                          % 3xm
m1 = size(planes1,2);
m2 = size(planes2,2);

Rvec = [(R1)*ones(m1,1); (R2)*ones(m2,1)];           % mx1

% Score s = |k·p|/R, large s => small K
score = @(k) abs(P.'*(k(:)/norm(k))) ./ Rvec;    % mx1

theta1_of_k = @(k) acosd( max(score(k)) );       % since best K = 1/max(score)
theta2_of_k = @(k) acosd( second_largest(score(k)) );

% We want maxima of theta1, theta2:
% maximize theta1 <=> minimize max(score)
% maximize theta2 <=> minimize second_largest(score)
obj1 = @(k) max(score(k));                       % minimize
obj2 = @(k) second_largest(score(k));            % minimize

K = fibonacci_sphere(nStarts);                   % 3xnStarts deterministic seeds
opts = optimset('Display','off','TolX',1e-10,'TolFun',1e-12,'MaxIter',5000);

best1 = inf; best2 = inf;
k1 = K(:,1); k2 = K(:,1);

for s = 1:nStarts
    k0 = K(:,s);

    x1 = fminsearch(obj1, k0, opts);
    x1 = x1(:); x1 = x1/norm(x1);
    v1 = obj1(x1);
    if v1 < best1, best1 = v1; k1 = x1; end

    x2 = fminsearch(obj2, k0, opts);
    x2 = x2(:); x2 = x2/norm(x2);
    v2 = obj2(x2);
    if v2 < best2, best2 = v2; k2 = x2; end
end

Theta1_max = acosd(best1);
Theta2_max = acosd(best2);

end

function s2 = second_largest(v)
v = sort(v,'descend');
s2 = v(2);
end

function K = fibonacci_sphere(n)
i = (0:n-1).';
phi = (1 + sqrt(5))/2;
theta = 2*pi*i/phi;
z = 1 - 2*(i + 0.5)/n;
r = sqrt(max(0,1 - z.^2));
x = r.*cos(theta);
y = r.*sin(theta);
K = [x.'; y.'; z.'];
end
