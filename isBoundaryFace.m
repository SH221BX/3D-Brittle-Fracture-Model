function tf = isBoundaryFace(V,bx,tol)
x0 = all(abs(V(:,1) - 0)  <= tol);
xb = all(abs(V(:,1) - bx) <= tol);
y0 = all(abs(V(:,2) - 0)  <= tol);
yb = all(abs(V(:,2) - bx) <= tol);
z0 = all(abs(V(:,3) - 0)  <= tol);
zb = all(abs(V(:,3) - bx) <= tol);
tf = x0 || xb || y0 || yb || z0 || zb;
end
