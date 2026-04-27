function pts = intersect_plane_with_faces(faces,n,p0,tol)
pts = [];
for i = 1:numel(faces)
    V = faces{i};
    m = size(V,1);
    for j = 1:m
        p1 = V(j,:);
        p2 = V(mod(j,m)+1,:);
        v = p2 - p1;
        den = dot(n,v);

        if abs(den) < tol
            if abs(dot(n,p1-p0)) < tol
                pts = [pts; p1; p2];
            end
            continue
        end

        t = dot(n,p0-p1)/den;
        if t >= -tol && t <= 1+tol
            x = p1 + t*v;
            if all(x >= min(p1,p2)-tol & x <= max(p1,p2)+tol)
                pts = [pts; x];
            end
        end
    end
end
pts = unique(round(pts,6),'rows','stable');
end