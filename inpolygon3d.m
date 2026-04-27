function inside = inpolygon3d(point, vertices)
    % Check if a point is inside a 3D polygon using the convex hull method
    k = convhull(vertices);
    inside = false;
    for i = 1:size(k, 1)
        tri = vertices(k(i, :), :);
        if is_point_in_triangle(point, tri)
            inside = true;
            break;
        end
    end
end

function inside = is_point_in_triangle(pt, tri)
    % Check if a point is inside a 3D triangle using barycentric coordinates
    v0 = tri(2, :) - tri(1, :);
    v1 = tri(3, :) - tri(1, :);
    v2 = pt - tri(1, :);

    dot00 = dot(v0, v0);
    dot01 = dot(v0, v1);
    dot02 = dot(v0, v2);
    dot11 = dot(v1, v1);
    dot12 = dot(v1, v2);

    invDenom = 1 / (dot00 * dot11 - dot01 * dot01);
    u = (dot11 * dot02 - dot01 * dot12) * invDenom;
    v = (dot00 * dot12 - dot01 * dot02) * invDenom;

    inside = (u >= 0) && (v >= 0) && (u + v <= 1);
end