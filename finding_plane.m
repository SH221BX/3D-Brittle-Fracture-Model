function result = finding_plane(startCell, voronoiCells, Pivot, cell_theta_x, cell_theta_y, cell_theta_z, Gamma)

tolerance = 1e-6;  result = [];
targetCell = voronoiCells(startCell);
cellCenter = targetCell.center;

[~, ~, rotated_planes, K_IC] = finding_basis(cell_theta_x(startCell), cell_theta_y(startCell), cell_theta_z(startCell), Gamma);
[~, sortIdx] = sort(K_IC, 'ascend');
rotated_planes = rotated_planes(sortIdx, :);
K_IC = K_IC(sortIdx, :);
normal = rotated_planes(1,:) / norm(rotated_planes(1,:));


pts = [];
for f = 1:length(targetCell.faces)
    faceVerts = targetCell.faces{f};
    nVerts = size(faceVerts, 1);
    for j = 1:nVerts
        v1 = faceVerts(j,:);
        v2 = faceVerts(mod(j, nVerts) + 1,:);
        l = v2 - v1;
        denom = dot(normal, l);
        if abs(denom) > tolerance
            t_val = dot(normal, (Pivot - v1)) / denom;
            if t_val >= 0 && t_val <= 1
                pts = [pts; v1 + t_val * l];
            end
        end
    end
end
pts = unique(round(pts,4), 'rows');

if size(pts,1) >= 3
    B = null(normal);
    pr = pts * B;
    cent2D = mean(pr,1);
    ang = atan2(pr(:,2)-cent2D(2), pr(:,1)-cent2D(1));
    [~, order] = sort(ang);
    sortedPts = pts(order,:);
    % fill3(sortedPts(:,1), sortedPts(:,2), sortedPts(:,3), 'c', 'FaceAlpha',0.3, 'EdgeColor','k');
else
    error('Cell %d does not form a valid intersection polygon.', startCell);
end

result.cell = startCell;
result.sortedPts = sortedPts;
result.K_IC = K_IC(1);

end
