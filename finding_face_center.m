function [sortedVertices, face_center] = finding_face_center(currentVertices)
    v1 = currentVertices(2, :) - currentVertices(1, :);
    v2 = currentVertices(3, :) - currentVertices(1, :);
    normal = cross(v1, v2);
    normal = normal / norm(normal); 
    face_center = mean(currentVertices, 1);
    basis_a = currentVertices(1, :) - face_center;
    basis_b = cross(normal, basis_a);
    projectedVertices = zeros(size(currentVertices, 1), 2);
    for j = 1:size(currentVertices, 1)
        projectedVertices(j, 1) = dot(currentVertices(j, :) - face_center, basis_a);
        projectedVertices(j, 2) = dot(currentVertices(j, :) - face_center, basis_b);
    end

    angles = atan2(projectedVertices(:,2), projectedVertices(:,1));
    [~, sortedIndices] = sort(angles, 'descend','MissingPlacement','first');
    sortedVertices = currentVertices(sortedIndices, :);
end