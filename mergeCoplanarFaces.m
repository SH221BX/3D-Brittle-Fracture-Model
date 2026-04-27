function [newV, newF] = mergeCoplanarFaces(V, F, tolerance)
    newV = V;  % Copying original vertices, assuming they remain unchanged.
    newF = {};  % Using a cell array to accommodate faces with varying numbers of vertices.

    normals = zeros(size(F, 1), 3);
    for i = 1:size(F, 1)
        normal = cross(V(F(i, 2), :) - V(F(i, 1), :), V(F(i, 3), :) - V(F(i, 1), :));
        normals(i, :) = normal / norm(normal);
    end

    % Identify and group faces with similar normals.
    groupedFaces = {};  % A cell array to hold arrays of face indices.
    usedFaces = false(size(F, 1), 1);

    for i = 1:size(F, 1)
        if usedFaces(i)
            continue;  % Skip faces that have been grouped already.
        end

        % Start a new group with the current face.
        currentGroup = i;
        for j = i+1:size(F, 1)
            if acos(dot(normals(i, :), normals(j, :))) < tolerance
                currentGroup = [currentGroup, j];
                usedFaces(j) = true;
            end
        end

        groupedFaces{end+1} = currentGroup;  % Add the new group to the list of groups.
    end

    % Merge the faces in each group.
    for i = 1:numel(groupedFaces)
        group = groupedFaces{i};
        uniqueVertices = unique(F(group, :),'stable');

        % The new face is simply a list of the unique vertices in the group.
        newF{end+1} = uniqueVertices;  % Add this list to the cell array of faces.
    end
end
