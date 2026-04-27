clear all; clc;
format short; warning off;
t_total = [];

if ~exist('M2_ArchiveEquR1','var') || ~isstruct(M2_ArchiveEquR1)
    M2_ArchiveEquR1 = struct('rindex',{},'SimStore',{});
end

doPlot = false;

for rindex = 3:1:3

    t_total = tic;
    disp(rindex);
    rng(rindex, "v4");

    bx        = 100;
    Gamma     = 2;    Gamma_G = 1;
    
    
    numPoints = 100;
    tol       = 1e-6;

    BX = [0 0 0; bx 0 0; 0 bx 0; bx bx 0;
          0 0 bx; bx 0 bx; 0 bx bx; bx bx bx];

    seeds = rand(numPoints-1, 3) .* bx;

    t_voro = tic;
    [V1, C1] = voronoi3d_cuboid(seeds, BX);

    for iteration = 1:10
        newSeeds = zeros(size(seeds));
        for i = 1:length(C1)
            if ~isempty(C1{i})
                newSeeds(i,:) = mean(V1(C1{i},:), 1);
            end
        end
        seeds    = newSeeds;
        [V1, C1] = voronoi3d_cuboid(seeds, BX);
    end
    t_voronoi_sec = toc(t_voro);

    numCells     = length(C1);
    cell_theta_x = -90 + 180 * rand(1, numCells);
    cell_theta_y = acosd(2 * rand(1, numCells) - 1) - 90;
    cell_theta_z = -90 + 180 * rand(1, numCells);

    if ~exist('SimStore','var') || ~isstruct(SimStore)
        SimStore = struct('rng_index',{},'Num_Seed',{},'z',{},'PathLog',{}, ...
            'Kmax',{},'K_T1',{},'K_T2',{},'K_GB',{},'srcMax',{},'idxMax',{}, ...
            'I_Cells',{},'M_time',{},'T_time',{},'RowData',{});
    end

    for zz = 50:1:50
        zSlice = zz;

        voronoiCells      = struct();
        intersectingCells = [];

        for k = 1:length(C1)
            if isempty(C1{k}), continue; end

            Vk = V1(C1{k}, :);
            Fk = convhull(Vk);
            [Vk, Fk] = mergeCoplanarFaces(Vk, Fk, 1e-3);

            cellCenter  = mean(Vk, 1);
            totalVolume = 0;

            voronoiCells(k).cellNumber  = k;
            voronoiCells(k).Vertices    = Vk;
            voronoiCells(k).FaceIdx     = Fk';
            voronoiCells(k).center      = cellCenter;
            voronoiCells(k).faces       = cell(length(Fk), 1);
            voronoiCells(k).faceCenters = zeros(length(Fk), 3);
            voronoiCells(k).faceNormals = zeros(length(Fk), 3);
            voronoiCells(k).faceKeys    = cell(length(Fk), 1);

            cellIntersects = false;

            for i = 1:length(Fk)
                currentFace     = Fk{i};
                currentVertices = Vk(currentFace, :);
                [sortedVertices, ~] = finding_face_center(currentVertices);

                voronoiCells(k).faces{i}         = sortedVertices;
                voronoiCells(k).faceCenters(i,:) = mean(sortedVertices, 1);

                nrm = cross(sortedVertices(2,:) - sortedVertices(1,:), ...
                            sortedVertices(3,:) - sortedVertices(1,:));
                nrm = nrm ./ norm(nrm);
                voronoiCells(k).faceNormals(i,:) = nrm;
                voronoiCells(k).faceKeys{i} = mat2str(sortrows(round(sortedVertices, 6)));

                for j = 1:size(sortedVertices,1)-2
                    v1 = sortedVertices(j,:);
                    v2 = sortedVertices(j+1,:);
                    v3 = sortedVertices(j+2,:);
                    totalVolume = totalVolume + ...
                        abs(dot((v1-cellCenter), cross((v2-cellCenter),(v3-cellCenter)))) / 6;
                end

                m = size(currentFace, 1);
                for j = 1:m
                    v1 = Vk(currentFace(j), :);
                    v2 = Vk(currentFace(mod(j,m)+1), :);
                    if (v1(3)-zSlice) * (v2(3)-zSlice) < 0
                        cellIntersects = true;
                        break
                    end
                end
            end

            voronoiCells(k).volume      = totalVolume;
            voronoiCells(k).scaleFactor = 0.25 * min(sqrt(sum((Vk-cellCenter).^2, 2)));
            voronoiCells(k).color       = rand(1,3);

            if cellIntersects
                intersectingCells(end+1) = k;
                patchColor =  voronoiCells(k).color ; patchAlpha = 0.0; patchedge = 0.25;
                patch('Vertices', sortedVertices, 'Faces', patchFaces, 'FaceColor', patchColor, 'FaceAlpha', patchAlpha, 'EdgeAlpha', patchedge, 'EdgeColor', 'k', 'LineStyle', '--','LineWidth',0.15);
                hold on;
            end
        end
        intersectingCells = unique(intersectingCells, 'stable');

        uniqueFaceMap = containers.Map('KeyType','char','ValueType','any');

        for i = 1:length(voronoiCells)
            if isempty(voronoiCells(i).faces), continue; end
            for fi = 1:length(voronoiCells(i).faces)
                key = voronoiCells(i).faceKeys{fi};
                if isKey(uniqueFaceMap, key)
                    uniqueFaceMap(key) = [uniqueFaceMap(key), voronoiCells(i).cellNumber];
                else
                    uniqueFaceMap(key) = voronoiCells(i).cellNumber;
                end
            end
        end

        sharedFaces = {};
        keys_sf = uniqueFaceMap.keys;
        for k = 1:length(keys_sf)
            cells_sf = unique(uniqueFaceMap(keys_sf{k}), 'sorted');
            sharedFaces{end+1} = struct('cells', cells_sf, 'vertices', eval(keys_sf{k}));
        end

        GB = {};
        for k = 1:numel(sharedFaces)
            cells_k = sharedFaces{k}.cells;
            if numel(cells_k) < 2,                         continue; end
            if ~all(ismember(cells_k, intersectingCells)), continue; end

            Vf = sharedFaces{k}.vertices;
            [Vord, fc] = finding_face_center(Vf);
            A = 0;
            for j = 2:size(Vord,1)-1
                A = A + 0.5*norm(cross(Vord(j,:)-Vord(1,:), Vord(j+1,:)-Vord(1,:)));
            end
            nrm = cross(Vord(2,:)-Vord(1,:), Vord(3,:)-Vord(1,:));
            nrm = nrm ./ norm(nrm);
            GB{end+1} = struct('cells',cells_k,'vertices',Vord,'centroid',fc,'normal',nrm,'area',A);
        end
        GB = GB';

        gbIndex = containers.Map('KeyType','char','ValueType','int32');
        gbKeys  = containers.Map('KeyType','char','ValueType','logical');
        for g = 1:numel(GB)
            key          = mat2str(sortrows(round(GB{g}.vertices, 6)));
            gbIndex(key) = g;
            gbKeys(key)  = true;
        end

        In_store     = struct('cell',{},'polygon',{},'walls',{},'neighbors',{},'GB_Idx',{});
        plane_normal = [0 0 1];

        for jj = 1:length(intersectingCells)
            targetCell  = voronoiCells(intersectingCells(jj));
            plane_point = targetCell.center;
            plane_point(3) = zSlice;

            intersectionPoints = [];
            for i = 1:length(targetCell.faces)
                faceVertices = targetCell.faces{i};
                m = size(faceVertices, 1);
                for j = 1:m
                    p1e = faceVertices(j,:);
                    p2e = faceVertices(mod(j,m)+1,:);
                    ls  = p2e - p1e;
                    den = dot(plane_normal, ls);
                    if abs(den) > 1e-3
                        t = dot(plane_normal, plane_point-p1e) / den;
                        if t >= 0 && t <= 1
                            pI = p1e + t*ls;
                            if all(pI >= min(p1e,p2e) & pI <= max(p1e,p2e))
                                intersectionPoints = [intersectionPoints; pI];
                            end
                        end
                    end
                end
            end

            intersectionPoints        = unique(round(intersectionPoints,4),'rows','stable');
            [sorted_intersections, ~] = finding_face_center(intersectionPoints);

            touchedFacesIdx = [];
            tol_touch = 1e-3;
            for i = 1:numel(targetCell.faces)
                Vface = targetCell.faces{i};
                m = size(Vface, 1);
                hits = 0;
                for q = 1:size(sorted_intersections,1)
                    p = sorted_intersections(q,:);
                    for j = 1:m
                        a = Vface(j,:); b = Vface(mod(j,m)+1,:);
                        v = b - a;    w = p - a;
                        if norm(v) < tol_touch, continue; end
                        tpar  = dot(w,v) / dot(v,v);
                        dperp = norm(cross(w,v)) / norm(v);
                        if tpar >= -tol_touch && tpar <= 1+tol_touch && dperp < tol_touch && ...
                                all(p >= min(a,b)-tol_touch & p <= max(a,b)+tol_touch)
                            hits = hits + 1;
                            if hits >= 2
                                touchedFacesIdx(end+1) = i;
                                break
                            end
                        end
                    end
                end
            end

            touchedFacesIdx = unique(touchedFacesIdx, 'stable');
            GB_walls = {};  GB_idx = [];  neighborCells = [];

            for ii = 1:numel(touchedFacesIdx)
                key = targetCell.faceKeys{touchedFacesIdx(ii)};
                if isKey(gbIndex, key)
                    gidx            = gbIndex(key);
                    GB_idx(end+1)   = gidx;
                    GB_walls{end+1} = GB{gidx};
                    cs              = GB{gidx}.cells;
                    neighborCells   = [neighborCells; setdiff(cs, targetCell.cellNumber)];
                end
            end

            [GB_idx, uix] = unique(GB_idx, 'stable');
            GB_walls      = GB_walls(uix);
            neighborCells = unique(neighborCells, 'stable');

            In_store(end+1) = struct('cell',targetCell.cellNumber, ...
                'polygon',sorted_intersections,'walls',{GB_walls}, ...
                'neighbors',neighborCells,'GB_Idx',GB_idx);
        end

        ID         = [];
        K_selected = [];
        PathLog    = struct('cell',{},'step',{},'mode',{},'K_sel',{}, ...
            'K_T1',{},'K_T2',{},'K_GB1',{},'K_GB2',{},'poly',{},'faces',{});

        for w = 1:length(intersectingCells)
            id         = intersectingCells(w);
            targetCell = voronoiCells(id);
            cellID     = targetCell.cellNumber;

            T2_plane = In_store(w).polygon;
            if isempty(T2_plane) || size(T2_plane,1) < 3, continue; end

            [~, ~, ~, K_eff] = finding_basis(cell_theta_x(id), cell_theta_y(id), cell_theta_z(id), Gamma);
            [~, sort_order]  = sort(K_eff, 'ascend');
            T2.K = K_eff(sort_order(2));
            T1.K = K_eff(sort_order(1));

            cellCenter    = mean(T2_plane, 1);
            cellCenter(3) = zSlice;

            T1_plane = finding_plane(cellID, voronoiCells, cellCenter, ...
                                     cell_theta_x, cell_theta_y, cell_theta_z, Gamma);

            n1 = cross(T1_plane.sortedPts(2,:)-T1_plane.sortedPts(1,:), ...
                       T1_plane.sortedPts(3,:)-T1_plane.sortedPts(1,:));
            n1 = n1 ./ norm(n1);
            p1 = mean(T1_plane.sortedPts, 1);

            n2 = cross(T2_plane(2,:)-T2_plane(1,:), T2_plane(3,:)-T2_plane(1,:));
            n2 = n2 ./ norm(n2);
            p2 = mean(T2_plane, 1);

            if n1(3) < 0,           n1 = -n1; end
            if dot(n2,[0 0 1]) < 0, n2 = -n2; end

            epsn = 1e-6;
            p2a  = p2 + epsn*n2;

            SIF_collection1 = [];
            Upper_GB        = {};

            for i = 1:numel(targetCell.faces)
                key = targetCell.faceKeys{i};
                if isKey(gbKeys, key), continue; end

                V0 = targetCell.faces{i};
                Wa = clipHalfPlane(V0, n2, p2a, true);
                if ~isempty(Wa) && size(Wa,1) >= 3
                    [Wa, ~] = finding_face_center(Wa);
                    Upper_GB{end+1} = Wa;
                    [~, ~, SIF1] = finding_gb_plane_angles(Wa, Gamma_G);
                    SIF_collection1(end+1) = SIF1;
                end
            end

            if isempty(Upper_GB), continue; end
            GB1.K = max(SIF_collection1);

            if numel(Upper_GB) == 1
                facesSel = Upper_GB;
                if GB1.K <= T2.K
                    Kmin    = SIF_collection1;
                    modeStr = 'GB';
                    polySel = [];
                else
                    Kmin    = T2.K;
                    modeStr = 'T2';
                    polySel = [];
                end
                ID(end+1)         = id;
                K_selected(end+1) = Kmin;
                PathLog(end+1)    = struct('cell',id,'step',w,'mode',modeStr, ...
                    'K_sel',Kmin,'K_T1',T1.K,'K_T2',T2.K,'K_GB1',GB1.K,'K_GB2',NaN, ...
                    'poly',polySel,'faces',{facesSel});
                continue
            end

            GB_wall = {};
            for i = 1:numel(targetCell.faces)
                key = targetCell.faceKeys{i};
                if isKey(gbIndex, key)
                    GB_wall{end+1} = GB{gbIndex(key)};
                end
            end

            margin = 5e-6;
            ek = @(a,b) mat2str(sortrows(round([a;b], 6)));

            allUpperEdges = {};
            allUpperKeys  = {};

            for k = 1:numel(GB_wall)
                Vg = GB_wall{k}.vertices;
                Wa = clipHalfPlane(Vg, n2, p2a, true);
                if isempty(Wa) || size(Wa,1) < 2, continue; end
                [Wa, ~] = finding_face_center(Wa);
                m = size(Wa, 1);
                for e = 1:m
                    a   = Wa(e,:);
                    b   = Wa(mod(e,m)+1,:);
                    mid = 0.5*(a+b);
                    if abs(dot(n2, mid-p2)) <= margin, continue; end
                    allUpperEdges{end+1} = [a; b];
                    allUpperKeys{end+1}  = ek(a, b);
                end
            end

            if isempty(allUpperEdges), continue; end

            [~, ~, ic] = unique(allUpperKeys);
            counts     = accumarray(ic(:), 1);
            UpperExposedEdges = {};
            for q = 1:numel(allUpperEdges)
                if counts(ic(q)) == 1
                    UpperExposedEdges{end+1} = allUpperEdges{q};
                end
            end

            if isempty(UpperExposedEdges), continue; end

            Pall  = cell2mat(UpperExposedEdges(:));
            z_avg = mean(Pall(:,3));

            pts_sec = intersect_plane_with_faces(targetCell.faces, [0 0 1], [0 0 z_avg], 1e-10);
            if size(pts_sec,1) >= 3
                [Psec, ~] = finding_face_center(pts_sec);
                centroid  = mean(Psec, 1);
            else
                continue
            end

            zu = Pall(:,3);
            [zu_sorted, ~] = sort(zu, 'ascend');
            if numel(zu_sorted) >= 2
                minZ_upper = zu_sorted(2);
            else
                minZ_upper = zu_sorted(1);
            end
            maxZ_upper = max(zu);

            zstart = 1.05 * minZ_upper;
            zmax   = 0.95 * maxZ_upper;
            Nz     = 10;
            zgrid  = linspace(zstart, zmax, Nz);

            p2s = p2; p2s(3) = zstart;
            p1s = p1; p1s(3) = zstart;

            planeSweep  = struct('z',{},'p',{},'poly',{},'touched_idx',{},'SIFs',{},'SIFmax',{},'touched_faces',{});
            planeSweep1 = struct('z',{},'p',{},'poly',{},'touched_idx',{},'SIFs',{},'SIFmax',{},'touched_faces',{});

            for kk = 1:Nz
                if abs(n2(3)) > 1e-12
                    pk = p2s + ((zgrid(kk)-zstart)/n2(3)) * n2;
                else
                    pk = p2s; pk(3) = zgrid(kk);
                end

                if abs(n1(3)) > 1e-12
                    pk1 = p1s + ((zgrid(kk)-zstart)/n1(3)) * n1;
                else
                    pk1 = p1s; pk1(3) = zgrid(kk);
                end

                pts  = intersect_plane_with_faces(targetCell.faces, n2,  pk,  tol);
                pts1 = intersect_plane_with_faces(targetCell.faces, n1, pk1, tol);

                poly  = [];  poly1 = [];
                if size(pts, 1) >= 3, poly  = order_plane_polygon(pts,  n2); end
                if size(pts1,1) >= 3, poly1 = order_plane_polygon(pts1, n1); end

                touched_idx = [];  touched_faces = {};
                for i = 1:numel(targetCell.faces)
                    V0 = targetCell.faces{i};
                    if isBoundaryFace(round(V0,6), bx, tol), continue; end
                    d = (V0 - pk) * n2';
                    if (min(d) <= tol) && (max(d) >= -tol)
                        touched_idx(end+1)   = i;
                        touched_faces{end+1} = V0;
                    end
                end
                SIFs_UT2 = T2.K;
                planeSweep(kk) = struct('z',zgrid(kk),'p',pk,'poly',poly, ...
                    'touched_idx',touched_idx,'SIFs',SIFs_UT2,'SIFmax',SIFs_UT2, ...
                    'touched_faces',{touched_faces});

                touched_idx1 = [];  touched_faces1 = {};  SIFs_UT1 = [];
                for i1 = 1:numel(targetCell.faces)
                    V0 = targetCell.faces{i1};
                    if isBoundaryFace(round(V0,6), bx, tol), continue; end
                    d = (V0 - pk1) * n1';
                    if (min(d) <= tol) && (max(d) >= -tol)
                        touched_idx1(end+1)   = i1;
                        touched_faces1{end+1} = V0;
                        [~,~,SIF_t1] = finding_gb_plane_angles(V0, Gamma);
                        SIFs_UT1(end+1,1) = SIF_t1;
                    end
                end
                SIFs_UT1 = [SIFs_UT1; T1.K];
                Smax1    = max(SIFs_UT1);
                planeSweep1(kk) = struct('z',zgrid(kk),'p',pk1,'poly',poly1, ...
                    'touched_idx',touched_idx1,'SIFs',SIFs_UT1,'SIFmax',Smax1, ...
                    'touched_faces',{touched_faces1});
            end

            Svec = arrayfun(@(s) s.SIFmax, planeSweep);
            [~, bestk] = min(Svec);
            best = planeSweep(bestk);

            normidx = @(v) sort(unique(v(~isnan(v))));
            ref  = normidx(best.touched_idx);
            keep = false(1, numel(planeSweep1));
            for mn = 1:numel(planeSweep1)
                cur      = normidx(planeSweep1(mn).touched_idx);
                keep(mn) = all(ismember(cur, ref));
            end
            planeSweep1 = planeSweep1(keep);

            best1 = [];
            if ~isempty(planeSweep1)
                Svec1 = arrayfun(@(s) s.SIFmax, planeSweep1);
                [~, bestk1] = min(Svec1);
                best1 = planeSweep1(bestk1);
            end

            if isempty(best1)
                Kvals_sel = [max(SIF_collection1), best.SIFmax];
            else
                Kvals_sel = [max(SIF_collection1), best.SIFmax, best1.SIFmax];
            end

            Kmin = min(Kvals_sel);
            ties = find(abs(Kvals_sel - Kmin) <= tol);

            if any(ties==1) && any(ties==2)
                idx = 2;
            else
                pref = inf(1, numel(Kvals_sel));
                [~, ii] = min(pref(ties));
                idx = ties(ii);
            end

            switch idx
                case 1
                    modeStr  = 'GB';
                    polySel  = [];
                    facesSel = Upper_GB;
                    ID(end+1)         = id;
                    K_selected(end+1) = max(SIF_collection1);

                case 2
                    triFaces = {};
                    for i = 1:numel(UpperExposedEdges)
                        seg  = UpperExposedEdges{i};
                        triV = [seg(1,:); seg(2,:); centroid];
                        [triV, ~] = finding_face_center(triV);
                        triFaces{end+1} = triV;
                    end
                    modeStr  = 'T2';
                    polySel  = best.poly;
                    facesSel = triFaces;
                    ID(end+1)         = id;
                    K_selected(end+1) = best.SIFmax;

                case 3
                    if isempty(best1), continue; end
                    modeStr  = 'T1';
                    polySel  = best1.poly;
                    facesSel = best1.touched_faces;
                    ID(end+1)         = id;
                    K_selected(end+1) = best1.SIFmax;
            end

            PathLog(end+1) = struct('cell',id,'step',w,'mode',modeStr,'K_sel',Kmin, ...
                'K_T1',T1.K,'K_T2',T2.K,'K_GB1',GB1.K,'K_GB2',NaN, ...
                'poly',polySel,'faces',{facesSel});

        end

        t_total_sec    = toc(t_total);
        t_fracture_sec = t_total_sec - t_voronoi_sec;

        if ~isempty(PathLog)
            modes = {PathLog.mode};
            isT1  = strcmp(modes,'T1');
            isT2  = strcmp(modes,'T2');
            isGB  = strcmp(modes,'GB');
            nT1 = sum(isT1);  nT2 = sum(isT2);  nGB = sum(isGB);
            tot  = numel(modes);
            fprintf('T1=%d (%.1f%%), T2=%d (%.1f%%), GB=%d (%.1f%%)\n', ...
                nT1,100*nT1/tot, nT2,100*nT2/tot, nGB,100*nGB/tot);
            K_T1 = K_selected(isT1);
            K_T2 = K_selected(isT2);
            K_GB = K_selected(isGB);
        else
            modes = {}; isT1=[]; isT2=[]; isGB=[]; K_T1=[]; K_T2=[]; K_GB=[];
        end

        disp(t_total_sec);

        if isempty(K_selected)
            Kmax = NaN;  idxMax = [];  srcMax = 'N/A';
        else
            Kmax   = max(K_selected);
            idxMax = find(abs(K_selected - Kmax) <= tol);
            srcMax = PathLog(idxMax(1)).mode;
        end

        cellsList = ID(:);
        n = numel(cellsList);

        if n == 0
            RowStore = struct('numRows',0,'cellsPerRow',[],'RowOrderCells',[], ...
                'RowOrderK',[],'rowBuckets',{{}},'visitOrder',[], ...
                'KR_curve',[],'a_over_L',[]);
            SimStore(end+1) = struct('rng_index',rindex,'Num_Seed',numPoints,'z',zSlice, ...
                'PathLog',PathLog,'Kmax',Kmax,'K_T1',K_T1,'K_T2',K_T2,'K_GB',K_GB, ...
                'srcMax',srcMax,'idxMax',idxMax,'I_Cells',intersectingCells, ...
                'M_time',t_voronoi_sec,'T_time',t_total_sec,'RowData',RowStore);
            continue
        end

        id2idx = containers.Map(num2cell(cellsList), num2cell(1:n));

        E = [];
        for k = 1:numel(sharedFaces)
            c = intersect(sharedFaces{k}.cells, cellsList');
            if numel(c) >= 2
                [I,J] = find(triu(true(numel(c)),1));
                E = [E; [c(I) c(J)]];
            end
        end

        Adj = sparse(n, n);
        for r = 1:size(E,1)
            a = id2idx(E(r,1));  b = id2idx(E(r,2));
            Adj(a,b) = 1;        Adj(b,a) = 1;
        end

        G = graph(Adj);

        C_mat = cell2mat(arrayfun(@(cid) voronoiCells(cid).center, cellsList, 'UniformOutput', false));
        xy = C_mat(:,1:2);

        ymin    = min(xy(:,2));  ymax = max(xy(:,2));
        numRows = max(1, round(sqrt(n)));
        edgesY  = linspace(ymin, ymax+(ymax==ymin), numRows+1);
        rowIdx  = discretize(xy(:,2), edgesY);

        rowBuckets = cell(numRows,1);
        for r = 1:numRows
            ids = find(rowIdx==r);
            if isempty(ids), rowBuckets{r} = []; continue; end
            [~, ord] = sort(xy(ids,1),'ascend');
            rowBuckets{r} = ids(ord);
        end

        if ~exist('chosenCell','var')
            X = [bx/2  0  bx/2];
            [~, ix] = min(vecnorm(seeds(cellsList,:) - X, 2, 2));
            chosenCell = cellsList(ix);
        end

        cur        = id2idx(chosenCell);
        visitOrder = cur;
        visited    = false(n,1);
        visited(cur) = true;

        for r = 1:numRows
            targets = rowBuckets{r};
            targets = targets(~visited(targets));
            for t = 1:numel(targets)
                dst = targets(t);
                if ~visited(dst)
                    if Adj(cur,dst) == 1
                        visitOrder(end+1) = dst;
                        visited(dst) = true;
                        cur = dst;
                    else
                        p = shortestpath(G, cur, dst);
                        p = p(~visited(p));
                        if ~isempty(p)
                            visitOrder = [visitOrder, p(:)'];
                            visited(p) = true;
                            cur = p(end);
                        end
                    end
                end
            end
        end

        if any(~visited)
            leftovers = find(~visited);
            for k = 1:numel(leftovers)
                dst = leftovers(k);
                p   = shortestpath(G, cur, dst);
                p   = p(~visited(p));
                if ~isempty(p)
                    visitOrder = [visitOrder, p(:)'];
                    visited(p) = true;
                    cur = p(end);
                end
            end
        end

        RowOrderCells = cellsList(visitOrder);

        Kvals_row = inf(n,1);
        for t = 1:numel(PathLog)
            cid = PathLog(t).cell;
            if isKey(id2idx, cid), Kvals_row(id2idx(cid)) = PathLog(t).K_sel; end
        end
        RowOrderK = Kvals_row(visitOrder);

        rowMaxK = nan(numRows,1);
        a_row   = nan(numRows,1);
        y0    = min(xy(:,2));
        yTop  = max(xy(:,2));
        Lspan = max(yTop - y0, eps);

        for r = 1:numRows
            ids_in_row = find(rowIdx(visitOrder) == r);
            if ~isempty(ids_in_row)
                rowMaxK(r) = max(RowOrderK(ids_in_row));
                rowCells   = RowOrderCells(ids_in_row);
                yc = arrayfun(@(cid) voronoiCells(cid).center(2), rowCells);
                a_row(r)   = mean(yc) - y0;
            end
        end

        mask = ~isnan(rowMaxK) & ~isnan(a_row);
        try
            KR_curve = cummax(rowMaxK(mask));
        catch
            KR_curve = rowMaxK(mask);
            for i = 2:numel(KR_curve), KR_curve(i) = max(KR_curve(i-1), KR_curve(i)); end
        end
        a_over_L = a_row(mask) / Lspan;

        cellsPerRow = cellfun(@numel, rowBuckets).';
        RowStore = struct('numRows',numRows,'cellsPerRow',cellsPerRow, ...
            'RowOrderCells',RowOrderCells,'RowOrderK',RowOrderK, ...
            'rowBuckets',{rowBuckets},'visitOrder',visitOrder, ...
            'KR_curve',KR_curve,'a_over_L',a_over_L);

        SimStore(end+1) = struct('rng_index',rindex,'Num_Seed',numPoints,'z',zSlice, ...
            'PathLog',PathLog,'Kmax',Kmax,'K_T1',K_T1,'K_T2',K_T2,'K_GB',K_GB, ...
            'srcMax',srcMax,'idxMax',idxMax,'I_Cells',intersectingCells, ...
            'M_time',t_voronoi_sec,'T_time',t_total_sec,'RowData',RowStore);

    end

    M2_ArchiveEquR1(end+1) = struct('rindex',rindex,'SimStore',SimStore);
    save('M2_ArchiveEquR1_all.mat','M2_ArchiveEquR1','-v7.3');
    SimStore = [];
    clearvars -except rindex M2_ArchiveEquR1 SimStore t_total

end