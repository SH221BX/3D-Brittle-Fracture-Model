function [allBasis,planes,rotated_planes,K_eff] = finding_basis(Ax, Ay, Az, Gamma)

    planes = [];
    hundred_type_planes = [1 0 0; 0 1 0; 0 0 1]; 
    one_eleven_type_planes = [1 1 1; -1 1 1; 1 -1 1; -1 -1 1];
    one_eleven_type_planes = []; 
    one_ten_type_planes = [   0   1  -1;
    0   1   1;
    1  -1   0;
    1   0  -1;
    1   0   1;
    1   1   0];

     %one_ten_type_planes  = [];

    planes = [hundred_type_planes; one_ten_type_planes; one_eleven_type_planes];
    theta = zeros(size(planes,1),1);
    phi = zeros(size(planes,1),1);
    K_eff = zeros(size(planes,1),1);

    rotation_theta_x = Ax;
    rotation_theta_y = Ay;
    rotation_theta_z = Az;

    R_theta_z = [cos(rotation_theta_z) -sin(rotation_theta_z) 0; sin(rotation_theta_z) cos(rotation_theta_z) 0; 0 0 1];
    R_theta_y = [cos(rotation_theta_y) 0 sin(rotation_theta_y); 0 1 0; -sin(rotation_theta_y) 0 cos(rotation_theta_y)];
    R_theta_x = [cos(rotation_theta_x) -sin(rotation_theta_x) 0; sin(rotation_theta_x) cos(rotation_theta_x) 0; 0 0 1];

    rotated_planes = (R_theta_z * R_theta_y * R_theta_x * planes')';

    z_axis = [0 0 1];
    start_index_110 = size(hundred_type_planes,1) + 1;
    end_index_110 = start_index_110 + size(one_ten_type_planes,1) - 1;
    start_index_111 = end_index_110 + 1;
    end_index_111 = start_index_111 + size(one_eleven_type_planes,1) - 1;

    x_cartesian = sin(deg2rad(phi)) .* cos(deg2rad(theta));
    y_cartesian = sin(deg2rad(phi)) .* sin(deg2rad(theta));
    z_cartesian = cos(deg2rad(phi));

    CE_100 = 3.00; CE_110 = 6.48; CE_111 = 8.21;
    CE = zeros(size(planes,1),1);
    CE(1:size(hundred_type_planes,1)) = CE_100;

    if ~isempty(one_ten_type_planes)
        CE(start_index_110:end_index_110) = CE_110;
    end
    if ~isempty(one_eleven_type_planes)
        CE(start_index_111:end_index_111) = CE_111;
    end

    E = 568*10^9; nu = 0.25;
    O = [0 0 0];
    tolerance = 1e-5; allBasis = {};

    for i = 1:size(rotated_planes,1)
        vec = rotated_planes(i,:);
        r = norm(vec);
         phi(i) = acosd(dot(z_axis,vec)/r);
         K_eff(i) = Gamma/abs(cosd(phi(i)));
        planeDir = rotated_planes(i,:);

        if i==1
            basis1 = rotated_planes(2,:);
            basis2 = rotated_planes(3,:);
        elseif i==2
            basis1 = rotated_planes(1,:);
            basis2 = rotated_planes(3,:);
        elseif i==3
            basis1 = rotated_planes(1,:);
            basis2 = rotated_planes(2,:);
        else
            R_z = [cosd(90) sind(90) 0; -sind(90) cosd(90) 0; 0 0 1];
            basis1 = planeDir*R_z;
            basis2 = cross(planeDir,basis1);
        end

        basis1 = basis1/norm(basis1); basis2 = basis2/norm(basis2);  allBasis{i} = [basis1; -basis2];
        planePoints = 12*[-basis1 - basis2; basis1 - basis2; basis1 + basis2; -basis1 + basis2];
        % fill3(planePoints(:,1),planePoints(:,2),planePoints(:,3),'r','FaceAlpha',.25); hold on;
        % xlim([-1.5 1.5]); ylim([-1.5 1.5]); zlim([-1.5 1.5]);
        % axis equal; grid off; box off; axis off;
    end

allBasis = allBasis';

end






