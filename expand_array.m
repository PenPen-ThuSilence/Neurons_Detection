function exp_array = expand_array(node_x, node_y, hn, xTarget, yTarget, ...
                                   CLOSED, MAX_X, MAX_Y, former_theta, ...
                                   thetas, fill_gap, theta_thre)
    %Function to return an expanded array
    %This function takes a node and returns the expanded list
    %of successors,with the calculated fn values.
    %The criteria being none of the successors are on the CLOSED list.
    
    Node = [node_x, node_y];
    Target = [xTarget, yTarget];
    
    % meshgrid
    [x_grid, y_grid] = produce_meshgrid(Node, 1, MAX_X, MAX_Y);
    point_num = numel(x_grid);
    exp_points = [reshape(x_grid, [point_num, 1]), reshape(y_grid, [point_num, 1])];
    
    % whether neighbour points in closed list
    Isclose = ismember(exp_points, CLOSED, 'rows');
    exp_points(Isclose, :) = [];
    
    exp_array = exp_points;
    exp_count = size(exp_array, 1);
    
    if ~exp_count
        [area_x, area_y] = produce_meshgrid(Node, fill_gap, MAX_X, MAX_Y);
        point_num = numel(area_x);
        link_points = [reshape(area_x, [point_num, 1]), reshape(area_y, [point_num, 1])];
        
        area_theta = atan2d(link_points(:,2) - Node(2), link_points(:,1) - Node(1));
        % get angle from neurite detector
        angle = thetas(Node(2), Node(1));
        % since angle means two opposite direction, 
        % we should find the direction to forward
        theta1 = angle; 
        theta2 = - sign(theta1) * (180 - abs(angle));
        theta_dis = abs([theta1, theta2] - former_theta);
        theta_dis(theta_dis > 180) = 360 - theta_dis(theta_dis > 180);
        [~, index_theta] = min(theta_dis);
        forward_theta = (index_theta == 1) * theta1 + (index_theta == 2) * theta2;
        % find those points which are in the forward direction
        theta_differ = abs(area_theta - forward_theta);
        in_direction = theta_differ < theta_thre;
        link_points = link_points(in_direction, :);
        % whether neighbour points in closed list
        Isclose = ismember(link_points, CLOSED, 'rows');
        link_points(Isclose, :) = [];
        if size(link_points, 1) > 0
            % continue with the point after gap
            dist = sqrt((link_points(:,2) - Node(2)).^2 + ...
                        (link_points(:,1) - Node(1)).^2);
            [~, index_forward] = min(dist);
            next_after_gap = link_points(index_forward, :);
            exp_array = next_after_gap;
            exp_count = 1;
        end
    end
    
    % calculate hn, gn and fn
    hn_new = hn + sqrt(sum((repmat(Node, exp_count, 1) - exp_array).^2, 2));
    
    gn_new = sqrt(sum((repmat(Target, exp_count, 1) - exp_array).^2, 2));
    
    fn_new = hn_new + gn_new;
    
    exp_array = [exp_array, hn_new, gn_new, fn_new];
end  