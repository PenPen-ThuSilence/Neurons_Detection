function [connected, synapse] = find_synapse(origins, Neurons, R, source,...
                            skeleton, thetas, theta_thre, fill_gap)
% current: current pixel on the path (start point)
% former: former pixel on the path (tell the forward direction)
% source: the neuron(index) this path originates from 
% Neurons & R: decide the terminal of path(get other neurons / dead path)

%% Variable Definition
% num of Neurons(R)
num = length(R);

% init target neuron
connected = zeros(num, 1);
synapse = cell(num, 1);

num_bifur = size(origins, 1);
% parameters for each bifur 
% current point, former point, path length, path
bifur_paras = cell(num_bifur, 1);
for i = 1 : num_bifur
    bifur_paras{i} = zeros(3, 2);
    bifur_paras{i}(3, 2) = 1;
    bifur_paras{i}(1, :) = origins(i, :); 
    bifur_paras{i}(2, :) = Neurons(source, :);
end

[M, N] = size(skeleton);
skeleton_save = skeleton;

% for compare
% if line is going through the source neuron, stop it
% to prevent mis-stop, set R(source) smaller
R_compare = 1.25 * R';
R_compare(source) = 1.15 * R(source);

while num_bifur > 0
    % get information from bifur parameters
    current = bifur_paras{1}(1, :);
    former = bifur_paras{1}(2, :);
    path = bifur_paras{1}(4:end, :);
    path_length = bifur_paras{1}(3, 1);
    refresh_skeleton = bifur_paras{1}(3, 2);
    if refresh_skeleton
        skeleton = skeleton_save;
        bifur_paras{1}(3, 2) = 0;
    end
    % advance along path 1
    while true
        % plot synapses
%         plot([current(1), former(1)],[current(2), former(2)],'LineWidth',2,'Color','blue');
        %% whether reach other neurons
        dis = Neurons - repmat(current, num, 1);
        dis = sqrt(sum(dis.^2, 2));
        if max(dis < R_compare)
            target = find(dis < R_compare);
            num_bifur = num_bifur - 1;
            bifur_paras(1) = [];
            if target == source
                break;
            end   
            if connected(target) == 0
                connected(target) = path_length;
                synapse{target} = path;
            elseif path_length < connected(target)
                connected(target) = path_length;
                synapse{target} = path;
            end
            break;
        end 
        % cover passed points 
        skeleton(current(2), current(1)) = false;
        path = [path; current];
        % direction of last movement
        former_theta = atan2d(current(2) - former(2), current(1) - former(1));
        %% whether there are bifurcations
        % get neighbours
        [neigh_x, neigh_y] = produce_meshgrid(current, 1, M, N);
        index = sub2ind([M, N], neigh_y, neigh_x);
        neigh_index = find(skeleton(index));
        % find feasible direction 
        neigh_x = neigh_x(neigh_index);
        neigh_y = neigh_y(neigh_index);
        % number of bifurcations
        current_bifur = length(neigh_index);
        %% no way to go
        if ~current_bifur
            % whether this is just a small gap ?
            [area_x, area_y] = produce_meshgrid(current, fill_gap, M, N);
            area_theta = atan2d(area_y - current(2), area_x - current(1));
            % get angle from neurite detector
            angle = thetas(current(2), current(1));
            % since angle means two opposite direction, 
            % we should find the direction to forward
            theta1 = angle; 
            theta2 = - sign(theta1) * (180 - abs(angle));
            theta_dis = abs([theta1, theta2] - former_theta);
            theta_dis(theta_dis > 180) = 360 - theta_dis(theta_dis > 180);
            [~, index_theta] = min(theta_dis);
            forward_theta = (index_theta == 1) * theta1 + (index_theta == 2) * theta2;
            theta_differ = abs(area_theta - forward_theta);
            in_direction = theta_differ < theta_thre;
            available = in_direction & skeleton(area_y(:, 1), area_x(1, :));
            if sum(available(:)) > 0
                % continue with the point after gap
                dist = sqrt((area_x - current(1)).^2 + (area_y - current(2)).^2);
                dist_available = dist ./ available;
                [dis_min, index_forward] = min(dist_available(:));
                next_after_gap = [area_x(index_forward), area_y(index_forward)];
                former = current;
                current = next_after_gap;
                path_length = path_length + sqrt(dis_min);
            else
                % end this path
                num_bifur = num_bifur - 1;
                bifur_paras(1) = [];
                break;
            end
        %% more way to go
        elseif current_bifur > 1
            theta_neigh = atan2d(neigh_y - current(2), neigh_x - current(1));
            % sort directions according to the difference between former
            % direction and forward direction.
            del_theta = abs(theta_neigh - former_theta);
            del_theta(del_theta > 180) = 360 - del_theta(del_theta > 180);
            [~, del_index] = sort(del_theta);
            % num of all bifurs
            num_bifur = num_bifur + current_bifur - 1;
            new_bifur_paras = cell(current_bifur, 1);
            % value of bifurs
            for i = 1:current_bifur
                new_start = [neigh_x(del_index(i)), neigh_y(del_index(i))];
                skeleton(new_start(2), new_start(1)) = false;
                new_bifur_paras{del_index(i)} = [new_start; current; ...
                                      path_length + 1, 0; [path; new_start]];
            end
            bifur_paras(1) = [];
            bifur_paras = [new_bifur_paras; bifur_paras];
            break;
        %% one way to go
        else
            former = current;
            current = [neigh_x, neigh_y];
            path_length = path_length + sqrt(sum((current - former) .^ 2));
        end
    end
end

    function [x_grid, y_grid] = produce_meshgrid(center, r, m, n)
        % produce meshgrid with center and range distance
        % x_range
        x = center(1);
        y = center(2);
        x_range = x - r : x + r;
        exceed_x = x_range < 1 | x_range > n;
        x_range(exceed_x) = [];
        % y_range
        y_range = y - r : y + r;
        exceed_y = y_range < 1 | y_range > m;
        y_range(exceed_y) = [];
        [x_grid, y_grid] = meshgrid(x_range, y_range);  
    end
end