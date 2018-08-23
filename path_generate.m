function [target, path] = path_generate(current, former, Neurons, R, source, skeleton, COVER)
% current: current pixel on the path (start point)
% former: former pixel on the path (tell the forward direction)
% source: the neuron(index) this path originates from 
% Neurons & R: decide the terminal of path(get other neurons / dead path)

%% Variable Definition
% num of Neurons(R)
num = length(R);
% suppose that finally get one path
num_path = 0;
% init path with start point
current_path = zeros(0, 2);
current_path(1, :) = current;
length_path = 1;
path = [];
% init target neuron
target = [];

fill_gap = 10;
theta_thre = 10;
[M, N] = size(skeleton);
%% Recursion: DFS
while true
    % plot current
    COVER(current(2), current(1)) = true;
    plot(current(1),current(2),'.','color','blue','MarkerSize', 5);
    theta = atan2d(current(2) - former(2), current(1) - former(1));
    % get neighbours
    [neigh_x, neigh_y] = meshgrid(current(1) + (-1:1), current(2) + (-1:1));
    exceed_index = neigh_x < 1 | neigh_x > M | neigh_y < 1 | neigh_y > N; 
    neigh_x(exceed_index) = [];
    neigh_y(exceed_index) = [];
    index = sub2ind([M, N], neigh_y, neigh_x);
    neigh_index = find(skeleton(index)...
        & ~COVER(index));
    % exclude itself
    % find feasible direction 
    neigh_x = neigh_x(neigh_index);
    neigh_y = neigh_y(neigh_index);
    % number of bifurcations
    bifur = length(neigh_x);
    %% if bifurcation > 1, split into a few paths
    if bifur > 1
        theta_neigh = atan2d(neigh_y - current(2), neigh_x - current(1));
        % choose the outward branch that has the mininum included angle with the
        % inward branche in critical areas
        del_theta = abs(theta_neigh - theta);
        del_theta(del_theta > 180) = 360 - del_theta(del_theta > 180);
        [~, del_index] = sort(del_theta);
        % path before bifurcation
        former = current;
        % cover all bifurcations to prevent confusion
        for b = 1:bifur
            COVER(neigh_y(b), neigh_x(b)) = true;
        end
        % depth-first search
        for b = 1:bifur
            forward = del_index(b);
            current = [neigh_x(forward), neigh_y(forward)];
            % recursion
            [new_target, path_follow] = ...
                path_generate(current, former, Neurons, R, source, skeleton, COVER);
            % child-path gets to other neurons
            if ~isempty(new_target)
                % get a new path, re-init path and target
                current_path = current_path(current_path(:,1) > 0, :);
                new_path = cell(length(path) + length(path_follow), 1);
                for i = 1:length(path_follow)
                    new_path{i} = [current_path; path_follow{i}];
                end
                if num_path > 0
                    for j = 1:num_path
                        new_path{i+j} = path{j};
                    end
                end
                path = new_path;
                num_path = num_path + length(path_follow);
                target = [new_target, target];
            end
        end
        break;
    else
        %% no bifurcation, go on with the path
        next = [neigh_x, neigh_y];
        if isempty(next)
            row = max(1, current(1) - fill_gap):...
                  min(N, current(1) + fill_gap);
            col = max(1, current(2) - fill_gap):...
                  min(M, current(2) + fill_gap);
            [area_x, area_y] = meshgrid(row, col);
            dist = sqrt((area_x - current(1)).^2 + (area_y - current(2).^2));    
            in_circle = dist <= fill_gap;
            thetas = atan2d(area_y - current(2), area_x - current(1));
            in_direction = abs(thetas - theta) < theta_thre;
            available = in_circle & in_direction & skeleton(col, row) & ~COVER(col, row);
            if sum(available(:)) > 0
                dist_available = dist ./ available;
                [~, index_forward] = min(dist_available(:));
                next = [area_x(index_forward), area_y(index_forward)];
            else
                break;
            end
        end
        % extend lines
        length_path = length_path + 1;
        current_path(length_path, :) = current;
        % change former and current points
        former = current;
        current = next;
        %% decide whether to get other neurons
        dis = Neurons - repmat(current, num, 1);
        dis = sqrt(sum(dis.^2, 2));
        % for compare
        % if line is going through the source neuron, stop it
        % to prevent mis-stop, set R(source) smaller
        R_compare = R;
        R_compare(source) = 1.5 * R_compare(source);
        if max(dis < R_compare)
            target = find(dis < R_compare);
            if target == source
                target = [];
                break;
            end
            path = cell(1,1);
            path{1} = current_path;
            break;
        end
    end
end
% end function
end