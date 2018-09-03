function [connected, synapse] = kalman_synapse(originations, Neurons, R, BW, source, skeleton, thetas)
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

num_bifur = size(originations, 1);
bifur_thre = 15;
% parameters for each bifur 
% current point, former point, path length, path
bifur_paras = cell(num_bifur, 1);
for i = 1 : num_bifur
    bifur_paras{i} = zeros(3, 2);
    bifur_paras{i}(1, :) = originations(i, :); 
    bifur_paras{i}(2, :) = Neurons(source, :);
end
% bifur_path = cell(num_bifur, 1);

bifur_dis = 1;
[M, N] = size(skeleton);

% for compare
% if line is going through the source neuron, stop it
% to prevent mis-stop, set R(source) smaller
R_compare = 1.2 * R';

move_step = 1.5;
move_var = 2;
measure_var = 1;

while num_bifur > 0
    % get information from bifur parameters
    current = bifur_paras{1}(1, :);
    former = bifur_paras{1}(2, :);
    path = bifur_paras{1}(4:end, :);
    path_length = bifur_paras{1}(3, 1);
    % init P, Q, R in kalman filter
    P = measure_var * 2;
    Q = move_var;
    R = measure_var;
    while true
        plot([current(1), former(1)],[current(2), former(2)],'LineWidth',2,'Color','blue');
        skeleton(current(2), current(1)) = false;
        path = [path; current];
        % direction of last movement
        former_theta = atan2d(current(2) - former(2), current(1) - former(1));
        % whether there are bifurcations
        % meshgrid
        [xx, yy] = produce_meshgrid(current, bifur_dis, M, N);
        theta_mesh = atan2d(yy - current(2), xx - current(1));
        % forward direction
        theta_mesh_dis = abs(theta_mesh - former_theta);
        theta_mesh_dis(theta_mesh_dis > 180) = 360 - theta_mesh_dis(theta_mesh_dis > 180);
        forward_mesh = theta_mesh_dis < 140;
        % available direction from skeleton
        available_mesh = skeleton(yy(:,1), xx(1,:));
        available_forward = forward_mesh & available_mesh;
        available_forward(5) = false;
        % path terminate
        if ~sum(available_forward(:))
            num_bifur = num_bifur - 1;
            bifur_paras(1) = [];
            break;
        end
        neurite_theta = thetas(yy(:,1), xx(1,:));
        neurite_available_forward = neurite_theta(available_forward);
        if (range(neurite_available_forward(:)) > bifur_thre)
            % num of bifurs at next step
            current_bifur = length(neurite_available_forward);
            % num of all bifurs
            num_bifur = num_bifur + current_bifur - 1;
            new_bifur_paras = cell(current_bifur, 1);
            % value of bifurs
            for i = 1:current_bifur
                bifur_theta = neurite_available_forward(i);
                bifur_index = find(neurite_theta == bifur_theta);
                new_start = [xx(bifur_index), yy(bifur_index)];
                skeleton(new_start(2), new_start(1)) = false;
                new_bifur_paras{i} = [new_start; current; ...
                                      path_length + 1, 0; [path; new_start]];
            end
            bifur_paras(1) = [];
            bifur_paras = [new_bifur_paras; bifur_paras];
            break;
        end
        % angle of current point from neurite detector
        angle = thetas(current(2), current(1));
        % since angle means two opposite direction, 
        % we should find the direction to forward
        theta1 = angle; 
        theta2 = - sign(theta1) * (180 - abs(angle));
        theta_dis = abs([theta1, theta2] - former_theta);
        theta_dis(theta_dis > 180) = 360 - theta_dis(theta_dis > 180);
        [~, index_theta] = min(theta_dis);
        forward_theta = (index_theta == 1) * theta1 + (index_theta == 2) * theta2;
        % predict next point: predict = current + step * theta
        predict = current + move_step * [cosd(forward_theta), sind(forward_theta)];
        predict = round(predict);
        P = P + Q;
        plot(predict(1),predict(2),'.','color','red','MarkerSize', 10);
        % from prediction, get measurement: nearest point in BW_thin
        measure = find_near_center(skeleton, predict);
        plot(measure(1),measure(2),'.','color','green','MarkerSize', 10);
        % kalman scale
        K = P / (P + R);
        % estimate = predict + K * residual
        estimate = predict + K * (measure - predict);
        estimate = round(estimate);
        P = (1 - K) * P;
        former = current;
        current = estimate;
        path_length = path_length + sum((former - current) .^ 2);
        % whether reach other neurons
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
    end
end

    function center_p = find_near_center(skeleton, current)
        r = 1;
        [m, n] = size(skeleton);
        while true
            [x_grid, y_grid] = produce_meshgrid(current, r, m, n);
            index_grid = sub2ind(size(skeleton), y_grid, x_grid);
            if sum(sum(skeleton(index_grid))) > 0
                index_skeleton = find(skeleton(index_grid));
                dis_grid = (x_grid(index_skeleton) - current(1)).^2 ...
                         + (y_grid(index_skeleton) - current(2)).^2;
                [~, center_index] = min(dis_grid);
                center_index = index_skeleton(center_index);
                center_p = [x_grid(center_index), y_grid(center_index)];
                break; 
            else
                r = r + 1;
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