%% find synapses
Neurons = final_Neurons;
num = size(Neurons, 1);
connected = false(num, num);
synapse = cell(num);
[m, n] = size(BW_thin);

circle_area = circle_points(Neurons, round(R*1.5), BW_thin);

% for k = 1:num
%     points = circle_area{k};
%     plot(points(:,1),points(:,2),'.','color','green', 'MarkerSize', 5); 
% end

for k = 1:num
    COVER = false(size(BW_thin));
    % primary queue: intersection points of circle and neuron
    queue_p = circle_area{k};
    index = sub2ind([m, n], queue_p(:,2), queue_p(:,1));
    queue_p = queue_p(BW_thin(index), :);
    queue_p = merge_random(queue_p, 10);
    % cover all primary start points
    index = sub2ind([m, n], queue_p(:,2), queue_p(:,1));
    COVER(index) = true;
    % show primary queue
    draw_circles_k(Neurons, R, BW_thin, k); hold on;
    % label neurons
    for j = 1:num
        text(Neurons(j,1),Neurons(j,2),int2str(j),'FontSize',15,'Color','yellow');
    end
    plot(queue_p(:,1),queue_p(:,2),'.','color','red','MarkerSize', 15);
    % find path from primary queue points
    for i = 1:length(queue_p)
        % iterate to get path
        current = queue_p(i, :);
        former = Neurons(k, :);
        [target, path] = ...
            path_generate(current, former, Neurons, R', k, BW_thin, COVER);
        % find shortest path for each target
        table = tabulate(target);
        if ~isempty(target)
            if max(table(:,2)) > 1 % the number of paths for certain target > 1
                index_mul = find(table(:,2) > 1);
                for j = 1:length(index_mul)
                    index_path = find(target == table(index_mul(j),1));
                    length_path = zeros(length(index_path),1);
                    for jj = 1:length(index_path)
                        length_path(jj) = size(path{index_path(jj)},1);
                    end
                    [~,min_index] = min(length_path);
                    path(index_path(1)) = path(index_path(min_index));
                end
            end
            [target, u_index] = unique(target);
            for j = 1:length(target)
                if connected(k, target(j)) && ...
                        size(path{u_index(j)},1) < size(synapse{k, target(j)},1)
                    % find shorter path
                    original_path = synapse{k, target(j)};
                    plot(original_path(:,1),original_path(:,2),'.','color','blue','MarkerSize', 3);
                    synapse{k, target(j)} = path{u_index(j)};
                    new_path = synapse{k, target(j)};
                    plot(new_path(:,1),new_path(:,2),'.','color','green','MarkerSize', 3);
                    continue;
                end
                connected(k, target(j)) = true;
                synapse{k, target(j)} = path{u_index(j)};
                path_j = path{u_index(j)};
                plot(path_j(:,1),path_j(:,2),'.','color','green','MarkerSize', 3);
            end
        end
    end
end