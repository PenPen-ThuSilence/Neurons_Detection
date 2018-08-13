function [connectivity, label] = queue_synapse(Neurons, roi, R)
% function to detect synapse
% with queues

% parameters
num = size(Neurons, 1);
neighbour = [-1 -1;-1 0;-1 1;0 -1;0 1;1 -1;1 0;1 1];
img = roi;
label = double(img > 0);
[m,n] = size(img);
colors = rand(1, 3, num);

% outputs
connectivity = false(num);
connectivity(1:num+1:end) = true;

% expand the circle to get intersection
circle_area = draw_circles(Neurons, round(2*R), label);

% init stacks
queues = cell(num, 1);
% initialize
imshow(label, []); hold on;
axis on;
for i = 1:num
    start_points = circle_area{i};
    index = sub2ind([m, n], start_points(:,2), start_points(:,1));
    start_points = start_points(label(index) == 1,:);
    label(index(label(index) == 1)) = i + 1;

    queues{i} = start_points;
end

finish_flag = zeros(num, 1);
while true
    if min(finish_flag == 1)
        fprintf('Synapse detection terminates!\n');
        break;
    end
    for i = 1:num

        if ~isempty(queues{i})
            % head of queue
            current = queues{i}(1, :);            
            queues{i}(1, :) = [];
            % label head
            label(current(2), current(1)) = i + 1;
            % find neighbours
            neighbours = repmat(current, 8, 1) + neighbour;
            % whether out of range
            index_out_of_range = neighbours(:,1) < 1 | neighbours(:,2) < 1 | ...
                                 neighbours(:,1) > n | neighbours(:,2) > m;
            neighbours(index_out_of_range, :) = [];
            % find bright neighbours
            neibours_index = sub2ind(size(label), neighbours(:,2), neighbours(:,1));
            nei_type = label(neibours_index);
            % type = 1, add to un
            add_index = find(nei_type == 1);
            label(neibours_index(add_index)) = i + 1;
            bright_neigh = neighbours(add_index, :);
            % add bright neighbours to queue
            queues{i} = [queues{i}; bright_neigh];
            % reach other neurons ? type > 1?
            % minus self-type
            nei_type = setdiff(unique(nei_type), i+1);
            % neighbours labeled by other neurons, get synapse between.
            if max(nei_type) > 1
                connectivity(i, max(nei_type) - 1) = true;
                connectivity(max(nei_type) - 1, i) = true;
            end
        else
            finish_flag(i) = 1;
        end
    end
end

% number the neurons
for i = 1:num
    [points_i_x, points_i_y] = find(label == i + 1);
    plot(points_i_y, points_i_x, '.', 'color', colors(:,:,i));
    text(Neurons(i,1),Neurons(i,2),int2str(i),'FontSize',15,'Color','red');
end