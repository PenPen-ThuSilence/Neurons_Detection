function [connected, synapse, breadth] = synapse_detection(BW, Neurons, R, ...
                                                theta_thre, fill_gap)
%% assign possibility and direction for synapse
BW_thin = bwmorph(BW, 'thin', 20);

sigma = 1;
[U, V, theta] = neurite_vector(BW_thin, sigma);

theta = theta .* BW_thin;

%% find synaspe with thinned image
num = size(Neurons, 1);
connected = zeros(num, num);
synapse = cell(num);
[m, n] = size(BW_thin);

circle_area = circle_points(Neurons, round(R*1.3), BW_thin);

for k = 1:num
    % primary queue: intersection points of circle and neuron
    start_points = circle_area{k};
    index = sub2ind([m, n], start_points(:,2), start_points(:,1));
    start_points = start_points(BW_thin(index), :);
    start_points = merge_random(start_points, 10);
    
    %% PLOT PROCESS
% %     % show primary queue
%     draw_circles_k(Neurons, R, BW_thin, k); hold on;
%     % label neurons
%     for j = 1:num
%         text(Neurons(j,1),Neurons(j,2),int2str(j),'FontSize',15,'Color','yellow');
%     end
%     plot(start_points(:,1),start_points(:,2),'.','color','red','MarkerSize', 15);
% %     quiver(U, V);

    %% find path from primary queue points
    [connected_k, synapse_k] = find_synapse(start_points, Neurons, R, k,...
                                            BW_thin, theta, theta_thre, fill_gap);
    connected(k, :) = connected_k;
    synapse(k, :) = synapse_k;
end

% make sure connectivity is symmetrical and least length
for i = 1:num-1
    for j = i+1:num
        if connected(i, j) && connected(j, i) && connected(i, j) > connected(j, i)...
                || connected(j, i) && ~connected(i, j)
            connected(i, j) = connected(j, i);
            synapse(i,j) = synapse(j, i);
        end    
    end
end

breadth = synapse_breadth(synapse, BW);