%% read the image
load samples;
load label_neurons;
image = samples{1};
label = label_neurons{1};

image = image(501:1000, 501:1000);
%% preprocess
% binarize
BW = image > 30;

% remove dots
dots = image == 255;
% exclude synapses
min_area = 25;
filtered = filterRegions_area(dots, min_area);
% minL = 15;
% filtered = filterRegions_MajorAxis(dots, minL);
dots_filtered = dots & ~filtered;

BW_p = BW & ~dots_filtered;

% dilate and erode to connect parts of neurons
se = strel('disk', 2);
BW_p = imclose(BW_p, se);

% fill holes in neurons which caused when removing dots
BW_filled = imfill(BW_p, 'holes');
fill_differ = BW_filled & ~BW_p;
% preventing from fill very large holes
min_area = 50;
fill = fill_differ & ~filterRegions_area(fill_differ, min_area);

BW_p = BW_p | fill;
% thin synapse
% BW_p = bwmorph(BW_p, 'thin');

% figure;imshowpair(BW_p, image > 40);
% title('fill holes');
%% potential neurons
% radius range
R_center = 40;
R_range = 20;
Neurons = zeros(0, 2);

parfor r = R_center - R_range:R_center + R_range
    % circular kernel: pixels inside circle(r) are 1, outside ones are 0.
    kernel = circle_kernel(r, 0);
    
    % conv2
    density = conv2(double(BW_p), kernel, 'same');
    
    % pixels where conv2 result are larger, which means bright areas
    neu = density > 0.4 * sum(kernel(:));
    
    % find centroids of these bright areas
    stats = regionprops(neu, 'Centroid');
    len = length(stats);
    Centroids = zeros(len, 2);
    for i = 1:len
        Centroids(i, :) = stats(i).Centroid;
    end
    Neurons = [Neurons; Centroids]; 
end

% centroids of bright areas are potential neurons
draw_neurons(BW_p, Neurons);

%% get neurons from centroids
% This section tries to select neurons from points we got from lines. The
% main idea is to calculate numbers of bright angles in an anulus. For each
% degree x in 0:360, if there are bright pixels in the anulus whose angles 
% from the center of the anulus equal to x,we say x is a bright degree. For a point 
% p, A(p,R,annulus) draws a anulus around p with radius = R and width =
% annulus * 2. If number of bright degrees in A(p,R,annulus) is more than 360 *
% threshold_angle, we say point p is a neuron and R is the neuron's size.
% For every point, we try R in dimension 
% [R - R_range, R + R_range] (that is 45:105, if using parameters below).

% Besides, we use another R_around to determine whether there is a bright
% area around the center of neuron, because most neurons are bright in the
% center.

threshold_angle = 0.7;
R_center = 45;
R_range = 25;
R_around = 25;
threshold_around = 0.3;

[final_Neurons, grades, R, around] = IsNeurons_new_4(BW_p, Neurons, ...
                    'threshold_angle', threshold_angle, ...
                    'merge_dis', 2, ...
                    'R', R_center, 'R_range', R_range, 'annulus', 5, ...
                    'R_around', R_around, 'threshold_around', threshold_around);

draw_circles(final_Neurons, R, BW_p);

% number neurons
for i = 1:length(grades)
    text(final_Neurons(i,1),final_Neurons(i,2),int2str(i),'FontSize',20,'Color','red');
end
%% assign possibility and direction for synapse
BW_thin = bwmorph(BW_p, 'thin', 20);

img = image .* uint8(BW_thin);

sigma = 1;
[U, V, theta] = neurite_vector(BW_thin, sigma);

theta = theta .* BW_thin;

% show
figure; imshow(BW_thin);
hold on;
quiver(U, V);

%% synaspe with skeleton and kalman filter
Neurons = final_Neurons;
num = size(Neurons, 1);
connected = zeros(num, num);
synapse = cell(num);
[m, n] = size(BW_thin);

circle_area = circle_points(Neurons, round(R*1.2), BW_thin);

% for k = 1:num
%     points = circle_area{k};
%     plot(points(:,1),points(:,2),'.','color','green', 'MarkerSize', 5); 
% end

for k = 3:num
    % primary queue: intersection points of circle and neuron
    start_points = circle_area{k};
    index = sub2ind([m, n], start_points(:,2), start_points(:,1));
    start_points = start_points(BW_thin(index), :);
    start_points = merge_random(start_points, 10);
    % show primary queue
    draw_circles_k(Neurons, R, BW_thin, k); hold on;
    % label neurons
    for j = 1:num
        text(Neurons(j,1),Neurons(j,2),int2str(j),'FontSize',15,'Color','yellow');
    end
    plot(start_points(:,1),start_points(:,2),'.','color','red','MarkerSize', 15);
    quiver(U, V);
    % find path from primary queue points
    [connected_k, synapse_k] = kalman_synapse(start_points, Neurons, R, k, BW_thin, theta);
    connected(k, :) = connected_k;
    synapse(k, :) = synapse_k;
end

breadth = synapse_breadth(synapse, BW_p);
draw_synapse(connected, synapse, image, breadth);
 %% find synapses
% Neurons = final_Neurons;
% num = size(Neurons, 1);
% connected = false(num, num);
% synapse = cell(num);
% [m, n] = size(BW_thin);
% 
% circle_area = circle_points(Neurons, round(R*1.5), BW_thin);
% 
% % for k = 1:num
% %     points = circle_area{k};
% %     plot(points(:,1),points(:,2),'.','color','green', 'MarkerSize', 5); 
% % end
% 
% for k = 1:num
%     COVER = false(size(BW_thin));
%     % primary queue: intersection points of circle and neuron
%     start_points = circle_area{k};
%     index = sub2ind([m, n], start_points(:,2), start_points(:,1));
%     start_points = start_points(BW_thin(index), :);
%     start_points = merge_random(start_points, 10);
%     % cover all primary start points
%     index = sub2ind([m, n], start_points(:,2), start_points(:,1));
%     COVER(index) = true;
%     % show primary queue
%     draw_circles_k(Neurons, R, BW_thin, k); hold on;
%     % label neurons
%     for j = 1:num
%         text(Neurons(j,1),Neurons(j,2),int2str(j),'FontSize',15,'Color','yellow');
%     end
%     plot(start_points(:,1),start_points(:,2),'.','color','red','MarkerSize', 15);
%     % find path from primary queue points
%     for i = 1:length(start_points)
%         % iterate to get path
%         current = start_points(i, :);
%         former = Neurons(k, :);
%         [target, path] = ...
%             path_generate(current, former, Neurons, R', k, BW_thin, COVER);
%         % find shortest path for each target
%         table = tabulate(target);
%         if ~isempty(target)
%             if max(table(:,2)) > 1 % the number of paths for certain target > 1
%                 index_mul = find(table(:,2) > 1);
%                 for j = 1:length(index_mul)
%                     index_path = find(target == table(index_mul(j),1));
%                     length_path = zeros(length(index_path),1);
%                     for jj = 1:length(index_path)
%                         length_path(jj) = size(path{index_path(jj)},1);
%                     end
%                     [~,min_index] = min(length_path);
%                     path(index_path(1)) = path(index_path(min_index));
%                 end
%             end
%             [target, u_index] = unique(target);
%             connected_k = connected(k, :);
%             synapse_k = synapse(k, :);
%             for j = 1:length(target)
%                 if connected_k(target(j)) && ...
%                         size(path{u_index(j)},1) < size(synapse_k{target(j)}, 1)
%                     % found shorter path
%                     original_path = synapse_k{target(j)};
%                     plot(original_path(:,1),original_path(:,2),'.','color','blue','MarkerSize', 3);
%                     synapse_k{target(j)} = path{u_index(j)};
%                     new_path = synapse_k{target(j)};
%                     plot(new_path(:,1),new_path(:,2),'.','color','green','MarkerSize', 3);
%                     continue;
%                 end
%                 connected_k(target(j)) = true;
%                 synapse_k{target(j)} = path{u_index(j)};
%                 path_j = path{u_index(j)};
%                 plot(path_j(:,1),path_j(:,2),'.','color','green','MarkerSize', 3);
%             end
%             connected(k, :) = connected_k;
%             synapse(k, :) = synapse_k;
%         end
%     end
% end