%% read image
% Input: filename of image
filename = 'test_sample.tif';
img = imread(filename);
%% Preprocess
scale = 0.5;
level_low = 0.30;
level_high = 0.98;
disk_size = 1;
label_thre = 40;
intensity_thre = 160;
extent_thre = 0.3;
area_thre = 600;
binary_thre = 40;
min_pixel_num = 50;

[image_removed, BW_p] = preprocess(img, scale, level_low, level_high, ...
                                     disk_size, label_thre, intensity_thre, ...
                                     extent_thre, area_thre, binary_thre, min_pixel_num);
% show                                 
figure;
imshow(image_removed);

figure;
imshow(BW_p);
title('binary image');
%% Neuron Detection
% get potential points which can be neurons with kernel
% parameters
R_min = 20;
R_max = 60;
density_thre = 0.35:0.05:1;

potential_neurons = intensive_points(BW_p, R_min, R_max, density_thre);
draw_neurons(BW_p, potential_neurons);
% get neurons from centroids
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

% parameters
threshold_angle = 0.7;
R_around_scale = 0.6;
threshold_around = 0.4;                 

[final_Neurons, grades, R, around] = IsNeurons_new_4(BW_p, potential_neurons, ...
                    'threshold_angle', threshold_angle, ...
                    'merge_dis', 2, ...
                    'R_max', R_max, 'R_min', R_min, 'annulus', 3, ...
                    'R_around', R_around_scale, 'threshold_around', threshold_around);

draw_circles(final_Neurons, R, image_removed);

% number neurons
for i = 1:length(grades)
    text(final_Neurons(i,1),final_Neurons(i,2),int2str(i),'FontSize',10,'Color','red');
end
%% assign possibility and direction for synapse
BW_thin = bwmorph(BW_p, 'thin', 20);

sigma = 8;
[U, V, theta] = neurite_vector(BW_p, sigma);

theta = theta .* BW_thin;
%% Synapse detection
% degree and distance threshold when connecting broken synapse
theta_thre = 10;
fill_gap = 10;
[connected, synapse] = synapse_detection(BW_p, BW_thin, final_Neurons,...
                                                R, theta, theta_thre, fill_gap);
                          
breadth = synapse_breadth(synapse, BW_p);

draw_synapse(connected, synapse, breadth);
% %% Astar for the shortest paths
% num = size(connected, 1);
% synapse_new = cell(num);
% 
% connected_tri = triu(connected, 0);
% 
% figure; imshow(BW_thin); hold on; 
% 
% for i = 1 : 1
%     synapse_i = synapse_new(i, :);
%     for j = i + 1 : num
%         if connected(i, j)
%             path = synapse{i, j};
%             Start = path(2, :);
%             optimal_path = Astar_Synapse(BW_p, i, R, final_Neurons, j, theta, ...
%                                     fill_gap, theta_thre, connected_tri(i, :));
% %             optimal_path = path_Astar(BW_thin, Start, R, final_Neurons, ...
% %                                 Target, theta, fill_gap, theta_thre);
%             synapse_i{j} = optimal_path;
%             synapse_new(i, :) = synapse_i;
%         end
%     end
% end
% 
% 
% draw_circles(final_Neurons, R, image_removed);
% 
% % number neurons
% for i = 1:length(grades)
%     text(final_Neurons(i,1),final_Neurons(i,2),int2str(i),'FontSize',10,'Color','red');
% end
% 
% draw_synapse(connected, synapse_new, breadth);