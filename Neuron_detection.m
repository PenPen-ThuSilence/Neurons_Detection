%% read image
% Input: filename of image
filename = 'images/cropped-t0.tif';
img = imread(filename);
img = img(5001:7000, 5001:7000);
%% Preprocess
scale = 0.5;
level_low = 0.3;
level_high = 0.98;
disk_size = 1;
label_thre = 20;
intensity_thre = 180;
extent_thre = 0.4;
area_thre = 400;
binary_thre = 35;
min_pixel_num = 40;

[image_removed, BW_p] = preprocess(img, scale, level_low, level_high, ...
                                     disk_size, label_thre, intensity_thre, ...
                                     extent_thre, area_thre, binary_thre, min_pixel_num);
% show                                 
figure;
imshow(image_removed);
figure;
imshow(BW_p);
%% get potential points which can be neurons with kernel
% parameters
R_min = 10;
R_max = 25;
density_thre = 0.45:0.05:1;

potential_neurons = intensive_points(BW_p, R_min, R_max, density_thre);
draw_neurons(BW_p, potential_neurons);
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

% parameters
threshold_angle = 0.65;
R_max = 22;
R_min = 10;
R_around = 6;
threshold_around = 0.2;

[final_Neurons, grades, R, around] = IsNeurons_new_4(BW_p, potential_neurons, ...
                    'threshold_angle', threshold_angle, ...
                    'merge_dis', 2, ...
                    'R_max', R_max, 'R_min', R_min, 'annulus', 3, ...
                    'R_around', R_around, 'threshold_around', threshold_around);

draw_circles(final_Neurons, R, image_removed);

% number neurons
for i = 1:length(grades)
    text(final_Neurons(i,1),final_Neurons(i,2),int2str(i),'FontSize',10,'Color','red');
end

%% Synapse detection
% degree and distance threshold when connecting broken synapse
theta_thre = 15;
fill_gap = 10;
[connected, synapse, breadth] = synapse_detection(BW_p, final_Neurons, R,...
                                                theta_thre, fill_gap);
                                            
draw_synapse(connected, synapse, image_removed, breadth);