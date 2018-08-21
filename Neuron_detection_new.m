%% detect centers of neurons
function [final_Neurons, R, BW_p] = Neuron_detection_new(image, threshold_binary, ... 
                        threshold_angle, R_center, R_range, R_around, ...
                        threshold_around, density_thre)
%% Input Parameters
% img: image to process
% threshold_angle: threshold of angles when judging neurons, reference
% value: 0.65
% R_center, R_range: range of radiuses of neurons to detect, reference
% value: R_center(75), R_range(25), which means radius range is 50:100
% R_around: soma radius when judging neurons, reference value: 20
% threshold_around: threshold of around

%% preprocess
% binarize
BW_p = image > threshold_binary;

% dilate and erode to connect parts of neurons
se = strel('disk', 4);
BW_p = imclose(BW_p, se);

% fill holes in neurons which caused when removing dots
BW_filled = imfill(BW_p, 'holes');
fill_differ = BW_filled & ~BW_p;
% preventing from fill very large holes
min_area = 200;
fill = fill_differ & ~filterRegions_area(fill_differ, min_area);

BW_p = BW_p | fill;
% figure;imshowpair(BW_p, image > threshold_binary);
% title('fill holes');
%% potential neurons
% radius range
Neurons = zeros(0, 2);

R_center = round(R_center);
R_range = round(R_range);

r_range = R_center - R_range : 3 : R_center + R_range;
num = length(r_range);

parfor j = 1:num
    r = r_range(j);
    % circular kernel: pixels inside circle(r) are 1, outside ones are 0.
    kernel = circle_kernel(r, 0);
    
    % conv2
    density = conv2(double(BW_p), kernel, 'same');
    
    % pixels where conv2 result are larger, which means bright areas
    neu = density > density_thre * sum(kernel(:));
    
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
% draw_neurons(BW_p, Neurons);
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

[final_Neurons, grades, R, around] = IsNeurons_new_4(BW_p, Neurons, ...
                    'threshold_angle', threshold_angle, ...
                    'merge_dis', 2, ...
                    'R', R_center, 'R_range', R_range, 'annulus', 2, ...
                    'R_around', R_around, 'threshold_around', threshold_around);

% draw_circles(final_Neurons, R, image);
% 
% % number neurons
% for i = 1:length(grades)
%     text(final_Neurons(i,1),final_Neurons(i,2),int2str(i),'FontSize',10,'Color','red');
% end