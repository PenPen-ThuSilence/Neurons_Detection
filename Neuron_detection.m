%% detect centers of neurons
function [Neurons, R, BW_p] = Neuron_detection(img, d0, num_peaks,...
                        min_line_length, fill_gap, extend_length, ...
                        threshold_angle, R_center, R_range, R_around, ...
                        threshold_around)
%% Input Parameters
% img: image to process
% d0: distance in butterworth filter, reference value: 12
% num_peaks: number of peaks in hough TF, reference value: 8000 (img: 4000 x 4000)
% min_line_length: minimum length of lines detected with hough TF,
% reference value: 55
% fill_gap: maximum length of gaps to fill when detecting lines, reference
% value: 10
% extend_length: length to extend lines, reference value: 180
% threshold_angle: threshold of angles when judging neurons, reference
% value: 0.65
% R_center, R_range: range of radiuses of neurons to detect, reference
% value: R_center(75), R_range(25), which means radius range is 50:100
% R_around: soma radius when judging neurons, reference value: 20
% threshold_around: threshold of around

% example:
% [Neurons, R, BW_p] = Neuron_detection(img, 12, 8000,...
%                        55, 10, 180, ...
%                        0.65, 75, 25, 20, ...
%                        0.25)
%% 1. Read the image
% read image and convert it to uint16
img = uint16(img);
[m,n] = size(img);

%% 2. high frequency emphasis filtering
% use butterworth filter to emphasize parts with high frequency
a = 0.5;
b = 4;
% fhfebtw = a * img_new + b * high-frequency parts
fhfebtw = high_f_emphasis(img, d0, a, b);

% figure; imshow(fhfebtw);

%% 3. remove background
% thresholding
img_remove_back = fhfebtw > 15000;
% erode and then dilate to remove remaining small dots in the background
se = strel('disk', 1);
img_remove_back = imopen(img_remove_back, se);
% figure; imshow(img_remove_back);

%% 4. detect white dots
% dots and synapses are emphasized after filtering
dots = fhfebtw > 65000;

% exclude synapses
minL = 20;
filtered = filterRegions_MajorAxis(dots, minL);
dots_filtered = dots & ~filtered;

%% 5. remove dots

% erode dots, because detected dots are smaller than original ones
se = strel('disk', 4);
dots_dialte = imdilate(dots_filtered, se);
BW = img_remove_back & ~dots_dialte;

% fill holes, because it can also remove something inside neurons
BW_filled = imfill(BW, 'holes');
fill_differ = BW_filled & ~BW;
% remove some filling of very large area
min_area = 300;
fill = fill_differ & ~filterRegions_area(fill_differ, min_area);

BW = BW | fill;

% dilate and erode to link structures of neurons
BW = imclose(BW, se);

% region area to remove remaining white dots, we fill holes and link
% structures of neurons to protect neurons here.
min_area = 200;
BW_p = filterRegions_area(BW, min_area);
% here we finish preprocessing parts
%% 6. synapse detection
% We can extend synapses to intersect within neurons, so first we need to
% detect lines to get synapses.

% Hough TF: (x,y) -> (rho, theta). 

% Based on hough tf, we can detect lines in the image which are mostly
% synapses. Before transfroms, we should do some preprocessing: 

% 1. Get edge with canny method. 
% 2. The edges of a synapse are obviously two parallel lines, here we use close compution to combine the two parallel lines. 
% 3. Then we use bwmorph to thin the edge to single-bit branches.
% 4. BUT single-bit branches are hard to detect, so we do another dilation.

lines = line_detect_BW(BW_p, round(num_peaks), ...
    'theta_space', 3, 'rho_space', 10, 'fill_gap', fill_gap, 'min_length', min_line_length);
% lines = line_detection(BW_p, num_peaks, 'dilate_size', 2, ...
%     'theta_space', 2, 'rho_space', 5, 'fill_gap', fill_gap, 'min_length', min_line_length);
% parameters here
% 
% 1. num_peak: how many peaks to use after Hough TF
% 2. dilate_size: (in step 4) dilation makes it easier to detect lines.
% 3. theta_space&rho_space: parameters that influence the range of TF
%    results, the fewer space, the more lines.
% 4. fill_gap: length limit of gap to connect two lines
% 5. min_length: length limit of lines(very short lines are ignored)

%% 7. extend lines to intersection
% extend_length: extended length in both directions
% short to miss intersections, long to complex

extended_lines = lines_extend(lines, extend_length);

% get intersection points of extended lines
points = line_intersect_points(extended_lines, size(BW_p,1), size(BW_p,2));

% draw_lines(BW_p, extended_lines); hold on;
% axis on, xlabel x, ylabel y;
% plot(points(:,1),points(:,2),'.','color','blue', 'MarkerSize', 8);
% title('intersections');
%% 6. get neurons from intersection points
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

[Neurons, grades, R, around] = IsNeurons_new_2(BW_p, points, ...
                    'threshold_angle', threshold_angle, ...
                    'merge_dis', 4, ...
                    'R', R_center, 'R_range', R_range, 'annulus', 3, ...
                    'R_around', R_around, 'threshold_around', threshold_around);

draw_circles(Neurons, R, BW_p);

% number neurons
for i = 1:length(grades)
    text(Neurons(i,1),Neurons(i,2),int2str(i),'FontSize',15,'Color','red');
end
% the method still needs improvements, because there are obvious 
% mistakes from the image. We need to add filters to avoid mistaking
% some cross of synapses as neurons and try to get more neurons.
% In addition, it's very time-consuming but it can be parallel.

%% 7. find synapse from neurons
% With neurons and R we get before, we use a simple method to find synapse.

% First, we initialize queues of each neurons with intersections of
% circle(neuron, R) and bright points in BW. 

% And we label points in queues of different neurons with 
% different labels. 

% Then for each queue, we get queue head out, 
% find its bright neighbours(conn = 8) and add them to queue.

% If there are labeled neighbours of different neurons,
% we get a synapse. 

% If queues are all empty, the algorithm terminates.

% For we start expanding from every neuron, 
% different labels meet and get synapse.

% [connectivity, label] = queue_synapse(Neurons, roi, R);

% print matrix of connections

% In the preprocessing, we may break off some synapses,  
% so we need to try to link broken synapses. 
% Besides, the queue-method is quite time-consuming 
% if applied to larger images because it can't be parallel

%% 8. link cracked synapses

% link_dis = 10;
% linked_connectivity = link(Neurons, label, link_dis, connectivity);