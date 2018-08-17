%% main function
load('train_imgs.mat');
load('train_labels.mat');

d0 = 10;
num_peaks = 6000;
min_line_length = 80;
fill_gap = 12;
extend_length = 150;
threshold_angle = 0.7;
R_center = 75;
R_range = 30;
R_around = 40;
threshold_around = 0.3;

F1 = detect_perform(train_imgs, train_labels, d0, num_peaks,...
                        min_line_length, fill_gap, extend_length, ...
                        threshold_angle, R_center, R_range, R_around, ...
                        threshold_around);