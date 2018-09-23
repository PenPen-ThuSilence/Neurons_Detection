function [Neurons, grades, R, around] = IsNeurons_new_4(I, points, varargin)
%% Input, Output, Init
% Input definition
ip = inputParser;

% Required
ip.addRequired('I');
ip.addRequired('points');

% Parameter
% Iterations to judge that certain point is not a neuron.
% ip.addParameter('threshold', 0.8);
ip.addParameter('threshold_angle', 0.65);
% merge dis
ip.addParameter('merge_dis', 3);
% R is adjustable in the range
ip.addParameter('R_min', 25);
ip.addParameter('R_max', 50);
% max move of points in each iteration
% ip.addParameter('max_move', 1);
% annulus: width of the circle
ip.addParameter('annulus', 3);
ip.addParameter('R_around', 15);
ip.addParameter('threshold_around', 0.3);

% Parse input
parse(ip, I, points, varargin{:});
ir = ip.Results;

% Output
Neurons = zeros(0,2);
R = [];
% grade of each neuron
grades = [];
around = [];

% size of I
[M,N] = size(I);

% BW
if length(unique(I)) <= 2
    BW = I;
else
    BW = (I > 2^11);
end

R_p = ir.R_min:ir.R_max;
R_p = round(R_p);

ir.R_around = round(ir.R_around);

    %% merge near points
% %     Neurons = merge_neurons(Neurons, ir.merge_dis);
%     for ii = 1:num-1
%         for jj = ii+1:num
%             dis = (points(ii, :) - points(jj, :)).^2;
%             dis = sqrt(sum(dis));
%             if dis < ir.merge_dis
%                 points(ii,:) = mean(points([ii,jj], :));
%                 points(jj,:) = NaN;
%             end
%         end
%     end
% % %     fprintf('Merge points: %d ---> %d, ', size(points,1), size(points,1) - sum(isnan(points(:,1))));
%     points(isnan(points(:,1)), :) = [];
%% judge points
points = round(points);
parfor k = 1:size(points,1)
    point_k = points(k,:);
    % get random R
    len = length(R_p);
    angles = zeros(len, 1);
    for j = 1:len
        % generate meshgrid
        % remove edge
        annulus = ir.annulus;
        row = max(1, point_k(1) - R_p(j) - annulus):...
            min(N, point_k(1) + R_p(j) + annulus);
        col = max(1, point_k(2) - R_p(j) - annulus)...
            :min(M, point_k(2) + R_p(j) + annulus);
        [area_x, area_y] = meshgrid(row, col);
        % calculate distance
        dist = sqrt((area_x - point_k(1)).^2 + (area_y - point_k(2)).^2);
        % number of points in the circular area of points(k)
        in_circle = dist >= R_p(j) - annulus & dist <= R_p(j) + annulus;
        % number of points in the area that are bright
        in_circle_bright = in_circle & BW(col, row);
        % angle that bright points appear
        slope = atan2d(area_y - point_k(2), area_x - point_k(1));
        slope_bright = slope(in_circle_bright);
        angle = length(unique(round(slope_bright)));
        angles(j) = angle / length(unique(round(slope)));
    end
    [grade, index] = max(angles);
    best_R = R_p(index);
    R_soma = ir.R_around * best_R;
    % whether exists soma around 
    row = max(1, point_k(1) - R_soma):...
        min(N, point_k(1) + R_soma);
    col = max(1, point_k(2) - R_soma)...
        :min(M, point_k(2) + R_soma);
    [area_x, area_y] = meshgrid(row, col);
    dist = sqrt((area_x - point_k(1)).^2 + (area_y - point_k(2)).^2);
    circle_around = dist < R_soma;
    circle_around_bright = circle_around & BW(col, row);
    around_p = sum(circle_around_bright(:)) / sum(circle_around(:));
    if max(angles) > ir.threshold_angle && around_p > ir.threshold_around
        Neurons = [Neurons; points(k,:)];
        grades = [grades, grade];
        R = [R, best_R];
        around = [around, around_p];
    end         
end

% merge near points
temp = 0;
merge_dis = ir.R_max * 1.25;
while true
    num = size(Neurons, 1);
    if temp == num
        break;
    end
    temp = num;
    for i = 1:num-1
        for j = i+1:num
            dis = (Neurons(i, :) - Neurons(j, :)).^2;
            dis = sqrt(sum(dis));
            if dis < merge_dis
                [~, choose] = max([grades(i), grades(j)]);
                if choose == 1
                    Neurons(j,:) = NaN;
                    grades(j) = NaN;
                    R(j) = NaN;
                    around(j) = NaN;
                else
                    Neurons(i,:) = NaN;
                    grades(i) = NaN;
                    R(i) = NaN;
                    around(i) = NaN;
                end
            end
        end
    end
    Neurons(isnan(Neurons(:,1)), :) = [];
    grades(isnan(grades)) = [];
    R(isnan(R)) = [];
    around(isnan(around)) = [];
end