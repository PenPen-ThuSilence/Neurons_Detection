function points = intensive_points(BW_p, R_min, R_max, density_thre)

points = zeros(0, 2);

% radius range
r_range = R_min : 2 : R_max;
num = length(r_range);

for j = 1:num
    r = r_range(j);
    % circular kernel: pixels inside circle(r) are 1, outside ones are 0.
    kernel = circle_kernel(r, 0);
    
    % conv2
    density = conv2(double(BW_p), kernel, 'same');
    
    % pixels where conv2 result are larger, which means bright areas
    for i = 1:length(density_thre)
        neu = density > density_thre(i) * sum(kernel(:));

        % find centroids of these bright areas
        stats = regionprops(neu, 'Centroid');
        len = length(stats);
        Centroids = zeros(len, 2);
        for ii = 1:len
            Centroids(ii, :) = stats(ii).Centroid;
        end
        % centroids of bright areas are potential neurons
        points = [points; Centroids]; 
    end
end