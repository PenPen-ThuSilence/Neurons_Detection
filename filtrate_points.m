%% filtrate points

index = sub2ind(size(image), points(:,2), points(:,1));

BW = zeros(size(image));
BW(index) = 1;

r = 90;

kernel = circle_kernel(r, 0);

density = conv2(BW, kernel, 'same');

max_density = ordfilt2(density, sum(kernel(:)),kernel);

max_count = max(density(:));

imshow(density > 0.1 * max_count);

%%
R_center = 75;
R_range = 35;
Neurons = zeros(0, 2);

for r = R_center - R_range:R_center + R_range
    kernel = circle_kernel(r, 0);

    density = conv2(double(BW_p), kernel, 'same');

    neu = density > 0.5 * sum(kernel(:));
    
    stats = regionprops(neu, 'Centroid');
    len = length(stats);
    Centroids = zeros(len, 2);
    for i = 1:len
        Centroids(i, :) = stats(i).Centroid;
    end
    Neurons = [Neurons; Centroids]; 
end

draw_neurons(img, Neurons);
%%
imshowpair(neu, BW_p);
stats = regionprops(neu, 'Centroid');