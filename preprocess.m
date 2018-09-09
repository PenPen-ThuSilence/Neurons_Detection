function [image_removed, BW_p] = preprocess(img, scale, level_low, level_high, ...
                                     disk_size, label_thre, intensity_thre, ...
                                     extent_thre, area_thre, binary_thre, min_pixel_num)
%% img_as_uint16
img = double(img);
img = (img - min(img(:))) / range(img(:)) * 65535;
img = uint16(img);

%% resize to reduce resolution
img_resize = imresize(img, scale, 'bilinear');

%% adjustment of histogram

image_adjusted = hist_adjust(img_resize, level_low, level_high);
image_adjusted = uint8(image_adjusted / 256);
%% enhance contrast
se = strel('disk', disk_size);
num = sum(se.Neighborhood(:));
% maximum in disk domain
max_filter_image = ordfilt2(image_adjusted, num, se.Neighborhood);
% minimum in disk domain
min_filter_image = ordfilt2(image_adjusted, 1, se.Neighborhood);
% enhance 
max_differ = max_filter_image - image_adjusted;
min_differ = image_adjusted - min_filter_image;
% if near maximum, assign pixel as maximum of domin. 
% The same goes with minimum.
image_enhanced = uint8(max_differ > min_differ) .* min_filter_image + ...
                 uint8(max_differ <= min_differ) .* max_filter_image;
%% remove unconnected neurons

% binarize and label regions
BW = image_enhanced > label_thre;
% label of regions
[L, num] = bwlabel(BW);
% regionprops
STATS = regionprops(L, 'Extent', 'Area');
STATS = cell2mat(struct2cell(STATS));
% median of regions
p_median = zeros(1, num);
p_index = regionprops(L, 'PixelIdxList');
for i = 1:num
    index_i = p_index(i).PixelIdxList;
    p_median(i) = median(image_enhanced(index_i));
end
% area of regions
p_area = STATS(1, :);
% extent of regions
p_extent = STATS(2, :);

image_removed = image_enhanced;
% remove unconnected neurons
remove_index = find(p_median > intensity_thre & p_area < area_thre ...
                    & p_extent > extent_thre);
                
image_removed(ismember(L, remove_index)) = 0;
%% preprocess for neuron detection

% binarize
BW = image_removed > binary_thre;

% dilate and erode to connect parts of neurons
se = strel('disk', 2);
BW_p = imclose(BW, se);

% fill holes in neurons which caused when removing dots
BW_filled = imfill(BW_p, 'holes');
fill_differ = BW_filled & ~BW_p;
% preventing from fill very large holes
min_area = 40;
fill = fill_differ & ~filterRegions_area(fill_differ, min_area);

BW_p = BW_p | fill;

% remove small regions in binary images which were dark in uint8 images
BW_open = bwareaopen(BW_p, min_pixel_num);

BW_p = BW_open;