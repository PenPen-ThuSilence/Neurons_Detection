% produce samples

img_0 = imread('images/preprocessed-t20.tif');

x_start = 2000;
x_end = 7000;
y_start = 4000;
y_end = 11000;

img = img_0(x_start:x_end, y_start:y_end);

region_x = 2500;
region_y = 2500;

sample_num = 10;

samples = cell(sample_num, 1);

for i = 1:sample_num
   x = round(rand(1) * (size(img, 1) - region_x - 1));
   y = round(rand(1) * (size(img, 2) - region_y - 1));
   sample = img(x : x + region_x - 1, y : y + region_y - 1);
   figure;
   imshow(sample);
   samples{i} = sample;
end

%%
for i = 1:length(samples)
    figure;
    imshow(samples{i});
end