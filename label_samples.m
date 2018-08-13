%% label samples

num = 24;
label_neurons = cell(num, 1);
imgs = cell(num, 1);
for i = 0:num-1
    filename = ['samples/sample_', int2str(i), '.tif'];
    img = imread(filename);
    img = (img - min(img(:))) / (max(img(:)) - min(img(:))) * 65535;
    img = uint16(img);
    img = hist_adjust(img);
    imgs{i+1} = img;
%     sample = imread(filename);
%     sample = hist_adjust(sample);
%     label_neurons{i+1} = readPoints(sample);
%     save label_neurons;
end


%% select
% sample without neurons are removed
load labels;
% 14, 17, 19
label_neurons(cellfun(@isempty,label_neurons)) = [];
imgs([14,17,19]) = [];

num = length(label_neurons);

index = randperm(num);
train_num = ceil(num / 2);
train_imgs = imgs(index(1:train_num));
train_labels = label_neurons(index(1:train_num));

test_num = num - train_num;
test_imgs = imgs(index(train_num+1 : num));
test_labels = label_neurons(index(train_num+1 : num));

save(train_imgs.mat, train_imgs);
save('train_labels.mat', train_labels);
save('test_imgs.mat', test_imgs);
save('test_labels.mat', test_labels);
%%
for i = 1:11
    draw_neurons(train_imgs{i}, train_labels{i});
end
%% histogram adjust
function img_new = hist_adjust(img)
    [m,n] = size(img);
    % for the grayscale of main parts in the image is within a small range,
    % adjust to improve contrast.
    hist = imhist(img, 2^16);
    left_removal = 0.001;
    flag = true;
    right_removal = 0.999;
    counts = m * n;
    for i = 1:2^16
        if sum(hist(1:i)) > left_removal * counts && flag
            left = i;
            flag = false;
        end
        if sum(hist(1:i)) > right_removal * counts
            right = i;
            break;
        end
    end
    img_new = imadjust(img, [left / 65536, right / 65536]);
end