%% label samples


num = length(samples);

label_neurons = cell(num, 1);
imgs = cell(num, 1);
for i = 1:num
    sample = samples{i};
    label_neurons{i} = readPoints(sample);
end


%% select
% sample without neurons are removed

num = length(label_neurons);

index = randperm(num);
train_num = ceil(num / 2);
train_imgs = imgs(index(1:train_num));
train_labels = label_neurons(index(1:train_num));

test_num = num - train_num;
test_imgs = imgs(index(train_num+1 : num));
test_labels = label_neurons(index(train_num+1 : num));

save('train_imgs.mat', 'train_imgs');
save('train_labels.mat', 'train_labels');
save('test_imgs.mat', 'test_imgs');
save('test_labels.mat', 'test_labels');
%%
for i = 1:9
    draw_neurons(samples{i}, label_neurons{i});
end
%% histogram adjust
function img_new = hist_adjust(img)
    [m,n] = size(img);
    % for the grayscale of main parts in the image is within a small range,
    % adjust to improve contrast.
    hist = imhist(img, 2^16);
    left_removal = 0.0005;
    flag = true;
    right_removal = 0.9995;
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