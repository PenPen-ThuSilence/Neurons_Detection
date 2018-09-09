%% histogram adjust
function img_new = hist_adjust(img, level_low, level_high)
    [m,n] = size(img);
    % for the grayscale of main parts in the image is within a small range,
    % adjust to improve contrast.
    hist = imhist(img, 2^16);
    flag = true;
    counts = m * n;
    for i = 1:2^16
        if sum(hist(1:i)) > level_low * counts && flag
            left = i;
            flag = false;
        elseif sum(hist(1:i)) > level_high * counts
            right = i;
            break;
        end
    end
    img_new = imadjust(img, [left / 65536, right / 65536]);
end