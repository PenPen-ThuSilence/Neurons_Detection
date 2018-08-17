function new_pts = modify_Points(image, points)
%readPoints   Read manually-defined points from image
%   POINTS = READPOINTS(IMAGE) displays the image in the current figure,
%   then records the position of each click of button 1 of the mouse in the
%   figure, and stops when another button is clicked. The track of points
%   is drawn as it goes along. The result is a 2 x NPOINTS matrix; each
%   column is [X; Y] for one point.
% 
new_pts = points;

imshow(image, []); hold on;     % display image
plot(new_pts(:,1), new_pts(:,2), 'ro');
k = size(new_pts, 1);

delete_thre = 50;

while 1
    [xi, yi, but] = ginput(1);      % get a point
    if ~isequal(but, 1)             % stop if not button 1
        break
    end
    dis = sqrt(sum((repmat([xi, yi], k ,1) - new_pts) .^ 2, 2));
    [min_dis, index] = min(dis);
    % delete pt
    if min_dis < delete_thre
        k = k - 1;
        new_pts(index, :) = [];
        imshow(image, []); hold on;
        plot(new_pts(:,1), new_pts(:,2), 'ro');
    else
        k = k + 1;
        new_pts(k,1) = xi;
        new_pts(k,2) = yi;
        imshow(image, []); hold on;
        plot(new_pts(:,1), new_pts(:,2), 'ro');
    end
end
hold off;
end