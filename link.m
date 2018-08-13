function connectivity = link(Neurons, label, link_dis, connectivity)

% link_dis: distance limit to link terminal vertexes

% draw dentrites
num = size(Neurons,1);
colors = rand(1,3,num);
imshow(label); hold on;
axis on;
for i = 1:num
    [points_i_x, points_i_y] = find(label == i + 1);
    plot(points_i_y, points_i_x, '.', 'color', colors(:,:,i));
    text(Neurons(i,1),Neurons(i,2),int2str(i),'FontSize',25,'Color','red');
end

BW = label > 1;
hold on;
hulls = regionprops(BW, 'ConvexHull');

for i = 1:length(hulls)
    plot(hulls(i).ConvexHull(:,1), hulls(i).ConvexHull(:,2), '.', 'color', rand(1,3));
end

% bound of image
[M,N] = size(label);

for i = 1:length(hulls)
    for j = 1:size(hulls(i).ConvexHull, 1)
        % from a vertex
        vertex = floor(hulls(i).ConvexHull(j, :));
        vertex = vertex + (vertex == 0);
        v_label = label(vertex(2), vertex(1)) - 1;
        if (v_label <= 0)
            continue;
        end
        % find all points whose distance from vertex < link_dis
        row = max(1, vertex(1) - link_dis):...
            min(N, vertex(1) + link_dis);
        col = max(1, vertex(2) - link_dis)...
            :min(M, vertex(2) + link_dis);
        [area_x, area_y] = meshgrid(row, col);
        % distance matrix
        dist = sqrt((area_x - vertex(1)).^2 + (area_y - vertex(2)).^2);
        % number of points in the circular area of points(k)
        in_circle = dist <= link_dis;
        % 1-D coordinates of pixels in the circle
        x_in_circle = area_x(in_circle);
        y_in_circle = area_y(in_circle);
        index = sub2ind([M,N], y_in_circle, x_in_circle);
        % all different labels in the circle area
        [labels_in_circle, IA] = unique(label(index));
        % ignore 0 (background) 1 (unlabeled bright pixels)
        labels_in_circle = labels_in_circle - 1;
        for k = 1:length(labels_in_circle)
            if labels_in_circle(k) < 1
                continue;
            elseif ~connectivity(v_label, labels_in_circle(k))
                connectivity(v_label, labels_in_circle(k)) = true;
                connectivity(labels_in_circle(k), v_label) = true;
                plot([vertex(1), x_in_circle(IA(k))], ...
                    [vertex(2), y_in_circle(IA(k))],'r','LineWidth', 5);
            end
        end
    end
end