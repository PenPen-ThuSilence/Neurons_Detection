function circle_area = draw_circles_k(neurons, R, I, k)

figure;
imshow(I), hold on;
axis on, xlabel x, ylabel y;
plot(neurons(:,1),neurons(:,2),'.','color','red', 'MarkerSize', 15);

[M, N] = size(I);
num = length(R);
circle_area = cell(num, 1);

for i = 1:num
    % generate meshgrid
    annulus = 1;
    row = max(1, neurons(i,1) - R(i) - annulus):...
        min(N, neurons(i,1) + R(i) + annulus);
    col = max(1, neurons(i,2) - R(i) - annulus)...
        :min(M, neurons(i,2) + R(i) + annulus);
    [area_x, area_y] = meshgrid(row, col);

    dist = sqrt((area_x - neurons(i,1)).^2 + (area_y - neurons(i,2)).^2);
    % number of points in the circular area of points(k)
    in_circle = dist >= R(i) - annulus & dist <= R(i) + annulus;
    if i == k
        plot(area_x(in_circle), area_y(in_circle), '.','color','green', 'MarkerSize', 1);
    else
        plot(area_x(in_circle), area_y(in_circle), '.','color','red', 'MarkerSize', 1);
    end
    circle_area{i} = [area_x(in_circle), area_y(in_circle)];
end