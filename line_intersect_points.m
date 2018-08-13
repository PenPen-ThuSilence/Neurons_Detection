%% get points where lines intersect
function points = line_intersect_points(lines, m, n)
% m, n define the bounds
% get number of line segment
num_l = length(lines);
points = [];
% get intersections of every two line segments
for i = 1 : num_l-1
    for j = i+1 : num_l
        if lines(i).theta == 0
            k1 = Inf;
            b1 = lines(i).point1(1);
        else
            if lines(i).theta == -90
                k1 = 0;
            else
                k1 = (lines(i).point1(2) - lines(i).point2(2))/...
                    (lines(i).point1(1) - lines(i).point2(1));
            end
            b1 = lines(i).point1(2) - k1 * lines(i).point1(1);
        end
        if lines(j).theta == 0
            k2 = Inf;
            b2 = lines(j).point1(1);
        else
            if lines(j).theta == -90
                k2 = 0;
            else
                k2 = (lines(j).point1(2) - lines(j).point2(2))/...
                    (lines(j).point1(1) - lines(j).point2(1));
            end
            b2 = lines(j).point1(2) - k2 * lines(j).point1(1);
        end
        p = round(intersection(k1,b1,k2,b2));
        if (lines(i).point1(1) - p(1)) * (lines(i).point2(1) - p(1)) <= 0 ...
                && (lines(j).point1(1) - p(1)) * (lines(j).point2(1) - p(1)) <= 0 ...
                && p(1) > 0 && p(1) <= n && p(2) > 0 && p(2) <= m

            points = [points; p];
        end
    end
end