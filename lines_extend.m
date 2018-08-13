function lines = lines_extend(lines, ex_l)

% number of lines
N = length(lines);

parfor k = 1:N
   xy = [lines(k).point1; lines(k).point2];
   if lines(k).theta == 0
       [~,low] = min(xy(:,2));
       xy(low, 2) = xy(low, 2) - ex_l;
       high = 3 - low;
       xy(high, 2) = xy(high, 2) + ex_l;
   else
       slope = (xy(1,2) - xy(2,2)) / (xy(1,1) - xy(2,1));
       ex_lx = ex_l / sqrt(1 + slope^2);
       [~,left] = min(xy(:,1));
       right = 3 - left;
       xy(left, :) = xy(left, :) - [ex_lx, ex_lx * slope];
       xy(right, :) = xy(right, :) + [ex_lx, ex_lx * slope];
   end
   lines(k).point1 = xy(1,:);
   lines(k).point2 = xy(2,:);
end
