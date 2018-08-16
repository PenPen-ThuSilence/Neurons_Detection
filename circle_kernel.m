function kernel = circle_kernel(Rmax, Rmin)

% init rows and columns
[row, col] = deal(-Rmax:Rmax);
% meshgrid
[x, y] = meshgrid(row, col);
% distance from center (0,0)
dist = x.^2 + y.^2;
% kernel
kernel = double(dist <= Rmax^2 & dist >= Rmin^2);

% figure;
% imshow(kernel);
% title(['Circular kernel with R = ' int2str(R)]);