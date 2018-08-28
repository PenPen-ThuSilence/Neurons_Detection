%%
image = imread('images/preprocessed-t0.tif');
image = image(4401:4700, 6451:6750);
%%
% same function as img_as_float
image = double(image);
image = image ./ max(image(:));

% gaussian kernel
sigma = 2;

gaussian_img = imgaussfilt(image, sigma, 'padding', 'replicate');

% second-order derivatives
[gx, gy] = gradient(gaussian_img);
[gxx, gxy] = gradient(gx);
[gyx, gyy] = gradient(gy);

[m, n] = size(image);
E = zeros(m, n);
U = E;
V = E;

for i = 1:m
   for j = 1:n 
        Hessian = [gxx(i,j) gxy(i,j); gyx(i,j) gyy(i,j)];
        [eigen_vector, eigen_value] = eig(Hessian, 'vector');
        [~, index] = max(abs(eigen_value));
        lambda_x = eigen_value(index);
        E(i, j) = (lambda_x < 0) * lambda_x;
        % eigenvector
        U(i, j) = eigen_vector(1, 3 - index);
        V(i, j) = eigen_vector(2, 3 - index);
   end
end

% eigenvector .* eigenvalue
U = U .* E;
V = V .* E;

figure; imshow(image);

hold on;
quiver(U, V);

theta = atan2d(V, U);
%%

BW_thin = bwmorph(BW_p, 'thin', 20);