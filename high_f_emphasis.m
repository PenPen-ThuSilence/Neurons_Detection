function fhfebtw = high_f_emphasis(Img, d0, a, b)
%% Input
% d0: radius
% H_hfe = a + b * H_up, a = 0.5, b = 2

f = double(Img);
[r,c] = size(f);

% FFT
F = fft2(f);
% shift
G = fftshift(F);
% center
mu = floor(r/2);
mv = floor(c/2);

n = 2;

[u, v] = meshgrid(1:r, 1:c);
d = sqrt((u - mu).^2 + (v - mv).^2);
Hlpbtw = 1 ./ (1 + 0.414 * (d ./ d0).^(2 * n));
Hhpbtw = 1 - Hlpbtw;
Hhfebtw = a + b * Hhpbtw;

Ghfebtw = G .* Hhfebtw;

% IFFT
ghfebtw = ifftshift(Ghfebtw);

% real part
fhfebtw = uint16(real(ifft2(ghfebtw)));
