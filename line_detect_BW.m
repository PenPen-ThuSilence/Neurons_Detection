function lines = line_detect_BW(BW, num_peaks, varargin)
%% Input Definition
    ip = inputParser;
    ip.addRequired('BW');
    ip.addRequired('num_peaks');
    ip.addParameter('theta_space', 1);
    ip.addParameter('rho_space', 5);
    ip.addParameter('fill_gap', 5);
    ip.addParameter('min_length',60);
    parse(ip, BW, num_peaks, varargin{:});
    ir = ip.Results;
%% find Peaks of Hough TF
    [H,theta,rho] = hough(BW, 'RhoResolution', ir.rho_space,...
        'Theta', -90:ir.theta_space:...
        90-mod(180, ir.theta_space)-ir.theta_space * (mod(180, ir.theta_space)==0));
    
    % Peaks of hough TF
    P = houghpeaks(H, num_peaks, 'threshold', 0);
    
%     figure;
%     imshow(imadjust(rescale(H)),[],...
%            'XData',theta,...
%            'YData',rho,...
%            'InitialMagnification','fit');
%     xlabel('\theta (degrees)');
%     ylabel('\rho');
%     axis on;
%     axis normal; 
%     hold on;
%     colormap(gca,hot);
%     x = theta(P(:,2));
%     y = rho(P(:,1));
%     plot(x,y,'.','color','green');
    
%% Houghlines
    lines = houghlines(BW, theta, rho, P,...
        'FillGap', ir.fill_gap, 'MinLength', ir.min_length);
    
end