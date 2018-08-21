function F1 = detect_perform_test(paras)
                
threshold_binary = paras(1);
threshold_angle = paras(2);
R_center = paras(3);
R_range = paras(4);
R_around = paras(5);
threshold_around = paras(6);
density_thre = paras(7);

imgs = load('samples.mat');
imgs = imgs.samples;
neuron_labels = load('label_neurons.mat');
neuron_labels = neuron_labels.label_neurons;

% evaluate the performance of neuron detection
num_imgs = 1;

TP = 0;
num_neurons_all = 0;
num_labels_all = 0;

% dis threshold for match neruons
match_threshold = 50;

for i = 1:num_imgs
    img = imgs{i};
    labels = neuron_labels{i};
    Neurons = Neuron_detection_new(img, threshold_binary, ... 
                        threshold_angle, R_center, R_range, R_around, ...
                        threshold_around, density_thre);
%     hold on;
%     axis on, xlabel x, ylabel y;
%     plot(labels(:,1),labels(:,2),'.','color','green', 'MarkerSize', 15);                
    
    num_labels = size(labels, 1);
    num_neurons = size(Neurons, 1);
    
    if isempty(num_neurons)
        continue;
    end
    
    num_neurons_all = num_neurons_all + num_neurons;
    num_labels_all = num_labels + num_labels_all;
    
    for j = 1:num_labels
        dis = repmat(labels(j, :), num_neurons, 1) - Neurons;
        dis = sqrt(sum(dis.^2, 2));
        if min(dis) < match_threshold
            TP = TP + 1;
        end
    end
end

% precision = TP / (TP + FP)
precision = TP / num_neurons_all;
% recall = TP / (TP + FN)
recall = TP / num_labels_all;

% F measure to evaluate precision and recall together
F1 = - 2 * (precision * recall) / (precision + recall);

end