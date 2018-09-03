function draw_synapse(connected, synapse, image)

num = length(connected);
figure;
imshow(image);
hold on;
title('synapse');

for i = 1 : num - 1
    for j = i + 1 : num
        if connected(i, j)
            len = size(synapse{i, j}, 1);
            path = synapse{i, j};
            for k = 1:len-1
                plot([path(k, 1), path(k+1, 1)],[path(k, 2), path(k+1, 2)],'LineWidth',2,'Color','blue');
            end
        end
    end
end
            