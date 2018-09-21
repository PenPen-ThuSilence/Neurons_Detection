function optimal_path = path_Astar(BW, Start, R, Neurons, Target, thetas, ...
                                    fill_gap, theta_thre)

[MAX_Y, MAX_X] = size(BW);
xStart = Start(1);
yStart = Start(2);
xTarget = Target(1);
yTarget = Target(2);

R_compare = 1.25 * R';
R_compare(source) = 1.15 * R(source);

OPEN = [];

% Put all background points on the Closed list
[close_y, close_x] = find(~BW);

CLOSED = [close_x, close_y];
CLOSED_COUNT = size(CLOSED, 1);
% set the starting node as the first node
xNode = xStart;
yNode = yStart;
OPEN_COUNT = 1;
% path cost = 0
path_cost = 0;
goal_distance = distance(xNode, yNode, xTarget, yTarget);
OPEN(OPEN_COUNT, :) = insert_open(xNode, yNode, xNode, yNode, ...
                                path_cost, goal_distance, goal_distance);
former_theta = thetas(yNode, xNode);
% remove from list
OPEN(OPEN_COUNT, 1) = 0;
CLOSED_COUNT = CLOSED_COUNT + 1;
CLOSED(CLOSED_COUNT, :) = [xNode, yNode];

NoPath = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% START ALGORITHM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
while((xNode ~= xTarget || yNode ~= yTarget) && NoPath == 1)
    exp_array = expand_array(xNode,yNode,path_cost,xTarget,yTarget,CLOSED,...
                                MAX_X,MAX_Y, former_theta, thetas, fill_gap, theta_thre);
    exp_count = size(exp_array,1);
    %UPDATE LIST OPEN WITH THE SUCCESSOR NODES
    %OPEN LIST FORMAT
    %--------------------------------------------------------------------------
    %IS ON LIST 1/0 |X val |Y val |Parent X val |Parent Y val |h(n) |g(n)|f(n)|
    %--------------------------------------------------------------------------
    %EXPANDED ARRAY FORMAT
    %--------------------------------
    %| X val | Y val | h(n) | g(n) | f(n) |
    %--------------------------------
    for i = 1 : exp_count
        % whether expanded point in OPEN list
        exp_point = exp_array(i, :);
        Isopen = ismember(exp_point(1:2), OPEN(:,2:3), 'rows');
        Index_in_open = find(Isopen);
        if ~isempty(Index_in_open)
            index_open = find(ismember(OPEN(:,2:3), exp_point(1:2), 'rows'));
            if exp_point(5) < OPEN(index_open, 8)
                OPEN(index_open, 4:8) = [xNode, yNode, exp_point(3:end)];
            end
        else
            OPEN_COUNT = OPEN_COUNT + 1;
            OPEN(OPEN_COUNT, :) = insert_open(exp_point(1), exp_point(2), ...
                                              xNode, yNode, exp_point(3), ... 
                                              exp_point(4), exp_point(5));
        end%End of insert new element into the OPEN list
    end%End of i for
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%END OF WHILE LOOP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Find out the node with the smallest fn 
    index_min_node = min_fn(OPEN, OPEN_COUNT, xTarget, yTarget);
    if (index_min_node ~= -1)    
    %Set xNode and yNode to the node with minimum fn
        xNode = OPEN(index_min_node, 2);
        yNode = OPEN(index_min_node, 3);
        former_theta = atan2d(yNode - OPEN(index_min_node, 5), ...
                        xNode - OPEN(index_min_node, 4));
        path_cost = OPEN(index_min_node, 6);%Update the cost of reaching the parent node
        %Move the Node to list CLOSED
        CLOSED_COUNT = CLOSED_COUNT + 1;
        CLOSED(CLOSED_COUNT, :) = [xNode, yNode];
        OPEN(index_min_node, 1) = 0;
        plot(xNode,yNode,'bo');
    else
        %No path exists to the Target!!
        NoPath = 0;%Exits the loop!
    end%End of index_min_node check
end%End of While Loop

%Once algorithm has run The optimal path is generated by starting of at the
%last node(if it is the target node) and then identifying its parent node
%until it reaches the start node.This is the optimal path

terminal = CLOSED(end, :);

optimal_path = terminal;
path_len = 1;

if terminal == Target % find path from start to target
node_index = find(ismember(OPEN(:, 2:3), terminal, 'rows'));
% parent node of terminal
parent = OPEN(node_index, 4:5);
% path
while parent ~= Start
    % add current point
    path_len = path_len + 1;
    optimal_path(path_len, :) = parent;
    % find parent of current point
    node_index = find(ismember(OPEN(:, 2:3), parent, 'rows'));
    parent = OPEN(node_index, 4:5);
end
%Plot the Optimal Path!
%     p = plot(optimal_path(end, 1)+.5,optimal_path(end, 2)+.5,'bo');
%     for i = path_len - 1 : -1 : 1
%         pause(.25);
%         set(p, 'XData', optimal_path(i,1) + .5, 'YData', optimal_path(i,2) + .5);
%         drawnow;
%     end
%     plot(optimal_path(:,1) + .5, optimal_path(:, 2) + .5);
else
    pause(1);
    h=msgbox('Sorry, No path exists to the Target!','warn');
    uiwait(h,5);
end