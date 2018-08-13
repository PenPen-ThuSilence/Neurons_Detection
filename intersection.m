function value = intersection(k1,b1,k2,b2)

if k1 == k2
    % if k1 = k2, then we can't get specific intersection anyway.
    x = nan;
    y = nan;
elseif k1 == Inf
    x = b1;
    y = k2 * x + b2;
elseif k2 == Inf
    x = b2;
    y = k1 * x + b1;
else
    x = (b2 - b1) / (k1 - k2);
    y = k1 * x + b1;
end
    
value = [x, y];