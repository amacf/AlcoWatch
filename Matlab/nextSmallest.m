% Returns the first element of items that is larger than x
% Otherwise return -1
function n = nextSmallest(items, x)
    i = 1;
    n = -1;
    if isempty(items)
        return;
    end
    while n < x
        n = items(i);
        i=i+1;
        if i==length(items)+1
            n = -1;
            return;
        end
    end
end