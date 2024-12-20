function [Xnew, b, a] = get_nearest_state_munkres(X, Xest)


dist = NaN(size(Xest, 2), size(X,2)); %costmatrix
b = NaN(1, size(Xest,2));
object_assign = NaN(size(b));

%build cost matrix
for m = 1:size(Xest,2)
    if ~isnan(Xest(1,m))
        for n = find(~isnan(X(1, :)))
            nu = Xest(:, m) - X(1:2, n);
            dist(m,n) = transpose(nu)*nu;
        end
    end
end
dist(isnan(dist)) = Inf;

if size(X,2) == 1
    b(:) = 0;
    [~, idx] = min(dist);
    b(idx) = 1;
else
    [assignment, cost] = munkres(dist.');
    assignment = assignment.';
    for m = 1:size(assignment, 1)
        assign = find(assignment(m, :));
        if any(assign)
            b(m) = assign;
        elseif isnan(Xest(1,m))
            b(m) = NaN;
        else
            b(m) = 0;
        end
    end
end

Xnew = Xest;
a = find(b>0);
Xnew(:,a) = X(:,b(b>0)); 

end