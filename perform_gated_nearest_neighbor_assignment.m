function b = perform_gated_nearest_neighbor_assignment(measurements, X, D, thresh)

mahalanobis_dist = NaN(size(measurements, 2), size(X,2)); %costmatrix
b = NaN(1, size(measurements,2));
object_assign = NaN(size(b));
invCov_Z = NaN(2,2,sum(~isnan(X(1, :))));

for n = find(~isnan(X(1, :)))
    Drot = rotate_extent(X(:, n), D(:,:,n));
    invCov_Z(:,:,n) = Drot^-1;
end

%build cost matrix
for m = 1:size(measurements,2)
    if ~isnan(measurements(1,m))
        for n = find(~isnan(X(1, :)))
            nu = measurements(:, m) - X(1:2, n);
            mahalanobis_dist(m,n) = (transpose(nu)*invCov_Z(:,:,n)*nu)^2;
            if mahalanobis_dist(m,n) > thresh
                mahalanobis_dist(m,n) = Inf;
            end
        end
    end
end
mahalanobis_dist(isnan(mahalanobis_dist)) = Inf;

if size(X,2) == 1
    b(:) = 0;
    [~, idx] = min(mahalanobis_dist);
    b(idx) = 1;
else
    [assignment, cost] = munkres(mahalanobis_dist.');
    assignment = assignment.';
    for m = 1:size(assignment, 1)
        assign = find(assignment(m, :));
        if any(assign)
            b(m) = assign;
        elseif isnan(measurements(1,m))
            b(m) = NaN;
        else
            b(m) = 0;
        end
    end
end

%
%         if isempty(mahalanobis_dist(m,:))
%             object_assign(m) = 0;
%         elseif min(mahalanobis_dist(m,:)) == Inf
%             object_assign(m) = 0;
%         else
%             [~, object_assign(m)] = min(mahalanobis_dist(m,:));
%         end
%     else
%         object_assign(m) = NaN;
%     end
%     b = object_assign;
% end

%%old version
% function b = perform_gated_nearest_neighbor_assignment(measurements, X, D, thresh)
%
% mahalanobis_dist = NaN(size(measurements, 2), sum(~isnan(X(1, :)))); %costmatrix
% b = NaN(1, size(measurements,2));
% object_assign = NaN(size(b));
% invCov_Z = NaN(2,2,sum(~isnan(X(1, :))));
%
% for n = find(~isnan(X(1, :)))
%     Drot = rotate_extent(X(:, n), D(:,:,n));
%     invCov_Z(:,:,n) = Drot^-1;
% end
%
% %build cost matrix
% for m = 1:size(measurements,2)
%     if ~isnan(measurements(1,m))
%         for n = find(~isnan(X(1, :)))
%             if ~isnan(X(1, n))
%                 nu = measurements(:, m) - X(1:2, n);
%                 mahalanobis_dist(m,n) = (transpose(nu)*invCov_Z(:,:,n)*nu)^2;
%                 if mahalanobis_dist(m,n) > thresh || isnan(mahalanobis_dist(m,n))
%                     mahalanobis_dist(m,n) = Inf;
%                 end
%             else
%                 mahalanobis_dist(m,n) = Inf;
%             end
%
%         end
%         if isempty(mahalanobis_dist(m,:))
%             object_assign(m) = 0;
%         elseif min(mahalanobis_dist(m,:)) == Inf
%             object_assign(m) = 0;
%         else
%             [~, object_assign(m)] = min(mahalanobis_dist(m,:));
%         end
%     else
%         object_assign(m) = NaN;
%     end
%     b = object_assign;
% end