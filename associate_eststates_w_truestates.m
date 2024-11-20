function assignment = associate_eststates_w_truestates(X, Xest)

error_location = 0;
if size(X,2) == 1
    I = 1;
else
[~, I] = find(~isnan(squeeze(X(1,:))), 1, 'last');
end

if size(Xest,2) == 1
    J = 1;
else
[~, J] = find(~isnan(squeeze(Xest(1,:))), 1, 'last');
end

cost_matrix = 10000 * ones(I, J);
for i = 1:I
    for j = 1:J
        if any(~isnan(X(1,i)) & ~isnan(Xest(1,j)))
            if ~isnan(X(1,i)) && ~isnan(Xest(1,j))
                num = norm(X(1:2, i) - Xest(1:2, j));
                exp_denom = 1;
            end
            cost_matrix(i,j) = num / exp(exp_denom);
        end
    end
end

assignment = munkres(cost_matrix.');
[omega, ~, ~] = find(assignment);

% da_remapped = NaN(size(da_est));
% for t = 1:T
%     [~, no_meas] = find(~isnan(squeeze(da_est(:,t).')), 1, 'last');
%     for m = 1:no_meas
%         if da_est(m,t) ~= 0
%             if ~any(find(assignment(da_est(m, t), :)))
%                 da_est(m, t) = I + 1;
%             else
%                 da_remapped(m, t) = find(assignment(da_est(m, t), :));
%             end
%         else
%             da_remapped(m, t) = da_est(m, t);
%         end
%     end
% end

end