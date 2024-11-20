function assignment = associate_esttracks_w_truetracks_notmunkres(X, Xest)

error_location = 0;
if size(X,3) == 1
    I = 1;
else
[~, I] = find(~isnan(squeeze(X(1,:,:))), 1, 'last');
end

if size(Xest,3) == 1
    J = 1;
else
[~, J] = find(~isnan(squeeze(Xest(1,:,:))), 1, 'last');
end

T = size(X,2);

cost_matrix = 10000 * ones(I, J);
for i = 1:I
    for j = 1:J
        if any(~isnan(X(1,:,i)) & ~isnan(Xest(1,:,j)))
            num = 0;
            exp_denom = 0;
            for t = 1:T
                if ~isnan(X(1,t,i)) && ~isnan(Xest(1,t,j))
                    num = norm(X(1:2, t, i) - Xest(1:2, t, j)) + num;
                    exp_denom = exp_denom + 1;
                end
            end
            cost_matrix(i,j) = num / exp(exp_denom);
        end
    end
end

[~, assignment] = min(cost_matrix);

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