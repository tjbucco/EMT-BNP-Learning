function [object_class_est] = match_objects_to_existingclasses(object_class_est, object_class_mcmc, object_indices_mcmc, DA_indices, t)

DA_to_mcmc = DA_indices(end - size(object_indices_mcmc,2) + 1:end);
object_indices_seq = unique(DA_to_mcmc);

if any(object_indices_seq > 0)
    object_indices_seq = object_indices_seq(object_indices_seq > 0);
    seq_to_mcmcindex = NaN(size(object_indices_seq));
    kprime = 0;
    for k = object_indices_seq
        kprime = kprime + 1;
        seq_to_mcmcindex(kprime) = object_indices_mcmc(find(DA_to_mcmc == k, 1, "first"));
        class_ind_mcmc = object_class_mcmc.Cest(seq_to_mcmcindex(kprime));
        if any(object_class_est.lambdaMest == object_class_mcmc.lambdaMest(class_ind_mcmc))
            class_ind_seq = find(object_class_est.lambdaMest == object_class_mcmc.lambdaMest(class_ind_mcmc), 1, "first");
            object_class_est.Cest(k) = class_ind_seq;
        else
            %if t == 1
            if k == 1
                object_class_est.Cest(k) = 1;
            else
                object_class_est.Cest(k) = max(object_class_est.Cest(1:k - 1)) + 1;
            end
            %else
            %    object_class_est.Cest(k) = max(object_class_est.Cest) + 1;
            %end
            object_class_est.Dest(:,:,object_class_est.Cest(k)) = object_class_mcmc.Dest(:,:,class_ind_mcmc);
            object_class_est.Qest(:,:,object_class_est.Cest(k)) = object_class_mcmc.Qest(:,:,class_ind_mcmc);
            object_class_est.lambdaMest(object_class_est.Cest(k)) = object_class_mcmc.lambdaMest(class_ind_mcmc);
        end
    end
end
object_class_est.Dest = object_class_est.Dest(:,:,1:max(object_class_est.Cest));
object_class_est.Qest = object_class_est.Qest(:,:,1:max(object_class_est.Cest));
object_class_est.lambdaMest = object_class_est.lambdaMest(1:max(object_class_est.Cest));
if isrow(object_class_est.lambdaMest)
    object_class_est.lambdaMest = object_class_est.lambdaMest';
end

end