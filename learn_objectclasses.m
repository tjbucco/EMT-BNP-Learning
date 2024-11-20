function [object_class_est] = learn_objectclasses(state_est, daoed_est, object_class_est_prev, measurements, data, hyper, script, mcmc, t)

%instantiate variables
C = zeros(daoed_est.total_no_of_objects, mcmc.L);
Dstar = NaN(2, 2, daoed_est.total_no_of_objects, mcmc.L);
Qstar = NaN(4, 4, daoed_est.total_no_of_objects, mcmc.L);

%update sufficient statistics
% if t > 1
%     existsteps = [object_class_est_prev.existsteps zeros(1, daoed_est.no_of_newobjectst)];
%     Msum = [object_class_est_prev.Msum zeros(1, daoed_est.no_of_newobjectst)];
%     Zsqsum = cat(2, object_class_est_prev.Zsqsum, zeros(2, daoed_est.no_of_newobjectst));
%     Vsqsum = cat(2,object_class_est_prev.Vsqsum, zeros(4, daoed_est.no_of_newobjectst));
% end
exist_notclutter_ind = (daoed_est.exist_ind_all > 0) & (daoed_est.notclutter_ind > 0);
existsteps = object_class_est_prev.existsteps;
Msum = object_class_est_prev.Msum;
Zsqsum = object_class_est_prev.Zsqsum;
Vsqsum = object_class_est_prev.Vsqsum;

for i = find(daoed_est.exist_ind_all > 0)
    rot = calcrotmat(state_est.Xest(:,i));
    premeasurements_i = rot.'* (measurements(:, daoed_est.DA == i) - state_est.Xest(1:2,i));
    existsteps(i) = object_class_est_prev.existsteps(i) + 1;
    Msum(i) = object_class_est_prev.Msum(i) + sum(daoed_est.DA == i);
    Zsqsum(:, i) = object_class_est_prev.Zsqsum(:, i) + sum(premeasurements_i.^2, 2);
    Vsqsum(:,i) = object_class_est_prev.Vsqsum(:,i) + sum(state_est.driving_noise_est(:, i).^2, 2);
end
if any(daoed_est.exist_ind_all > 0)
    object_class_est.existsteps = existsteps; object_class_est.Msum = Msum; object_class_est.Zsqsum = Zsqsum; object_class_est.Vsqsum = Vsqsum;
else
    object_class_est = object_class_est_prev;
    return
end
%Gibbs sampler algorithm - see bernds' thesis
%%%%%%%%%
%initialization:
% for ell = 1:mcmcL
%
%     %set each object to its own class
%     C(1, :, ell) = 1:daoed_est.total_no_of_objects;
%
%     %sample D*_c and Q*_c from base pdf
%     lsqinv = gamrnd(hyper.al0(:, ones(1, max(C(ell)))), 1./hyper.bl0(:, ones(1, max(C(ell)))), 2, max(C(ell)));
%     lsq = 1./lsqinv;
%     for c = 1:max(C(ell))
%         Dstar(:, :, c) = diag(lsq(:, c));
%     end
%
%     ssqinv = gamrnd(hyper.av0(:, ones(1, max(C(ell)))), 1./hyper.bv0(:, ones(1, max(C(ell)))), 2, max(C(ell)));
%     ssq = 1./ssqinv;
%     for c = 1:max(C(ell))
%         Qstar(:, :, c) = diag(ssq([1, 1, 2, 2], c));
%     end
% end

%normal iteration:
%first sample only, for all objects i:
%set each object to its own class
C(:, 1) = 1:daoed_est.total_no_of_objects;

%sample D*_c and Q*_c from base pdf
lsqinv = gamrnd(hyper.al0(:, ones(1, max(C(:,1)))), 1./hyper.bl0(:, ones(1, max(C(:,1)))), 2, max(C(:,1)));
lsq = 1./lsqinv;
for c = 1:max(C(:,1))
    Dstar(:, :, c, 1) = diag(lsq(:, c));
end

ssqinv = gamrnd(hyper.av0(:, ones(1, max(C(:,1)))), 1./hyper.bv0(:, ones(1, max(C(:,1)))), 2, max(C(:,1)));
ssq = 1./ssqinv;
for c = 1:max(C(:,1))
    Qstar(:, :, c, 1) = diag(ssq([1, 1, 2, 2], c));
end

%samples 2 ... K, for all objects i:
for ell = 2:mcmc.L
    Qtent = [];
    Dtent = [];
    for i = 1:daoed_est.total_no_of_objects

        %sample c_i^{n, k} = c w/ prob eq(5.108) and eq(5.110)
        Cwoi = [C(1:(i - 1), ell).' C((i + 1):end, ell - 1).'];
        if isempty(Cwoi)
            Cwoi = 0;
        end
        intial_ind = i >= daoed_est.total_no_of_objects - daoed_est.no_of_newobjectst;
        Qstarellminus1 = Qstar(:, :, :, ell - 1);
        Qstarellminus1 = cat(3, Qstarellminus1, Qtent);
        Dstarellminus1 = Dstar(:, :, :, ell - 1);
        Dstarellminus1 = cat(3, Dstarellminus1, Dtent);
        if intial_ind || t == 1
            a = calc_indicator_probs_seq(Qstarellminus1, Dstarellminus1, hyper.Sigmax1, Zsqsum(:, i), daoed_est.total_no_of_objects, C(:, ell - 1), Cwoi, hyper, existsteps(i), Msum(i), Vsqsum(:, i), intial_ind, state_est.Xest(:, i));
        else
            a = calc_indicator_probs_seq(Qstarellminus1, Dstarellminus1, hyper.Sigmax1, Zsqsum(:, i), daoed_est.total_no_of_objects, C(:, ell - 1), Cwoi, hyper, existsteps(i), Msum(i), Vsqsum(:, i), intial_ind);
        end
        C(i, ell) = find(mnrnd(1, a)) - 1;

        %if c_i^{n, k} is new class, sample D*tent and Q*tent from eq(5.113)
        if C(i, ell) == 0
            C(i, ell) = max(C(:, ell)) + 1;
            [Dtent, Qtent] = sample_class_parameters(object_class_est, C(:, ell), mcmc, hyper, ell, C(i, ell));
        end
    end

    %cull unused classes
    [~, ~, C(:, ell)] = unique(C(:, ell));

    %once c_i^{n, k} is sampled for all objects, sample D*_c and Q*_c from eq(5.113)
    for c = 1:max(C(:, ell))
        [Dstar(:, :, c, ell), Qstar(:, :, c, ell)] = sample_class_parameters(object_class_est, C(:, ell), mcmc, hyper, ell, c);
    end



    %     if ~mod(ell, mcmc.DPsi) && ell > mcmc.bis
    %         C(1, :, (ell - mcmc.bis)/mcmc.DPsi) = Ccs;
    %         Dstar(:, :, :, (ell - mcmc.bis)/mcmc.DPsi) = Dstarcs;
    %         Sigmavstar(:, :, :, (ell - mcmc.bis)/mcmc.DPsi) = Sigmavstarcs;
    %         lambdaMstar(1, :, (ell - mcmc.bis)/mcmc.DPsi) = lambdaMstarcs;
    %     end


end

[unq_samps,~,ic_unq] = unique(C(:, mcmc.bis + 1:end).', 'rows');
map_i = mode(ic_unq);
object_class_est.Cest = unq_samps(map_i,:).';
object_class_est.Dest = squeeze(nanmean(Dstar(:, :, :, mcmc.bis + 1:end), 4));
object_class_est.Qest = squeeze(nanmean(Qstar(:, :, :, mcmc.bis + 1:end), 4));

end