function [object_class_est, C, Dstar, Qstar, lambdaMstar, Ccs, Dstarcs, Qstarcs, lambdaMstarcs] = cluster_objects(init_tracking, Xcs, object_class_est_prev, hyper, factorials, mcmc, ell, C, Dstar, Qstar, lambdaMstar)

object_class_est = object_class_est_prev;
existsteps = zeros(1, init_tracking.maxNactive);
Msum = zeros(1, init_tracking.maxNactive);
denom_factorial = zeros(1, init_tracking.maxNactive);
Zsqsum = zeros(2, init_tracking.maxNactive);
Vsqsum = zeros(4, init_tracking.maxNactive);

for i = 1:init_tracking.maxNactive
    denom_factorial_exp = 0;
    for t = init_tracking.Ts(i):init_tracking.Te(i)
        rot = calcrotmat(Xcs(:,t, i));
        premeasurements_i = rot.'* (init_tracking.Z(:, :, t, i) - Xcs(1:2, t, i));
        existsteps(i) = existsteps(i) + 1;
        Msum(i) = Msum(i) + init_tracking.M(t,i);
        denom_factorial_exp = denom_factorial_exp + log(factorials(max(min(init_tracking.M(t,i), 150),1)));
        Zsqsum(:, i) = Zsqsum(:, i) + nansum(premeasurements_i.^2, 2);
        if t == init_tracking.Ts(i)
            driving_noise_est = [0;0;0;0];
        else
            driving_noise_est = Xcs(:,t,i) - hyper.F*Xcs(:,t-1,i);
        end
        Vsqsum(:,i) = Vsqsum(:,i) + nansum(driving_noise_est.^2, 2);
    end
    denom_factorial(i) = denom_factorial_exp;
end
if init_tracking.maxNactive > 0
    object_class_est.existsteps = existsteps; object_class_est.Msum = Msum; object_class_est.Zsqsum = Zsqsum; object_class_est.Vsqsum = Vsqsum;
else
    object_class_est = object_class_est_prev;
    return
end
%Gibbs sampler algorithm - see bernds' thesis
%%%%%%%%%
%initialization:
% for ell = 1:mcmc.L
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

%% normal iteration:
if ell == 1
%% first sample only, for all objects i:
%set each object to its own class
C(:, 1) = 1:init_tracking.maxNactive;

%sample D*_c and Q*_c and lambda*_c from base pdf
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

lambdaMstar(:, 1) = gamrnd(hyper.aM0(:, ones(1, max(C(:,1)))), 1./hyper.bM0(:, ones(1, max(C(:,1)))), 1, max(C(:,1)));

Ccs = C(:,1); Dstarcs = Dstar(:,:,:,1); 
Qstarcs = Qstar(:,:,:,1); lambdaMstarcs = lambdaMstar(:,1);

%% samples 2 ... K, for all objects i:
else
    Qtent = [];
    Dtent = [];
    lambdatent = [];
    for i = 1:init_tracking.maxNactive

        %sample c_i^{n, k} = c w/ prob eq(5.108) and eq(5.110)
        Cwoi = [C(1:(i - 1), ell).' C((i + 1):end, ell - 1).'];
        if isempty(Cwoi)
            Cwoi = 0;
        end
        Qstarellminus1 = Qstar(:, :, :, ell - 1);
        Qstarellminus1 = cat(3, Qstarellminus1, Qtent);
        Dstarellminus1 = Dstar(:, :, :, ell - 1);
        Dstarellminus1 = cat(3, Dstarellminus1, Dtent);
        lambdastarellminus1 = lambdaMstar(:, ell - 1);
        lambdastarellminus1 = cat(1, lambdastarellminus1, lambdatent);
        a = calc_indicator_probs(Qstarellminus1, Dstarellminus1, lambdastarellminus1, hyper.Sigmax1, Zsqsum(:, i), init_tracking.maxNactive, C(:, ell - 1), Cwoi, hyper, existsteps(i), Msum(i), Vsqsum(:, i), denom_factorial(i), Xcs(:,init_tracking.Ts(i), i));
        C(i, ell) = find(mnrnd(1, a)) - 1;

        %if c_i^{n, k} is new class, sample D*tent and Q*tent from eq(5.113)
        if C(i, ell) == 0
            C(i, ell) = max(C(:, ell)) + 1;
            [Dtent, Qtent, lambdatent] = sample_class_parameters(object_class_est, C(:, ell), mcmc, hyper, ell, C(i, ell));
        end
    end

    %cull unused classes
    [~, ~, C(:, ell)] = unique(C(:, ell));

    %once c_i^{n, k} is sampled for all objects, sample D*_c and Q*_c from eq(5.113)
    for c = 1:max(C(:, ell))
        [Dstar(:, :, c, ell), Qstar(:, :, c, ell), lambdaMstar(c, ell)] = sample_class_parameters(object_class_est, C(:, ell), mcmc, hyper, ell, c);
    end



    %     if ~mod(ell, mcmc.DPsi) && ell > mcmc.bis
    %         C(1, :, (ell - mcmc.bis)/mcmc.DPsi) = Ccs;
    %         Dstar(:, :, :, (ell - mcmc.bis)/mcmc.DPsi) = Dstarcs;
    %         Sigmavstar(:, :, :, (ell - mcmc.bis)/mcmc.DPsi) = Sigmavstarcs;
    %         lambdaMstar(1, :, (ell - mcmc.bis)/mcmc.DPsi) = lambdaMstarcs;
    %     end

    Ccs = C(:, ell);
    Dstarcs = Dstar(:, :, :, ell);
    Qstarcs = Qstar(:, :, :, ell);
    lambdaMstarcs = lambdaMstar(:, ell);

end

end