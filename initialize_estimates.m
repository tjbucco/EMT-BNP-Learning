function [state_est, Ccs, Dstarcs, Qstarcs, lambdaMstarcs] = initialize_estimates(init_tracking, data, hyper, script, mcmc)

total_no_of_objects = init_tracking.maxNactive;
maxTactive = init_tracking.maxTactive;

Xcs = NaN(4, script.maxT, total_no_of_objects);
for n = 1:total_no_of_objects
    if nansum(init_tracking.M(:, n))
        Xcs(1:2, :, n) = squeeze(nanmean(init_tracking.Z(:, :, :, n), 2));

        ti = find(init_tracking.M(:, n) == 0);

        ti = ti(ti > find(init_tracking.M(:, n) > 0, 1, 'first'));
        ti = ti(ti < find(init_tracking.M(:, n) > 0, 1, 'last')).';

        for t = ti
            for tp = t:init_tracking.Te(n)
                if ~any(ti == tp)
                    Xcs(1:2, t, n) = Xcs(1:2, t - 1, n) + (Xcs(1:2, tp, n) - Xcs(1:2, t - 1, n))/(tp - t + 1);

                    break
                end
            end
        end


        if (find(init_tracking.M(:, n) > 0, 1, 'last') - find(init_tracking.M(:, n) > 0, 1, 'first')) > 1
            ni = find(init_tracking.M(:, n) == 0);
            ni = ni(ni < find(init_tracking.M(:, n) > 0, 1, 'first'));
            ni = ni(ni >= init_tracking.Ts(n)).';

            for t = fliplr(ni)
                Xcs(1:2, t, n) = Xcs(1:2, t + 1, n) - (Xcs(1:2, t + 2, n) - Xcs(1:2, t + 1, n));
            end

            ni = find(init_tracking.M(:, n) == 0);
            ni = ni(ni > find(init_tracking.M(:, n) > 0, 1, 'last'));
            ni = ni(ni <= init_tracking.Te(n)).';

            for t = ni
                Xcs(1:2, t, n) = Xcs(1:2, t - 1, n) + (Xcs(1:2, t - 1, n) - Xcs(1:2, t - 2, n));
            end
        else
            ni = find(init_tracking.M(:, n) == 0);
            ni = ni(ni >= init_tracking.Ts(n));
            ni = ni(ni <= init_tracking.Te(n)).';

            for t = ni
                Xcs(1:2, t, n) = Xcs(1:2, find(init_tracking.M(:, n) > 0, 1, 'first'), n);
            end
        end

        for t = init_tracking.Ts(n):(init_tracking.Te(n) - 1)
            Xcs(3:4, t, n) = (Xcs(1:2, t + 1, n) - Xcs(1:2, t, n))./[hyper.F(1, 3); hyper.F(2, 4)];
        end

        if (init_tracking.Te(n) - init_tracking.Ts(n) + 1) > 1
            Xcs(3:4, init_tracking.Te(n), n) = Xcs(3:4, init_tracking.Te(n) - 1, n);
        else
            Xcs(3:4, init_tracking.Te(n), n) = zeros(2, 1);
        end
    else
        Xcs(:, init_tracking.Ts(n):init_tracking.Te(n), n) = zeros(4, init_tracking.Te(n)-init_tracking.Ts(n)+1);
    end

    Xcs(:,:, n) = csmcsampleauxpfas(init_tracking.Z(:, :, :, n), Xcs(:, :, n), init_tracking.Ts(n), init_tracking.Te(n), hyper.D0, hyper.F, hyper.Sigmav0, hyper.Sigmax1, mcmc.Rinit);
end

Ysq = NaN(2, maxTactive, total_no_of_objects);

%C = NaN(1, total_no_of_objects, (mcmc.L - mcmc.bis)/mcmc.DPsi);

Ccs = NaN(1, total_no_of_objects);
Ccs(1, :) = ones(1, total_no_of_objects);

% Dstar = NaN(2, 2, total_no_of_objects, (mcmc.L - mcmc.bis)/mcmc.DPsi);
% D = NaN(2, 2, total_no_of_objects, (mcmc.L - mcmc.bis)/mcmc.DPsi);

Dstarcs = NaN(2, 2, total_no_of_objects);
Dstarcs(:, :, 1) = hyper.D0;

% Sigmavstar = NaN(4, 4, total_no_of_objects, (mcmc.L - mcmc.bis)/mcmc.DPsi);
% Sigmav = NaN(4, 4, total_no_of_objects, (mcmc.L - mcmc.bis)/mcmc.DPsi);

Qstarcs = NaN(4, 4, total_no_of_objects);
Qstarcs(:, :, 1) = hyper.Sigmav0;

% lambdaMstar = NaN(1, total_no_of_objects, (mcmc.L - mcmc.bis)/mcmc.DPsi);
% lambdaM = NaN(1, total_no_of_objects, (mcmc.L - mcmc.bis)/mcmc.DPsi);

lambdaMstarcs = NaN(1, total_no_of_objects);
lambdaMstarcs(1, 1) = hyper.ElambdaM;

state_est = struct("Xcs", Xcs, "Ysq", Ysq);
%[Xcs, Ysq, C, Ccs, Dstar, D, Sigmavstar, Sigmav, Sigmavstarcs, Sigmavstarcs, lambdaMstar, lambdaM, lambdaMstarcs]
end