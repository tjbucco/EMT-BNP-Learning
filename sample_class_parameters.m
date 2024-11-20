function [Dstarcs, Qstarcs, lambdaMstarcs] = sample_class_parameters(object_class_est, Ccs, mcmc, hyper, j, c)

Mcsum = sum(object_class_est.Msum(Ccs == c));
Zcsqsum = sum(object_class_est.Zsqsum(:, Ccs == c), 2);
nocsum = sum(object_class_est.existsteps(Ccs == c)) + 1;
nocsum = sum(object_class_est.existsteps(Ccs == c));
ncsum = sum(Ccs == c);
Vcsqsum = sum(object_class_est.Vsqsum(:, Ccs == c), 2);
Vcsqsum_pos = sum(Vcsqsum(1:2));
Vcsqsum_velo = sum(Vcsqsum(3:4));
if j >= mcmc.Dfs
    alstarc = hyper.al0 + 1/2*Mcsum*[1; 1];
    blstarc = hyper.bl0 + 1/2*Zcsqsum;

    lsqstarcinv = gamrnd(alstarc, 1./blstarc, 2, 1);
    lsqstarc = 1./lsqstarcinv;
    Dstarcs = diag(lsqstarc);
else
    Dstarcs = hyper.D0;
end

if j >= mcmc.Qfs
    avstarc = hyper.av0 + 1/2*(nocsum - ncsum)*[2; 2];
    avstarc = hyper.av0 + 1/2*(nocsum - ncsum)*[2; 2];
    bvstarc = hyper.bv0 + 1/2*[Vcsqsum_pos; Vcsqsum_velo];

    ssqstarcinv = gamrnd(avstarc, 1./bvstarc, 2, 1);
    ssqstarc = 1./ssqstarcinv;
    Qstarcs = diag(ssqstarc([1, 1, 2, 2]));
else
    Qstarcs = hyper.Sigmav0;
end

if j >= mcmc.Dfs
    aMstarc = hyper.aM0 + Mcsum;
    bMstarc = hyper.bM0 + (nocsum - ncsum) + 1;

    lambdaMstarcs = gamrnd(aMstarc, 1./bMstarc);
else
    lambdaMstarcs = hyper.ElambdaM;
end

end