function [Xcs] = track_objects_PGAS(Z, Xcs, Ts, Te, Ccs, Dstarcs, F, Sigmavstarcs, Sigmax1, R, maxNactive)

for n = 1:maxNactive
    Xcs(:, :, n) = csmcsampleauxpfas(Z(:, :, :, n), Xcs(:, :, n), Ts(n), Te(n), Dstarcs(:, :, Ccs(n)), F, Sigmavstarcs(:, :, Ccs(n)), Sigmax1, R);
end

%state_est.Xcs = Xcs;
