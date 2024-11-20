GOSPA_PerfEMT = mean(OSPA_Perf_GWD_nolabels_pert, 2);
% GOSPA_PerfEMT = mean(OSPA_Perf_GWD_nolabels_pert_new, 2);

GOSPA_FEMT = mean(OSPA_FUSION_GWD_nolabels_pert, 2);
% GOSPA_FEMT = mean(OSPA_FUSION_GWD_nolabels_pert_new, 2);
%GOSPA_FEMT = mean(OSPA_FUSION_GWD_pert, 2);

GOSPA_BEMT = mean(OSPA_NC_GWD_nolabels_pert, 2);
% GOSPA_BEMT = mean(OSPA_NC_GWD_nolabels_pert_new, 2);
% GOSPA_BEMT = mean(OSPA_NC_GWD_pert, 2);

GOSPA_DPEMT = mean(OSPA_DP_GWD_nolabels_pert, 2);
% GOSPA_DPEMT = mean(OSPA_DP_GWD_nolabels_pert_new, 2);
% GOSPA_DPEMT = mean(OSPA_DP_GWD_pert, 2);

timeaxis = (1:maxT).';

GOSPA_results = [["time"; timeaxis], ["PRO"; GOSPA_DPEMT], ["BL"; GOSPA_BEMT], ["SEP"; real(GOSPA_FEMT)], ["KNO"; real(GOSPA_PerfEMT)]];
writematrix(GOSPA_results, "OSPAGW_results_journal.dat", 'Delimiter', 'tab');


clear plot

figure

plot(timeaxis.', GOSPA_DPEMT)
hold on
plot(timeaxis.', GOSPA_BEMT)
hold on
plot(timeaxis.', GOSPA_PerfEMT)
hold on
plot(timeaxis.', GOSPA_FEMT)
hold off
grid

ax = gca;
%ax.LineStyleOrder = ["-"; "--"];
% legend("DPEMT", "BEMT", "PerfEMT")
legend("PRO", "BL", "KNO", "SEP")