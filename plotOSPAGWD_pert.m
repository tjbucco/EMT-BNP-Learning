GOSPA_PerfEMT = mean(OSPA_Perf_GWD_nolabels_pert, 2);
GOSPA_BEMT = mean(OSPA_NC_GWD_nolabels_pert, 2);
GOSPA_DPEMT = mean(OSPA_DP_GWD_nolabels_pert, 2);

timeaxis = (1:maxT).';

clear plot

figure

plot(timeaxis.', GOSPA_DPEMT)
hold on
plot(timeaxis.', GOSPA_BEMT)
hold on
plot(timeaxis.', GOSPA_PerfEMT)
hold off
grid

ax = gca;
legend("DPEMT", "BEMT", "PerfEMT")