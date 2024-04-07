clear; close all;
path(path, './lib');
path(path, './mat_files');
path(path, './ftp_plot_scripts');
run('./lib/cvx-w64/cvx/cvx_setup.m')


run("./ftp_plot_scripts/fig1_beamtraining_overhead.m")
run("./ftp_plot_scripts/design_b_superres.m")
run("./ftp_plot_scripts/design_c_peaks.m")
run("./ftp_plot_scripts/design_d_beam.m")
run("./ftp_plot_scripts/micro_benchmark_angle.m")
run("./ftp_plot_scripts/micro_benchmark_gain.m")
run("./ftp_plot_scripts/snr_perf.m")
run("./ftp_plot_scripts/beamtraining_overhead.m")
run("./ftp_plot_scripts/super_resolution.m")
run("./ftp_plot_scripts/ofdm_beam.m")
run("./ftp_plot_scripts/ofdm_perf.m")