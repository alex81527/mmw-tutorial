run('./cvx-w64/cvx/cvx_setup.m')

load("./exp_data/swift2.0/ofdm-5gNR/data1.mat");
BPU = data.bpu;
PA =  data.pa;
PA.N_BEAM = 124;
create_cs_codebook(124, PA);
cb = "./codebooks/cs.mat";
[H30, r30, SNR30, n30, BPU, PA, maxk_pos30, maxk_pks30]=debug_peaks(data.rx_cs_probes, BPU, PA);
if BPU.DO_OFDM == 1
    [maxk_pos11, maxk_pks11] = super_resolution(H30,0);
else
    maxk_pos11 = maxk_pos30;
    maxk_pks11 = maxk_pks30;
end
[v_multipath,v2_multipath,v_dominantpath,path_ang,path_gain,unique_pos] = cs_multipath_algo(cb, BPU, PA, maxk_pos11, maxk_pks11, SNR30, data.rx_cs_probes, r30,1);

% PA.N_BEAM = length(data.az);
% [H29, r29, SNR29, n29, BPU, PA, maxk_pos29, maxk_pks29]=debug_peaks(data.rx_ex_probes, BPU, PA);
% figure; plot(data.az, r29);

v2 = v_multipath./v_multipath(1).*exp(1j*PA.PHASE_CAL);
v3 = v_dominantpath./v_dominantpath(1).*exp(1j*PA.PHASE_CAL);
v4 = v2_multipath./v2_multipath(1).*exp(1j*PA.PHASE_CAL);

PA.N_BEAM = 124;
[H31, r31, SNR31, n31, BPU, PA, maxk_pos31, maxk_pks31]=debug_peaks(data.rx_aco_probes, BPU, PA);
[rel_phase,rel_phase_avg,rel_phase_std] = get_rel_phase(r31);
aco_activated_ant = [1:1+124/4];
aco_sv = exp(1j*[0;-rel_phase_avg(aco_activated_ant(2:end)-1)]);
% aco_sv = sum(steervec(pa.getElementPosition()/PA.LAM, [-30 45;0 0]),2).*exp(1j*PA.PHASE_CAL);

sv = [v2 aco_sv];
legs = ["SWIFT", "ACO"];
% legs = [sprintf("SWIFT,%.1f",data.snr_cs_multipath(end)), ...
%     sprintf("CS,%.1f",data.snr_cs_dominantpath(end)), ...
%     sprintf("ACO,%.1f (ad,%.1f)",data.snr_aco(end),data.snr_11ad(end))];
% sv2beam(sv,[-90:1:90],[-50:1:50],PA,legs)

if BPU.DO_OFDM == 1
    PA.N_BEAM = length(data.aco_cbsize);
    [H40, r40, SNR40, n40, BPU, PA, maxk_pos40, maxk_pks40]=debug_peaks(data.rx_cs_multipath, BPU, PA);
    [maxk_pos41, maxk_pks41] = super_resolution(H40,1,[12]);
%     [maxk_pos41, maxk_pks41] = super_resolution(H40,0,logspace( -2, 1, 20 ));
end
pa = get_phased_array(PA.FREQ);
%% super res
close all;
fontsize = 15;
fig_size = [1 1 6 3];
linewidth = 2.5;
markersz = 8;

H = H40;
fftsize = size(H,1);
b = ifft(H);
hh = [b(end-fftsize/2+1:end); b(1:fftsize/2)];
fs = 200e6;
res = 1/fs/10;
Analog_BW = 160e6;
t = [0:fftsize-1]./fs;
[M,I] = max(abs(hh));
intial_tau = I/fs;
        
t_search = intial_tau + [-20e-9:res:20e-9];
[T2,Tau2] = ndgrid(t,t_search);
h_sinc_mat2 = sinc(Analog_BW*(T2-Tau2)); % get sincs at tau TOF

A = h_sinc_mat2;
b = hh./max(abs(hh)); %h_sinc_noise(1:length(h_sinc));
n = size(A,2);
gamma = logspace( -1, 0, 20 );
x_final = zeros(n,length(gamma));

mse_error = zeros(2,1);
k = 12;
cvx_begin quiet
    variable x(n) complex
    minimize( norm(A*x-b)+gamma(k)*norm(x,1) )
cvx_end

h_est = x;
[~,loc,~,prominence] = findpeaks(abs(h_est)./max(abs(h_est)),"MinPeakProminence", 0.1, "MinPeakHeight", 0.1);
[M,I] = maxk(prominence,2);

loc = loc(I);
tau = t_search(loc);
[T,Tau] = ndgrid(t,tau);
h_sinc_mat = sinc(Analog_BW*(T-Tau));
h_est2 = lsqminnorm(h_sinc_mat, hh, 1e-10,"warn");
h_sinc_mat_wt3 = h_sinc_mat*diag(h_est2);
h_sinc2 = sum(h_sinc_mat_wt3,2);

% for plotting
t2 = [0:0.1:fftsize-1]./fs;
[T2,Tau2] = ndgrid(t2,tau);
h_sinc_mat2 = sinc(Analog_BW*(T2-Tau2)); % get sincs at tau TOF
h_sinc_mat_wt5 = h_sinc_mat2*diag(h_est2);
h_sinc2_plot = sum(h_sinc_mat_wt5,2);

fig = figure('Units','inches', 'Position', fig_size);
shift = fftsize/2/fs;
scale = max(abs(hh));
stem((tau(1)-shift)*1e9,abs(h_est2(1))./scale, 'b-.', 'linewidth',linewidth,'MarkerSize',markersz); hold on;
plot((t2-shift)*1e9,abs(h_sinc_mat_wt5(:,1))./scale, 'b-.', 'linewidth',linewidth,'MarkerSize',markersz); hold on;
if length(tau)>1
    stem((tau(2)-shift)*1e9,abs(h_est2(2))./scale, 'r:', 'linewidth',linewidth,'MarkerSize',markersz); hold on;
    plot((t2-shift)*1e9,abs(h_sinc_mat_wt5(:,2))./scale, 'r:', 'linewidth',linewidth,'MarkerSize',markersz); hold on;
end
%                 plot(t2*1e9,abs(h_sinc2_plot), 'mo-', 'linewidth',1); hold on;
plot((t-shift)*1e9,abs(h_sinc2)./scale, 'gs-', 'linewidth',linewidth); hold on;
plot((t-shift)*1e9,abs(hh)./scale, 'ko-', 'linewidth',linewidth);
hold off; grid on; grid minor;
%                 title(sprintf("MSE=%.3f",mse_error(jj)));
% ylim([-40,10]);
xlim(([t_search(1)-5e-9 t_search(end)+2e-9] - shift)*1e9);
xlabel('ToF (ns)'); ylabel('Normalized magnitude')
set(gca,'fontsize',fontsize)
if length(tau)>1
    legend('','Est. path 1','', 'Est. path 2', 'Fitted CIR', 'Measured CIR')
    %                     legend('Est. path 1','Est. path 2', 'Fitted CIR', 'Measured CIR')
else
    legend('','Est. path 1','Fitted CIR','Mixed CIR')
end
exportgraphics(fig,"./figures/ofdm_superres.pdf",'Resolution',300);

%%%%
% a = conj(steervec(pa.getElementPosition()/PA.LAM, path_ang)).'*(v2./exp(1j*PA.PHASE_CAL));
% fprintf("TX beam complex gain's phase:"); disp(angle(a.'));
% fprintf("Estimated path complex gain's phase:"); disp(angle(path_gain.'));
% fprintf("Measured peaks' phase:"); disp(angle(h_est2));
% fprintf("Difference: "); disp(angle([path_gain.*a./h_est2]).')
%% beam plot
% close all;
fontsize = 18;
fig_size = [1 1 6 3];
linewidth = 2.5;
markersz = 8;
sysname = "CAMEO";

az=[-90:1:90];
el=[-50:1:50];

sv2 = exp(-1j*2*pi/4.*sv2psh(v2./exp(1j*PA.PHASE_CAL)));
[PAT_1,AZ_ANG,EL_ANG] = pattern(pa,PA.FREQ,az,el,...
    'PropagationSpeed',physconst('LightSpeed'),'Type','power','Normalize',true,...
    'CoordinateSystem','polar','Weights',sv2);
sv3 = exp(-1j*2*pi/4.*sv2psh(aco_sv./exp(1j*PA.PHASE_CAL)));
[PAT_2,AZ_ANG,EL_ANG] = pattern(pa,PA.FREQ,az,el,...
    'PropagationSpeed',physconst('LightSpeed'),'Type','power','Normalize',true,...
    'CoordinateSystem','polar','Weights',sv3);



% find angles from ACO
PAT_ACO = PAT_2;
PAT_CAMEO = PAT_1;
% path ang 1
[C,I] = max(PAT_ACO(:));
[I1,I2] = ind2sub(size(PAT_ACO),I);
fprintf("max az=%d, el=%d\n",az(I2), el(I1));
% path ang 2
PAT_ACO = squeeze(PAT_2(:,1:find(az==-20)));
[C,I] = max(PAT_ACO(:));
[I3,I4] = ind2sub(size(PAT_ACO),I);
fprintf("max az=%d, el=%d\n",az(I4), el(I3));



fig = figure('Units','inches', 'Position', fig_size);
color = 'brg';
PAT_azcut = squeeze(PAT_CAMEO(find(EL_ANG==0),:));
PAT_azcut_db = db(PAT_azcut, 'power');
polarplot(deg2rad(az), PAT_azcut_db, ['b-'], 'LineWidth', linewidth, "MarkerSize", markersz);hold on;
polarplot(deg2rad(az(I2))*ones(1,2), [-35 0], ['k--'], 'LineWidth', linewidth, "MarkerSize", markersz);hold on;
polarplot(deg2rad(az(I4))*ones(1,2), [-35 0], ['k--'], 'LineWidth', linewidth, "MarkerSize", markersz);hold on;
ax = gca;
% ax.RTickLabel = {""}; % remove ticklabels
% subtitle("Normalized gain (dB)", "Position",[0,-47]);
ax.ThetaDir = 'clockwise';
set(gca,'ThetaZeroLocation','top','FontSize',fontsize)
set(gca,'fontsize',fontsize)
% thetaticks(-90:30:90);
thetalim([-90 90]);
rlim([-35 0]);
legend([sysname, "Ground truth"], 'Location', 'northoutside', 'NumColumns',2, 'Fontsize',fontsize);
exportgraphics(fig,"./figures/ofdm_beam.pdf",'Resolution',300);


%% SNR barplot
% path(path, "./legendflex-pkg-master/setgetpos_V1.2");
% path(path, "./legendflex-pkg-master/matlab-hatchfill2-master");
fontsize = 22;
fig_size = [1 1 6 4];
linewidth = 2.5;
markersz = 8;
color = ['c','m','m','m','g','g']; % CAMEO, ACO, UbiG
% close all;
fig = figure('Units','inches', 'Position', fig_size);
hold on;
y3 = [data.snr_cs_multipath]; 
y4 = [data.snr_11ad];
bars(1) = bar(1,[y3]);
bars(2) = bar(2,[y4]);

% set(gca,'YScale','log');
% ylim([5e-7 1e2]);
% change facecolor
for ii=1:2
    bars(ii).FaceColor = color(ii);
end
    
% apply hatch pattern to bars
% hatchfill2(bars(1), 'single', 'HatchAngle', 45, 'HatchDensity', 40, 'HatchColor', 'black');
% hatchfill2(bars(2), 'single', 'HatchAngle', 45, 'HatchDensity', 20, 'HatchColor', 'black');
% hatchfill2(bars(2), 'cross', 'HatchAngle', 45, 'HatchDensity', 40, 'HatchColor', 'black');
% hatchfill2(bars(4), 'cross', 'HatchAngle', 45, 'HatchDensity', 20, 'HatchColor', 'black');
% hatchfill2(bars(5), 'single', 'HatchAngle', -45, 'HatchDensity', 40, 'HatchColor', 'black');
% hatchfill2(bars(6), 'single', 'HatchAngle', -45, 'HatchDensity', 20, 'HatchColor', 'black');
% legendData = {'CAMEO','Baseline'};
% [legend_h, object_h, plot_h, text_str] = legendflex(bars, legendData, 'FontSize', fontsize, ...
%     'anchor', [3 3], 'buffer', [-10 -10],'ncol',1);
% apply hatch pattern to legends
% hatchfill2(object_h(2+1), 'single', 'HatchAngle', 45, 'HatchDensity', 40/3, 'HatchColor', 'black');
% hatchfill2(object_h(6+2), 'single', 'HatchAngle', 45, 'HatchDensity', 20/3, 'HatchColor', 'black');
% hatchfill2(object_h(2+2), 'cross', 'HatchAngle', 45, 'HatchDensity', 40/3, 'HatchColor', 'black');
% hatchfill2(object_h(6+4), 'cross', 'HatchAngle', 45, 'HatchDensity', 20/3, 'HatchColor', 'black');
% hatchfill2(object_h(6+5), 'single', 'HatchAngle', -45, 'HatchDensity', 40/3, 'HatchColor', 'black');
% hatchfill2(object_h(6+6), 'single', 'HatchAngle', -45, 'HatchDensity', 20/3, 'HatchColor', 'black');

grid on;
ylabel("SNR (dB)");
set(gca,'xtick', [1,2]);
set(gca,'xticklabel', [sysname,"Baseline"]);
set(gca, 'FontSize', fontsize);
% set(legend_h,'fontsize',13);
set(gca, 'XMinorTick','on', 'XMinorGrid','on', 'YMinorTick','on', 'YMinorGrid','on');
exportgraphics(fig,"./figures/ofdm_snr.pdf",'Resolution',300);