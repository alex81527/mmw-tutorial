path(path, "./plot_scripts");
scenario = "NLOS";  %"NLOS"
if scenario == "LOS"
    datasize = 49;
    save_folder = "./exp_data/swift2.0/los_multipath";
    good_idx = setdiff([1:datasize],[3,10,11,12,13,18,25,42]); % LOS dataset 51+
elseif scenario == "NLOS"
    datasize = 30;
    save_folder = "./exp_data/swift2.0/nlos_multipath";
    good_idx = setdiff([1:datasize],[1]); % NLOS dataset 51+
end
offset = 50;
load(sprintf("%s/data%d.mat",save_folder,offset+1));
snr_cs_multipath    = zeros(datasize, length(data.aco_cbsize));
snr_cs_multipath_v2 = zeros(datasize, length(data.aco_cbsize));
snr_cs_dominantpath = zeros(datasize, length(data.aco_cbsize));
snr_aco             = zeros(datasize, length(data.aco_cbsize));
snr_11ad            = zeros(datasize, length(data.aco_cbsize));
tpt_cs_multipath    = zeros(datasize, length(data.aco_cbsize));
tpt_cs_multipath_v2 = zeros(datasize, length(data.aco_cbsize));
tpt_cs_dominantpath = zeros(datasize, length(data.aco_cbsize));
tpt_aco             = zeros(datasize, length(data.aco_cbsize));
tpt_11ad            = zeros(datasize, length(data.aco_cbsize));
for ii=1:datasize
    load(sprintf("%s/data%d.mat",save_folder, offset+ii));
    snr_cs_multipath(ii,:)    = data.snr_cs_multipath;
    snr_cs_multipath_v2(ii,:) = data.snr_cs_multipath_v2;
    snr_cs_dominantpath(ii,:) = data.snr_cs_dominantpath;
    snr_aco(ii,:)             = data.snr_aco;
    snr_11ad(ii,:)            = data.snr_11ad;

    tpt_cs_multipath(ii,:)    = data.tpt_cs_multipath;
    tpt_cs_multipath_v2(ii,:) = data.tpt_cs_multipath_v2;
    tpt_cs_dominantpath(ii,:) = data.tpt_cs_dominantpath;
    tpt_aco(ii,:)             = data.tpt_aco;
    tpt_11ad(ii,:)            = data.tpt_11ad;
    fprintf("#%d, bpu tx: %d, bpu rx: %d, pa tx: %d, pa rx: %d, snr: %d\n", ...
        ii,data.bpu.TX_RF_GAIN, data.bpu.RX_RF_GAIN, data.pa.TX_IF_GAIN, data.pa.RX_IF_GAIN, round(data.snr_cs_multipath(end)));
end

%% snr cdf
leg=["CAMEO", "ACO", "11ad"]; 
linestyles = ["-", "--", "-.", ":",'-']; 
if scenario == "LOS"
    export_fname = "./figures/los_snr_cdf.pdf";
elseif scenario == "NLOS"
    export_fname = "./figures/nlos_snr_cdf.pdf";
end
% good_idx = [1:datasize];
% size(good_idx)
close all;
fig_size = [1 1 6 4];
fontsize = 22;
p = [snr_cs_multipath(good_idx,3) snr_aco(good_idx,end) snr_11ad(good_idx,end)]; 
my_cdfplot(p,"SNR (dB)",leg,linestyles,[],fontsize,fig_size,export_fname);

% save("./plot_scripts/fig9c.mat", "p");
% my_lineplot(good_idx,[ p(:,1)-p(:,3)],"data #",...
%     "SNR over ACO (dB)",["SWIFT2.0"],[]);
% 
% my_lineplot(good_idx,[p(:,1)],"data #", "SNR (dB)",["SWIFT2.0"],[]);
%% Throughput cdf
leg=["CAMEO", "ACO", "11ad"]; 
linestyles = ["-", "--", "-.", ":",'-']; 
if scenario == "LOS"
    export_fname = "./figures/los_tpt_cdf.pdf";
elseif scenario == "NLOS"
    export_fname = "./figures/nlos_tpt_cdf.pdf";
end
% p = [tpt_cs_multipath(good_idx,end) tpt_cs_multipath_v2(good_idx,end) tpt_aco(good_idx,end) tpt_11ad(good_idx,end)]; 
% my_cdfplot(p(:,[1 3 4])./1e6,"Throughput (Mbps)",leg([1 3 4]),linestyles,[],"./figures/tpt_pdf.pdf");
close all;
fig_size = [1 1 6 4];
fontsize = 22;
p = [snr_cs_multipath(good_idx,3) snr_aco(good_idx,end) snr_11ad(good_idx,end)]; 
[rate, mcs]= get_11ad_datarate(p);
rate = rate/(1760e6/160e6);
my_cdfplot(rate./1e6,"Throughput (Mbps)",leg,linestyles,[],fontsize,fig_size,export_fname);

% p=rate;
% save("./plot_scripts/fig9d.mat", "p");

%% Required # probes vs SNR
% good_idx = [1:datasize];
% good_idx = setdiff([1:datasize],[3,10,11,12,13,18,25,42]); % LOS dataset 51+
good_idx = setdiff([1:datasize],[1]); % NLOS dataset 51+
num_probe = data.aco_cbsize;
snr_point = [3:20];

mean_snr1 = median(snr_cs_multipath(good_idx,:));
mean_snr2 = median(snr_aco(good_idx,:));
mean_snr3 = median(snr_11ad(good_idx,:));
ppp = zeros(length(snr_point),3);
for ii=1:length(snr_point)
    x1 = find(mean_snr1 > snr_point(ii), 1);
    x2 = find(mean_snr2 > snr_point(ii), 1);
    x3 = find(mean_snr3 > snr_point(ii), 1);

    if ~isempty(x1)
        ppp(ii,1) = num_probe(x1);
    else
        ppp(ii,1) = num_probe(end);
    end
    if ~isempty(x2)
        ppp(ii,2) = num_probe(x2);
    else
        ppp(ii,2) = num_probe(end);
    end
    if ~isempty(x3)
        ppp(ii,3) = num_probe(x3);
    else
        ppp(ii,3) = num_probe(end);
    end
end
leg=["SWIFT2.0","ACO", "11ad"]; 
linestyles = ["-", "--", "-.", ":",'-']; 
export_fname = [];
my_lineplot(snr_point,ppp,"SNR", "Required number of probes",leg,export_fname);

% figure;
% plot(num_probe, mean(snr_cs_multipath(good_idx,:)).'); hold on;
% plot(num_probe, mean(snr_aco(good_idx,:)).');
% figure; plot(good_idx, snr_aco(good_idx,1));
% g2 = setdiff(good_idx, [1 2 9 24 30 31 44 45 49]);
% g2 = setdiff(good_idx, [find(snr_aco(good_idx,1)>snr_aco(good_idx,4))]);
% g2 = find(snr_aco(good_idx,end)<=28 & snr_aco(good_idx,1)<=15);
% figure;
% plot(max(snr_aco(good_idx,:),[],2) - min(snr_aco(good_idx,:),[],2)); 
% plot(max(snr_11ad(good_idx,:),[],2) - min(snr_11ad(good_idx,:),[],2)); 


% snr_point = [10,15,20,25,30];
% pp = -1*ones(length(good_idx), length(snr_point),3);
% 
% for ii=1:length(good_idx)
%     for jj=1:length(snr_point)
%         x1 = find(snr_cs_multipath(good_idx(ii),:)>snr_point(jj), 1);
%         x2 = find(snr_aco(good_idx(ii),:)>snr_point(jj), 1);
%         x3 = find(snr_11ad(good_idx(ii),:)>snr_point(jj), 1);
%         if ~isempty(x1)
%             pp(ii,jj,1) = num_probe(x1);
%         end
%         if ~isempty(x2)
%             pp(ii,jj,2) = num_probe(x2);
%         end
%         if ~isempty(x3)
%             pp(ii,jj,3) = num_probe(x3);
%         end
%     end
% end
% 
% ppp = zeros(length(snr_point),3);
% for ii=1:length(snr_point)
%     valid_idx = find(pp(:,ii,1)>0);
%     ppp(ii,1) = mean(pp(valid_idx,ii,1));
% 
%     valid_idx = find(pp(:,ii,2)>0);
%     ppp(ii,2) = mean(pp(valid_idx,ii,2));
% 
%     valid_idx = find(pp(:,ii,3)>0);
%     ppp(ii,3) = mean(pp(valid_idx,ii,3));
% end
% 
% leg=["SWIFT2.0","ACO", "11ad"]; 
% linestyles = ["-", "--", "-.", ":",'-']; 
% export_fname = [];
% my_lineplot(snr_point,ppp,"SNR", "# probes",leg,export_fname);