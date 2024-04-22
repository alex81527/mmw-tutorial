clear;
load("./coal_plot_scripts/sim/data2.mat");
size(data2.sv)
pa = get_phased_array(60.48e9);

az = [-55:55];
gtx = steervec(pa.getElementPosition()/(physconst('LightSpeed')/60.48e9), ...
    [az ; zeros(1,length(az))]);

dot = gtx'*data2.sv;
[M, I] = max(abs(dot),[],1);
% figure;  plot(az(I))

% 2.
[N2,edges2] = histcounts(az(I),length(az));
% figure; plot(edges2(1:end-1), N2); 

cutoff = 12;
Nsamples = [100:100:1000];
min_cnt = min(N2(N2>=7));

new_idx = zeros(length(Nsamples), 1000);
L = zeros(length(Nsamples), 1);

for jj=1:length(Nsamples)
    y = az(I(1:Nsamples(jj)));
    tmp = [];
    for ii=1:length(N2)
        t = find(y>edges2(ii) & y<edges2(ii+1));
%         if length(t) >= min_cnt
%             tmp = [tmp t(1:min_cnt)];
%         end
        tmp = [tmp t(1:min(cutoff, length(t)))];
    end
    L(jj) = length(tmp);
    new_idx(jj, 1:L(jj)) = tmp;
end

% [N2,edges2] = histcounts(az(I(new_idx)),length(az));
% figure; plot(edges2(1:end-1), N2); 
for jj=1:length(Nsamples)
    if L(jj)>0
        fprintf("%d samples, mean az: %.2f vs %.2f\n", Nsamples(jj), mean(data2.az(1:Nsamples(jj))), mean(data2.az(I(new_idx(jj, 1:L(jj))))));
    end
end


A = zeros(32, length(data2.az));
for ii=1:length(data2.az)
    A(:,ii) = conj(steervec(pa.getElementPosition()/(physconst('LightSpeed')/60.48e9), ...
        [data2.az(ii); data2.el(ii)]));
end

rng(0);
path(path, "C:\Users\Weihan\Documents\GitHub\mmw-tutorial\coal_plot_scripts\sim");
[res,  az_err, el_err] = helper(data2, randperm(length(data2.az)), A);
% [res2,  az_err2, el_err2] = helper(data2, new_idx(randperm(length(new_idx))), A);

az_err2 = zeros(1, length(Nsamples));
res2 = zeros(32, length(Nsamples));
for jj=1:length(Nsamples)
    [res2_,  az_err2_, el_err2_] = helper(data2, I(new_idx(jj, 1:L(jj))), A);
    az_err2(jj) = az_err2_(end);
    res2(:,jj) = res2_(:,end);
end
figure;
plot(100*[1:length(az_err)], abs(az_err), 'bo-'); hold on;
plot(100*[1:length(az_err2)], abs(az_err2), 'rs-'); 
%
load("cal32_new_taoffice.mat"); 
gnd_truth = conj(exp(1j*calibration_vec.')./exp(1j*calibration_vec(26)));
% PA.PHASE_CAL(1:32) = calibration_vec;
figure;
plot(angle(conj(gnd_truth))); hold on;
plot(angle(res(:,end))); hold on;
plot(angle(res2(:,end))); hold on;

% for ii=1:size(res2,2)
%     cal_refant_26 = angle(res2(:,ii));
%     save(sprintf("./mat_files/cal32_%d_with_selection.mat", ii*100), "cal_refant_26");
% end

% plot([1:size(res, 2)]*100, rms(angle(res./gnd_truth), 1)); hold on;
% plot([1:size(res2, 2)]*100, rms(angle(res2./gnd_truth), 1)); hold on;