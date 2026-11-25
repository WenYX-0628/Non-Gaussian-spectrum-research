function yu2007_fig4()
%YU2007_FIG4  Yu et al. (2007) Fig.4：高斯/矩形图案的 PRIB 随谱宽变化
%
%   图案：64 个速度 bin，v in [-50,50] m/s；矩形的宽度取 sqrt(12)*sigma，
%   使其与高斯图案具有相同的二阶矩。PRIB 使用论文 3b 节给出的
%   |X(f)| 替代序列方法（避免偶对称图案相位退化）。
%
%   说明：论文未给出可直接复现的 PRIB 预处理细节，本图绝对数值与
%   Fig.4 不完全一致（见 README）；但高斯与矩形在中等谱宽下的分离
%   趋势是可复现的。

here = fileparts(mfilename('fullpath'));
addpath(here);
figdir = fullfile(here,'figs'); datadir = fullfile(here,'data');
if ~exist(figdir,'dir'), mkdir(figdir); end
if ~exist(datadir,'dir'), mkdir(datadir); end

N = 64; va = 50; v = linspace(-va,va,N);
sigs = 1:1:25;
Pg = zeros(size(sigs)); Pr = zeros(size(sigs));
for k = 1:numel(sigs)
    sg = sigs(k);
    g = exp(-v.^2/(2*sg^2));
    w = sqrt(12)*sg;
    r = double(abs(v) <= w/2);
    gdb = 10*log10(max(g,1e-6));
    rdb = 10*log10(max(r,1e-6));
    Pg(k) = abs(ysim_prib(gdb,0.9,256,'absfft'));
    Pr(k) = abs(ysim_prib(rdb,0.9,256,'absfft'));
end

fig = figure('Color','w','Position',[100 100 900 380]);
tiledlayout(1,2,'Padding','compact','TileSpacing','compact');

nexttile; hold on; box on; grid on;
sg0 = 12; g0 = exp(-v.^2/(2*sg0^2)); w0=sqrt(12)*sg0;
r0 = double(abs(v)<=w0/2);
plot(v,g0,'k-','LineWidth',1.5);
plot(v,r0,'k--','LineWidth',1.4);
xlim([-50 50]); ylim([0 1.05]);
xlabel('V_r (m s^{-1})'); ylabel('Amplitude');
title('(a) Gaussian and rectangular patterns (\sigma_v=12 m s^{-1})');
legend({'Gaussian','Rectangular'},'Location','northeast'); legend boxoff;

nexttile; hold on; box on; grid on;
plot(sigs,Pg,'ko-','LineWidth',1.1,'MarkerSize',3,'MarkerFaceColor','k');
plot(sigs,Pr,'ks--','LineWidth',1.1,'MarkerSize',3);
xlabel('\sigma_v (m s^{-1})'); ylabel('P (deg)');
xlim([0 25]); ylim([0 35]);
title('(b) PRIB vs spectrum width');
legend({'Gaussian','Rectangular'},'Location','northeast'); legend boxoff;

sgtitle('Yu et al. (2007) Fig. 4 reproduction','FontWeight','bold');
saveas(fig,fullfile(figdir,'yu2007_fig4.png'));
save(fullfile(datadir,'yu2007_fig4.mat'),'sigs','Pg','Pr','N','va');
end
