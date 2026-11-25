function yu2007_fig678(nreal)
%YU2007_FIG678  Yu et al. (2007) Fig.6-8：谱宽、谱平坦度、PRIB 的距离廓线
%
%   yu2007_fig678()    % 每个点 50 次实现（论文设置）
%   yu2007_fig678(10)  % 快速测试
%
%   9 个面板 = 3 个方位偏移 (-0.4,-0.2,0 deg) x 3 个距离偏移
%   (-120,0,+120 m)。横轴 rn = r0'*theta_b/rt，0.25:0.5:10。
%   每个点对 Td（doughnut）和 Tu（uniform）各做 nreal 次独立实现，
%   并由各自的自协方差矩重构高斯谱 Gd/Gu，统计：
%     sigma_vn = sigma_v/(2 va)   （自协方差法）
%     sigma_s  = std(10log10 S)   （dB 谱标准差）
%     P        = PRIB 相位，按论文对 50 次实现取算术平均
%
%   注意：PRIB 单次实现的相位在对称谱上会在 0/180° 附近跳变；论文
%   Fig.8 的 50 次平均曲线用算术平均才能得到 20-70° 的量级（本实现
%   与之一致）。因此这里使用 direct 方法并做算术平均。

if nargin < 1 || isempty(nreal), nreal = 50; end
here = fileparts(mfilename('fullpath'));
addpath(here);
figdir = fullfile(here,'figs'); datadir = fullfile(here,'data');
if ~exist(figdir,'dir'), mkdir(figdir); end
if ~exist(datadir,'dir'), mkdir(datadir); end

va = 35; M = 64; rt = 200; thb = pi/180;
azList = [-0.4 -0.2 0.0];
drList = [-120 0 120];
rnList = 0.25:0.5:10;
nr = numel(rnList);

R = struct();
R.rn = rnList; R.azList=azList; R.drList=drList;
R.svn_td = nan(3,3,nr); R.svn_tu = nan(3,3,nr);
R.ss_td  = nan(3,3,nr); R.ss_tu  = nan(3,3,nr);
R.ss_gd  = nan(3,3,nr); R.ss_gu  = nan(3,3,nr);
R.p_td   = nan(3,3,nr); R.p_tu   = nan(3,3,nr);
R.p_gd   = nan(3,3,nr); R.p_gu   = nan(3,3,nr);

rng(2007,'twister');
tic;
for ia = 1:3
    for idr = 1:3
        for ir = 1:nr
            r0p = rnList(ir)*rt/thb;
            o = struct('va',va,'M',M,'beamwidth',1,'sigma_r',87.5, ...
                'gridN',241,'vel_jitter',1.0,'smooth_bins',1.5, ...
                'r0',r0p+drList(idr),'az0',azList(ia), ...
                'vortices',struct('x',0,'y',r0p,'rt',rt,'vtm',50,'vrm',0));
            acc = zeros(1,14);
            for typ = 1:2
                if typ==1
                    o.refl = struct('type','doughnut','Wz',60,'r0z',220,'zm',1);
                else
                    o.refl = struct('type','uniform','Wz',1e6,'r0z',0,'zm',1);
                end
                [S,v] = y07_sim_volume(o);
                Sn = S/max(S);
                for k = 1:nreal
                    iq = ysim_iq_from_psd(Sn, M, 40, 1);
                    [Sb,vb] = ysim_periodogram(iq, va, M, 'hann');
                    f = ysim_spectrum_features(Sb, vb, struct('method','direct'));
                    am = ysim_ac_moments(iq, va);
                    svn = am.sv/(2*va);

                    % 由 AC 矩重构的高斯谱：解析谱 + 同 SNR 的噪声底
                    % 用谱质心而不是混叠后的 AC 平均速度来放置重构高斯，
                    % 避免高斯被平移到 Nyquist 边缘造成虚假的大动态范围
                    sg = min(max(am.sv,10),25);
                    vc = f.vm;
                    Sg = exp(-(vb-vc).^2/(2*sg^2))/(sqrt(2*pi)*sg);
                    Pnpsd = (sum(Sg)*10^(-40/10))/(2*va);
                    fg = ysim_spectrum_features(Sg+Pnpsd, vb, struct('method','direct'));

                    if typ==1
                        acc = acc + [svn, f.sigma_s, f.PRIB, fg.sigma_s, fg.PRIB, 0,0,0,0,0,0,0,0,0];
                    else
                        acc = acc + [0,0,0,0,0, svn, f.sigma_s, f.PRIB, fg.sigma_s, fg.PRIB,0,0,0,0];
                    end
                end
            end
            acc = acc/nreal;
            R.svn_td(ia,idr,ir)=acc(1); R.ss_td(ia,idr,ir)=acc(2); R.p_td(ia,idr,ir)=acc(3);
            R.ss_gd(ia,idr,ir)=acc(4); R.p_gd(ia,idr,ir)=acc(5);
            R.svn_tu(ia,idr,ir)=acc(6); R.ss_tu(ia,idr,ir)=acc(7); R.p_tu(ia,idr,ir)=acc(8);
            R.ss_gu(ia,idr,ir)=acc(9); R.p_gu(ia,idr,ir)=acc(10);
        end
        fprintf('panel az=%.1f dr=%+d done (%.1f s)\n',azList(ia),drList(idr),toc);
    end
end
save(fullfile(datadir,'yu2007_fig678.mat'),'R','nreal');

% ---------------- 绘图 ----------------
plot_profiles(R,'svn','\sigma_{vn}',[0 0.4],'yu2007_fig6.png','Yu et al. (2007) Fig. 6 reproduction');
plot_profiles(R,'ss','\sigma_s',[0 15],'yu2007_fig7.png','Yu et al. (2007) Fig. 7 reproduction');
plot_profiles(R,'p','P (deg)',[0 80],'yu2007_fig8.png','Yu et al. (2007) Fig. 8 reproduction');
end

function plot_profiles(R,kind,ylab,ylimv,fname,ttl)
fig = figure('Color','w','Position',[60 40 900 760]);
tiledlayout(3,3,'Padding','compact','TileSpacing','compact');
cols = {[0 0 0],[0.45 0.45 0.45],[0.1 0.3 0.9],[0.85 0.1 0.1]};
for ia=1:3
    for idr=1:3
        nexttile; hold on; box on; grid on;
        rn = R.rn;
        switch kind
            case 'svn'
                plot(rn,squeeze(R.svn_td(ia,idr,:)),'-','Color',cols{1},'LineWidth',1.4);
                plot(rn,squeeze(R.svn_tu(ia,idr,:)),'--','Color',cols{2},'LineWidth',1.4);
            case 'ss'
                plot(rn,squeeze(R.ss_td(ia,idr,:)),'-','Color',cols{1},'LineWidth',1.4);
                plot(rn,squeeze(R.ss_gd(ia,idr,:)),'-.','Color',cols{3},'LineWidth',1.2);
                plot(rn,squeeze(R.ss_tu(ia,idr,:)),':','Color',cols{2},'LineWidth',1.4);
                plot(rn,squeeze(R.ss_gu(ia,idr,:)),'--','Color',cols{4},'LineWidth',1.2);
            case 'p'
                plot(rn,squeeze(R.p_td(ia,idr,:)),'-','Color',cols{1},'LineWidth',1.4);
                plot(rn,squeeze(R.p_gd(ia,idr,:)),'-.','Color',cols{3},'LineWidth',1.2);
                plot(rn,squeeze(R.p_tu(ia,idr,:)),':','Color',cols{2},'LineWidth',1.4);
                plot(rn,squeeze(R.p_gu(ia,idr,:)),'--','Color',cols{4},'LineWidth',1.2);
        end
        xlim([0 10]); ylim(ylimv);
        set(gca,'XTick',0:2:10,'FontSize',7);
        if ia==1, title(sprintf('\\phi_0=%.1f^o',R.azList(ia)),'FontSize',8); end
        if idr==1, ylabel(sprintf('%+d m',R.drList(idr)),'FontSize',8); end
        if ia==3, xlabel('r_n = r_0''\theta_b/r_t','FontSize',7); end
        if idr>1, set(gca,'YTickLabel',[]); end
        if ia<3, set(gca,'XTickLabel',[]); end
        if ia==3 && idr==3
            if strcmp(kind,'svn')
                legend({'T_d','T_u'},'Location','best','FontSize',6); legend boxoff;
            else
                legend({'T_d','G_d','T_u','G_u'},'Location','best','FontSize',6); legend boxoff;
            end
        end
    end
end
sgtitle(ttl,'FontWeight','bold');
saveas(fig,fullfile(fileparts(mfilename('fullpath')),'figs',fname));
end
