function yu2007_fig5()
%YU2007_FIG5  Yu et al. (2007) Fig.5：9 个分辨率体积内的 Td/Tu/Gd 谱
%
%   龙卷位于 (r0'=60 km, phi0'=0)，体积中心在方位 -0.4/-0.2/0°、
%   距离 r0'-120/0/+120 m。Td=环状反射率(Wz=60)，Tu=均匀反射率，
%   Gd=由 Td 的自协方差矩重构的高斯谱。SNR=40 dB。

here = fileparts(mfilename('fullpath'));
addpath(here);
figdir = fullfile(here,'figs'); datadir = fullfile(here,'data');
if ~exist(figdir,'dir'), mkdir(figdir); end
if ~exist(datadir,'dir'), mkdir(datadir); end
rng(2007,'twister');

va = 35; M = 64; r0p = 60000;
azList = [-0.4 -0.2 0.0];
drList = [-120 0 120];
base = struct('va',va,'M',M,'beamwidth',1,'sigma_r',87.5,'gridN',301, ...
    'vel_jitter',2.0,'smooth_bins',2.0,'r0',r0p, ...
    'vortices',struct('x',0,'y',r0p,'rt',200,'vtm',50,'vrm',0));

fig = figure('Color','w','Position',[50 40 1000 900]);
tiledlayout(3,3,'Padding','tight','TileSpacing','tight');
Sd = cell(3,3); Su = cell(3,3); G = cell(3,3); vax=[];
for ir = 1:3
    for ia = 1:3
        o = base;
        o.az0 = azList(ia); o.r0 = r0p + drList(ir);
        % 龙卷固定在 (60 km, 0)，体积中心由 r0/az0 决定
        o.vortices.x = 0; o.vortices.y = r0p;

        o.refl = struct('type','doughnut','Wz',60,'r0z',220,'zm',1);
        [S,v] = y07_sim_volume(o);
        Sn = S/max(S);
        % 论文 Fig.5 的谱有一定起伏；用 2 个独立 trial 平均
        iq = ysim_iq_from_psd(Sn, M, 40, 2, 50000+100*ir+ia);
        [Sb1, vb] = ysim_periodogram(iq(:,1), va, M, 'hann');
        [Sb2, ~]  = ysim_periodogram(iq(:,2), va, M, 'hann');
        Sb = (Sb1+Sb2)/2;
        Sd{ir,ia} = Sb/max(Sb); vax=vb;
        am = ysim_ac_moments(iq(:,1), va);

        o.refl = struct('type','uniform','Wz',1e6,'r0z',0,'zm',1);
        [S,v] = y07_sim_volume(o);
        Sn = S/max(S);
        iq2 = ysim_iq_from_psd(Sn, M, 40, 2, 60000+100*ir+ia);
        [Sb21, ~] = ysim_periodogram(iq2(:,1), va, M, 'hann');
        [Sb22, ~] = ysim_periodogram(iq2(:,2), va, M, 'hann');
        Sb2 = (Sb21+Sb22)/2;
        Su{ir,ia} = Sb2/max(Sb2);

        % 由 Td 的 AC 矩重构高斯谱（解析重构，不再叠加新噪声）
        sg = max(am.sv,0.3);
        Sg = exp(-(v-am.v).^2/(2*sg^2))/(sqrt(2*pi)*sg);
        G{ir,ia} = Sg/max(Sg);

        nexttile((ir-1)*3+ia); hold on; box on;
        plot(vb,10*log10(max(Sd{ir,ia},10^(-80/10)))+25,'k-','LineWidth',1.1);
        plot(vb,10*log10(max(Su{ir,ia},10^(-80/10)))+25,'k:','LineWidth',1.1);
        plot(vb,10*log10(max(G{ir,ia},10^(-80/10)))+25,'k-.','LineWidth',1.0);
        xlim([-va va]); ylim([-40 45]);
        set(gca,'XTick',-30:30:30,'YTick',-40:20:40,'FontSize',7);
        if ir==1, title(sprintf('\\phi_0=%.1f^o',azList(ia)),'FontSize',9); end
        if ia==1, ylabel(sprintf('%+d m',drList(ir)),'FontSize',8); end
        if ir==3, xlabel('v_r (m s^{-1})','FontSize',8); end
        if ia>1, set(gca,'YTickLabel',[]); end
        if ir<3, set(gca,'XTickLabel',[]); end
        if ir==3 && ia==3
            legend({'T_d','T_u','G_d'},'Location','southwest','FontSize',7); legend boxoff;
        end
    end
end
sgtitle('Yu et al. (2007) Fig. 5 reproduction','FontWeight','bold');
saveas(fig,fullfile(figdir,'yu2007_fig5.png'));
save(fullfile(datadir,'yu2007_fig5.mat'),'Sd','Su','G','vax','azList','drList');
end
