function yu2009_fig1(doSim)
%YU2009_FIG1  Yu et al. (2009) Fig.1：双高斯谱的自协方差法理论偏差与数值验证
%
%   yu2009_fig1()      % 理论 + 200 次实现数值实验
%   yu2009_fig1(false) % 仅理论曲线
%
%   四组实验（Table 1）：
%     I   S2=10 dBm, s2=2 m/s
%     II  S2=10 dBm, s2=6 m/s
%     III S2= 1 dBm, s2=6 m/s
%     IV  S2=-20 dBm, s2=4 m/s
%   公共参数：S1=10 dBm, v1=10 m/s, s1=2 m/s, va=25 m/s，v2=-10:1:10。

if nargin < 1, doSim = true; end
here = fileparts(mfilename('fullpath'));
addpath(here);
figdir = fullfile(here,'figs'); datadir = fullfile(here,'data');
if ~exist(figdir,'dir'), mkdir(figdir); end
if ~exist(datadir,'dir'), mkdir(datadir); end

S1 = 10; v1 = 10; s1 = 2; va = 25;
expt = struct( ...
    'name', {'I','II','III','IV'}, ...
    'S2db', {10, 10, 1, -20}, ...
    's2',   {2, 6, 6, 4}, ...
    'color',{[0.85 0.1 0.1],[0.1 0.25 0.85],[0.05 0.55 0.15],[0 0 0]});
v2 = -10:1:10;
dv = v2 - v1;
ndv = numel(dv);

res = struct();
for e = 1:4
    S2 = 10^(expt(e).S2db/10);
    th_v = zeros(1,ndv); th_s = zeros(1,ndv);
    th_sk = zeros(1,ndv); th_kt = zeros(1,ndv);
    bias_v = zeros(1,ndv); bias_s = zeros(1,ndv);
    thb_v = zeros(1,ndv); thb_s = zeros(1,ndv);
    sd_v = zeros(1,ndv); sd_s = zeros(1,ndv);
    for k = 1:ndv
        o = y09_theory_ac(S1,S2,v1,v2(k),s1,expt(e).s2,va);
        th_v(k) = o.vr; th_s(k) = o.sv;
        th_sk(k) = o.sk; th_kt(k) = o.kt;
        thb_v(k) = o.bias_v; thb_s(k) = o.bias_sv;
        if doSim
            sim = y09_sim_bias(S1,S2,v1,v2(k),s1,expt(e).s2, ...
                               va,64,30,200,8);
            bias_v(k) = mean(sim.bias_v);
            bias_s(k) = mean(sim.bias_sv);
            sd_v(k)   = std(sim.bias_v);
            sd_s(k)   = std(sim.bias_sv);
        end
    end
    res(e).name = expt(e).name;
    res(e).dv = dv; res(e).th_v = th_v; res(e).th_s = th_s;
    res(e).th_sk = th_sk; res(e).th_kt = th_kt;
    res(e).bias_v = bias_v; res(e).bias_s = bias_s;
    res(e).thb_v = thb_v; res(e).thb_s = thb_s;
    res(e).sd_v = sd_v; res(e).sd_s = sd_s;
    res(e).color = expt(e).color;
    fprintf('Expt %s: max|bias_v|=%.2f m/s, bias_s range [%.2f, %.2f] m/s\n', ...
            expt(e).name, max(abs(bias_v)), min(bias_s), max(bias_s));
end

save(fullfile(datadir,'yu2009_fig1.mat'),'res','expt','S1','v1','s1','va','v2');

% ---------------- 绘图 ----------------
fig = figure('Color','w','Position',[80 80 980 760]);
tiledlayout(2,2,'Padding','compact','TileSpacing','compact');

nexttile; hold on; box on; grid on;
for e=1:4
    plot(dv,res(e).th_v,'-','Color',res(e).color,'LineWidth',1.3);
end
xlabel('\Delta v_r (m s^{-1})'); ylabel('\langle v_r \rangle (m s^{-1})');
xlim([-20 0]); title('(a) Theoretical mean velocity');
lgd = arrayfun(@(e)sprintf('Expt %s',res(e).name),1:4,'UniformOutput',false);
legend(lgd,'Location','southeast'); legend boxoff;

nexttile; hold on; box on; grid on;
for e=1:4
    plot(dv,res(e).th_s,'-','Color',res(e).color,'LineWidth',1.3);
end
xlabel('\Delta v_r (m s^{-1})'); ylabel('\sigma_v (m s^{-1})');
xlim([-20 0]); title('(b) Theoretical spectrum width');
legend(lgd,'Location','northwest'); legend boxoff;

nexttile; hold on; box on; grid on;
for e=1:4
    plot(res(e).th_sk,res(e).thb_v,'-','Color',res(e).color,'LineWidth',1.2);
    if doSim
        plot(res(e).th_sk,res(e).bias_v,'o','Color',res(e).color, ...
             'MarkerFaceColor',res(e).color,'MarkerSize',4);
    end
end
xlabel('Skewness'); ylabel('v_r bias (m s^{-1})');
title('(c) Velocity bias vs skewness');
xline(0,':','Color',[0.4 0.4 0.4]);

nexttile; hold on; box on; grid on;
for e=1:4
    plot(res(e).th_kt,res(e).thb_s,'-','Color',res(e).color,'LineWidth',1.2);
    if doSim
        plot(res(e).th_kt,res(e).bias_s,'o','Color',res(e).color, ...
             'MarkerFaceColor',res(e).color,'MarkerSize',4);
    end
end
xlabel('Kurtosis'); ylabel('\sigma_v bias (m s^{-1})');
title('(d) Width bias vs kurtosis');
xline(3,':','Color',[0.4 0.4 0.4]);

sgtitle('Yu et al. (2009) Fig. 1 reproduction','FontWeight','bold');
saveas(fig,fullfile(figdir,'yu2009_fig1.png'));
end
