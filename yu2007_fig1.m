function out = yu2007_fig1()
%YU2007_FIG1  Yu et al. (2007) Fig.1：距离、背景风、水平切变、反射率结构对龙卷谱的影响
%
%   四个面板：
%     (a) r0 = 5, 25, 75 km
%     (b) 背景风 v0 = 0, 5, 15 m/s
%     (c) 水平切变 dv/dx = 0, 0.02, 0.08 s^-1
%     (d) 反射率宽度 Wz = 20, 60, inf m
%
%   默认：Combined Rankine 龙卷，VTm=50 m/s，rt=200 m，r0z=220 m，
%   va=65 m/s，1° 波束，250 m 距离分辨率。

here = fileparts(mfilename('fullpath'));
addpath(here);
figdir = fullfile(here,'figs'); datadir = fullfile(here,'data');
if ~exist(figdir,'dir'), mkdir(figdir); end
if ~exist(datadir,'dir'), mkdir(datadir); end
rng(2007,'twister');

va = 65; M = 128;
base = struct('va',va,'M',M,'beamwidth',1,'sigma_r',87.5,'gridN',301, ...
    'vel_jitter',0.5,'smooth_bins',1.0,'u0',0,'v0',0,'shear_dvdx',0, ...
    'vortices',struct('x',0,'y',0,'rt',200,'vtm',50,'vrm',0), ...
    'refl',struct('type','doughnut','Wz',60,'r0z',220,'zm',1));

cases = {};
tags = {};
% (a) range
for r = [5 25 75]
    o = base; o.r0 = r*1000; o.vortices.y = o.r0;
    cases{end+1} = o; tags{end+1} = sprintf('r_0 = %d km',r);
end
% (b) background wind
for v0 = [0 5 15]
    o = base; o.r0 = 25000; o.vortices.y = o.r0; o.v0 = v0;
    cases{end+1} = o; tags{end+1} = sprintf('v_0 = %d m s^{-1}',v0);
end
% (c) horizontal shear
for sh = [0 0.02 0.08]
    o = base; o.r0 = 25000; o.vortices.y = o.r0; o.shear_dvdx = sh;
    cases{end+1} = o; tags{end+1} = sprintf('dv/dx = %.2f s^{-1}',sh);
end
% (d) reflectivity width：Wz=20 用高斯型（实心核，钟形谱），Wz=60 用
% doughnut（环状，平顶+中心凹陷），Wz=inf 为均匀反射率。论文未明确列出
% r0z 的切换，这是按 Fig.1d 三条曲线形态反推的等效设置。
for wz = [20 60 1e5]
    o = base; o.r0 = 25000; o.vortices.y = o.r0;
    if wz < 40
        o.refl = struct('type','gaussian','Wz',wz,'r0z',0,'zm',1);
    elseif wz < 1000
        o.refl = struct('type','doughnut','Wz',wz,'r0z',220,'zm',1);
    else
        o.refl = struct('type','uniform','Wz',wz,'r0z',0,'zm',1);
    end
    cases{end+1} = o;
    if wz > 1000, tags{end+1} = 'W_z = \infty'; else, tags{end+1} = sprintf('W_z = %d m',wz); end
end

Sdb = cell(1,12); vv = cell(1,12);
for k = 1:12
    [S,v] = y07_sim_volume(cases{k});
    S = S/max(S);
    Sdb{k} = 10*log10(max(S,10^(-75/10))) + 25;   % 峰值 25 dB，噪底 -50 dB
    vv{k} = v;
end

fig = figure('Color','w','Position',[80 80 860 700]);
tiledlayout(2,2,'Padding','compact','TileSpacing','compact');
styles = {'-','--','-.'};
colors = {[0 0 0],[0.15 0.15 0.15],[0.35 0.35 0.35]};
labs = {'(a) Range','(b) Background wind','(c) Horizontal shear','(d) Reflectivity width'};
for p = 1:4
    nexttile; hold on; box on; grid on;
    for j = 1:3
        k = (p-1)*3 + j;
        plot(vv{k},Sdb{k},styles{j},'Color',colors{j},'LineWidth',1.3);
    end
    xlim([-va va]); ylim([-55 35]);
    xlabel('Radial velocity (m s^{-1})'); ylabel('Doppler spectrum (dB)');
    title(labs{p});
    legend(tags((p-1)*3+(1:3)),'Location','southeast','FontSize',7); legend boxoff;
end
sgtitle('Yu et al. (2007) Fig. 1 reproduction','FontWeight','bold');
saveas(fig,fullfile(figdir,'yu2007_fig1.png'));

out = struct('Sdb',{Sdb},'v',{vv},'tags',{tags});
save(fullfile(datadir,'yu2007_fig1.mat'),'out');
end
