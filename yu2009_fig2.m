function yu2009_fig2()
%YU2009_FIG2  Yu et al. (2009) Fig.2：观测谱、自协方差重构谱与双高斯拟合谱
%
%   四组实验（列）× 三个速度间隔 dv=-18,-14,-2 m/s（行）。
%   每个面板：浅灰粗线=含噪观测谱；蓝色虚线=自协方差法重构；红色实线=双高斯拟合。

here = fileparts(mfilename('fullpath'));
addpath(here);
figdir = fullfile(here,'figs'); datadir = fullfile(here,'data');
if ~exist(figdir,'dir'), mkdir(figdir); end
if ~exist(datadir,'dir'), mkdir(datadir); end

S1 = 10; v1 = 10; s1 = 2; va = 25;
expt = struct('name',{'I','II','III','IV'}, 'S2db',{10,10,1,-20}, ...
              's2',{2,6,6,4}, 'color',{[0.85 0.1 0.1],[0.1 0.25 0.85],[0.05 0.55 0.15],[0 0 0]});
dvs = [-18 -14 -2];
ndv = numel(dvs); ne = 4;

F = cell(ne,ndv);
fig = figure('Color','w','Position',[60 60 1150 780]);
tiledlayout(ndv,ne,'Padding','compact','TileSpacing','compact');
for e = 1:ne
    S2 = 10^(expt(e).S2db/10);
    for k = 1:ndv
        v2 = v1 + dvs(k);
        [Sp,v,info] = y09_sim_spectrum(S1,S2,v1,v2,s1,expt(e).s2, ...
                                       va,64,30,1,8,1000+100*e+k);
        fit = y09_fit_dual(Sp,v,va,struct('nStart',12));
        F{e,k} = struct('Sp',Sp,'v',v,'info',info,'fit',fit);
    end
end
save(fullfile(datadir,'yu2009_fig2.mat'),'F','expt','dvs');

for e = 1:ne
    for k = 1:ndv
        nexttile((k-1)*ne+e); hold on; box on; grid on;
        d = F{e,k};
        Pnpsd = d.info.Pn/(2*va);
        Pac = max(sum(d.Sp)*d.info.dv - d.info.Pn, 1e-9);
        Sac = Pac*exp(-(d.v-d.info.vac).^2/(2*d.info.sac^2))/(sqrt(2*pi)*d.info.sac) + Pnpsd;
        yobs = 10*log10(max(d.Sp,eps));
        plot(d.v,yobs,'-','Color',[0.65 0.65 0.65],'LineWidth',2.2);
        plot(d.v,10*log10(max(Sac,eps)),'--','Color',[0.1 0.3 0.9],'LineWidth',1.2);
        plot(d.v,d.fit.yfit,'-','Color',[0.85 0.1 0.1],'LineWidth',1.4);
        xlim([-va va]); ylim([-45 25]);
        if k == 1
            title(sprintf('Expt %s',expt(e).name));
        end
        if e == 1, ylabel(sprintf('\\Delta v=%d',dvs(k))); end
        if k == ndv, xlabel('v_r (m s^{-1})'); end
        text(-va+2, 18, sprintf('sk=%.2f\nkt=%.2f',d.info.sk,d.info.kt), ...
             'FontSize',7,'VerticalAlignment','top');
        if e==1 && k==1
            legend({'Observed','Autocov.','Dual-Gaussian'},'Location','northeast','FontSize',7);
            legend boxoff;
        end
    end
end
sgtitle('Yu et al. (2009) Fig. 2 reproduction','FontWeight','bold');
saveas(fig,fullfile(figdir,'yu2009_fig2.png'));
end
