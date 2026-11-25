function yu2009_fig3(nreal)
%YU2009_FIG3  Yu et al. (2009) Fig.3：两种重构方法的 RMSE 概率密度
%
%   yu2009_fig3()    % 默认每个 dv 30 次实现
%   yu2009_fig3(20)
%
%   对四组实验、dv=-20:1:0 的每个间隔，重复 nreal 次独立实现，比较
%   - 自协方差法 + 高斯重构
%   - 双高斯 dB 域拟合
%   的 RMSE（dB）。计算量较大，若 Parallel Computing Toolbox 可用则并行。

if nargin < 1 || isempty(nreal), nreal = 30; end
here = fileparts(mfilename('fullpath'));
addpath(here);
figdir = fullfile(here,'figs'); datadir = fullfile(here,'data');
if ~exist(figdir,'dir'), mkdir(figdir); end
if ~exist(datadir,'dir'), mkdir(datadir); end

S1 = 10; v1 = 10; s1 = 2; va = 25; M = 64; snr_db = 30;
expt = struct('name',{'I','II','III','IV'}, 'S2db',{10,10,1,-20}, ...
              's2',{2,6,6,4}, 'color',{[0.85 0.1 0.1],[0.1 0.25 0.85],[0.05 0.55 0.15],[0 0 0]});
v2list = -10:1:10;
ndv = numel(v2list);

usePar = false;
try
    if license('test','Distrib_Computing_Toolbox')
        p = gcp('nocreate');
        if isempty(p), p = parpool('local',4); end %#ok<NASGU>
        usePar = true;
    end
catch
    usePar = false;
end
fprintf('Parallel: %d\n',usePar);

EAC = cell(1,4); EDG = cell(1,4);
for e = 1:4
    S2 = 10^(expt(e).S2db/10);
    eac = zeros(nreal,ndv); edg = zeros(nreal,ndv);
    for k = 1:ndv
        v2 = v2list(k);
        if usePar
            parfor r = 1:nreal
                [eac(r,k), edg(r,k)] = y09_rmse_one(S1,S2,v1,v2,s1,expt(e).s2, ...
                    va,M,snr_db, 700000 + 10000*e + 100*k + r, 8);
            end
        else
            for r = 1:nreal
                [eac(r,k), edg(r,k)] = y09_rmse_one(S1,S2,v1,v2,s1,expt(e).s2, ...
                    va,M,snr_db, 700000 + 10000*e + 100*k + r, 8);
            end
        end
    end
    EAC{e} = eac; EDG{e} = edg;
    fprintf('Expt %s: mean RMSE AC=%.2f dB, DG=%.2f dB\n', ...
            expt(e).name, mean(eac(:)), mean(edg(:)));
end
save(fullfile(datadir,'yu2009_fig3.mat'),'EAC','EDG','expt','v2list','nreal');

% ---------------- PDF ----------------
edges = 0:0.5:20;
ctr = edges(1:end-1)+diff(edges)/2;
fig = figure('Color','w','Position',[100 100 720 520]);
hold on; box on; grid on;
for e = 1:4
    for typ = 1:2
        if typ==1, x = EAC{e}(:); else, x = EDG{e}(:); end
        cnt = histcounts(x, edges, 'Normalization','pdf');
        cnt = smoothdata(cnt,'gaussian',3);
        if typ==1
            plot(ctr, cnt, '--', 'Color', expt(e).color, 'LineWidth',1.2);
        else
            plot(ctr, cnt, '-',  'Color', expt(e).color, 'LineWidth',1.6);
        end
    end
end
xlabel('RMSE (dB)'); ylabel('pdf');
xlim([0 20]); ylim([0 1.15]);
title('Yu et al. (2009) Fig. 3 reproduction');
lgd = [arrayfun(@(e)sprintf('%s AC',expt(e).name),1:4,'UniformOutput',false), ...
       arrayfun(@(e)sprintf('%s DG',expt(e).name),1:4,'UniformOutput',false)];
legend(lgd,'Location','northeast','FontSize',7); legend boxoff;
saveas(fig,fullfile(figdir,'yu2009_fig3.png'));
end
