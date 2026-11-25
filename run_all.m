function run_all(quick)
%RUN_ALL  运行 Yu et al. (2007, 2009) 全部仿真图件
%
%   run_all        % 完整运行（Fig.3 用 30 次实现，Fig.6-8 用 50 次）
%   run_all(true)  % 快速测试（减少实现次数）

if nargin < 1 || isempty(quick), quick = false; end
here = fileparts(mfilename('fullpath'));
addpath(here);
if ~exist(fullfile(here,'figs'),'dir'), mkdir(fullfile(here,'figs')); end
if ~exist(fullfile(here,'data'),'dir'), mkdir(fullfile(here,'data')); end

t0 = tic;
fprintf('=== Yu et al. (2007) ===\n');
fprintf('Fig.1 ...\n');    yu2007_fig1;
fprintf('Fig.2 ...\n');    yu2007_fig2;
fprintf('Fig.4 ...\n');    yu2007_fig4;
fprintf('Fig.5 ...\n');    yu2007_fig5;
if quick
    fprintf('Fig.6-8 ... (quick: 10 realizations)\n'); yu2007_fig678(10);
else
    fprintf('Fig.6-8 ... (50 realizations)\n');         yu2007_fig678(50);
end

fprintf('=== Yu et al. (2009) ===\n');
fprintf('Fig.1 ...\n');    yu2009_fig1(true);
fprintf('Fig.2 ...\n');    yu2009_fig2;
if quick
    fprintf('Fig.3 ... (quick: 8 realizations)\n');  yu2009_fig3(8);
else
    fprintf('Fig.3 ... (30 realizations)\n');         yu2009_fig3(30);
end
fprintf('Fig.8 ...\n');    yu2009_fig8;

fprintf('All done in %.1f s. Figures: %s\n', toc(t0), fullfile(here,'figs'));
end
