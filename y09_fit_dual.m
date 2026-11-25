function fit = y09_fit_dual(Sobs, v, va, opt)
%Y09_FIT_DUAL 双高斯谱非线性拟合（dB 域最小二乘，多初值）
%
%   fit = y09_fit_dual(Sobs, v, va, opt)
%
%   实现 Yu et al. (2009) 式(5)：在 dB 域最小化观测谱与
%       S_model(v) = S1 N(v;v1,s1) + S2 N(v;v2,s2) + Pn/(2va)
%   的均方误差；N 为归一化高斯 pdf。未知量为 6 个谱矩 + 噪声功率。
%
%   输入
%     Sobs : 观测功率谱密度（功率/(m/s)），长度 M
%     v    : 速度轴 [m/s]
%     va   : 最大不模糊速度 [m/s]
%     opt  : 可选
%        .nStart  初值个数上限（默认 24）
%        .verbose 是否打印（默认 false）
%
%   输出 fit
%     .S1,.S2 (mW), .v1,.v2 (m/s), .s1,.s2 (m/s), .Pn (mW)
%     .mse (dB^2), .yfit (dB), .success

if nargin < 4 || isempty(opt), opt = struct(); end
if ~isfield(opt,'nStart'), opt.nStart = 24; end
if ~isfield(opt,'verbose'), opt.verbose = false; end

Sobs = Sobs(:); v = v(:);
dv = mean(diff(v));
y = 10*log10(max(Sobs, eps));
M = numel(Sobs);
Ptot = sum(Sobs)*abs(dv);

% 噪声初值：最低 30% 谱值的中位数（换算成总噪声功率）
q = sort(Sobs);
nlow = max(3, round(0.3*M));
Pn0 = max(median(q(1:nlow))*2*va, 1e-9);
Sig0 = max(Ptot - Pn0, 1e-9);

[~, ipk] = max(Sobs);
vpk = v(ipk);

% 初值集合
v2c = unique([-va, -0.5*va, 0, 0.5*va, va, vpk-8, vpk+8, vpk-15, vpk+15]);
v2c = min(max(v2c, -va+0.5), va-0.5);
spairs = [1.5 1.5; 1.5 4; 4 1.5; 4 4; 2 6; 6 2];

starts = zeros(0,7);
for i = 1:numel(v2c)
    for j = 1:size(spairs,1)
        starts(end+1,:) = [log10(Sig0/2), log10(Sig0/2), vpk, v2c(i), ...
                           log10(spairs(j,1)), log10(spairs(j,2)), log10(Pn0)]; %#ok<AGROW>
    end
end
if size(starts,1) > opt.nStart
    idx = round(linspace(1,size(starts,1),opt.nStart));
    starts = starts(idx,:);
end

lb = [-60 -60 -va+0.2 -va+0.2 log10(0.3) log10(0.3) -60];
ub = [ 30  30  va-0.2  va-0.2 log10(15)  log10(15)   30];

opts = optimoptions('lsqnonlin','Display','off', ...
    'MaxIterations',100,'MaxFunctionEvaluations',1200, ...
    'FunctionTolerance',1e-9,'StepTolerance',1e-10);

best.mse = inf; best.p = []; best.exitflag = -1;
for s = 1:size(starts,1)
    p0 = min(max(starts(s,:), lb+1e-6), ub-1e-6);
    [p,resnorm,~,ef] = lsqnonlin(@(pp)resid(pp,v,y,va), p0, lb, ub, opts);
    mse = resnorm/M;
    if mse < best.mse
        best.mse = mse; best.p = p; best.exitflag = ef;
    end
end

p = best.p;
S1 = 10^p(1); S2 = 10^p(2);
v1 = p(3); v2 = p(4);
s1 = 10^p(5); s2 = 10^p(6);
Pn = 10^p(7);
Sm = y09_model_spectrum(S1,S2,v1,v2,s1,s2,v) + Pn/(2*va);
yfit = 10*log10(max(Sm,eps));

fit = struct('S1',S1,'S2',S2,'v1',v1,'v2',v2,'s1',s1,'s2',s2, ...
             'Pn',Pn,'mse',best.mse,'yfit',yfit,'exitflag',best.exitflag, ...
             'Ptot',S1+S2);
% 按平均速度排序，便于解释（S1 为平均速度较大的分量）
if fit.v1 < fit.v2
    [fit.v1,fit.v2] = deal(fit.v2,fit.v1);
    [fit.s1,fit.s2] = deal(fit.s2,fit.s1);
    [fit.S1,fit.S2] = deal(fit.S2,fit.S1);
end
end

function r = resid(p, v, y, va)
S1 = 10^p(1); S2 = 10^p(2);
s1 = 10^p(5); s2 = 10^p(6);
Pn = 10^p(7);
Sm = y09_model_spectrum(S1,S2,p(3),p(4),s1,s2,v) + Pn/(2*va);
r = y - 10*log10(max(Sm,eps));
end
