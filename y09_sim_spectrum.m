function [Spsd, v, info] = y09_sim_spectrum(S1,S2,v1,v2,s1,s2,va,M,snr_db,nreal,ntrial,seed)
%Y09_SIM_SPECTRUM 生成双高斯谱的含噪周期图（可多实现、多 trial 平均）
%
%   [Spsd,v,info] = y09_sim_spectrum(...)
%
%   每次实现生成 ntrial 条独立 IQ，分别做 Hann 周期图后平均；返回
%   功率谱密度 Spsd（功率/(m/s)），以及真值、噪声功率等。

if nargin < 8  || isempty(M),      M = 64;  end
if nargin < 9  || isempty(snr_db), snr_db = 30; end
if nargin < 10 || isempty(nreal),  nreal = 1; end
if nargin < 11 || isempty(ntrial), ntrial = 8; end
if nargin < 12, seed = []; end
if ~isempty(seed), rng(seed,'twister'); end

v = (-va:2*va/M:va-2*va/M).';
S = y09_model_spectrum(S1,S2,v1,v2,s1,s2,v);
dv = 2*va/M;

% ysim_iq_from_psd 的输入是“每 bin 功率”，而模型给出的是功率谱密度
% （功率/(m/s)），故乘以 bin 宽度。
iq = ysim_iq_from_psd(S*dv, M, snr_db, nreal*ntrial);
Sbin = zeros(M, nreal);
for it = 1:ntrial
    idx = (it-1)*nreal + (1:nreal);
    [Stmp, vt] = ysim_periodogram(iq(:,idx), va, M, 'hann');
    Sbin = Sbin + Stmp;
end
Sbin = Sbin/ntrial;
v = vt;
% 周期图 Sbin 的 sum = M * 平均功率；功率谱密度（积分=平均功率）需除以
% M*dv = 2*va。
Spsd = Sbin/(M*dv);

% 同一批 IQ 的 8-trial 平均自相关，用于自协方差法对比
[R0all, R1all] = ysim_acf_lag01(iq);
R0all = reshape(R0all, ntrial, nreal);
R1all = reshape(R1all, ntrial, nreal);
R0 = mean(R0all,1); R1 = mean(R1all,1);
vac = -(va/pi)*angle(R1);
sac = (sqrt(2)*va/pi)*sqrt(max(log(max(R0,eps)./max(abs(R1),eps)),0));

th = y09_dual_moments(S1,S2,v1,v2,s1,s2);
info = struct('Sbin',Sbin,'Spsd',Spsd,'v',v,'dv',dv, ...
              'vtrue',th.vr,'svtrue',th.sv,'sk',th.sk,'kt',th.kt, ...
              'Ps',sum(S)*dv,'Pn',sum(S)*dv/10^(snr_db/10), ...
              'va',va,'M',M,'snr_db',snr_db, ...
              'vac',vac,'sac',sac,'R0',R0,'R1',R1);
if nreal == 1
    info.Sbin = Sbin;
end
end
