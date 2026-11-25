function out = y09_sim_bias(S1,S2,v1,v2,s1,s2,va,M,snr_db,nreal,ntrial,seed)
%Y09_SIM_BIAS 双高斯谱的 IQ 数值实验：自协方差法偏差
%
%   out = y09_sim_bias(S1,S2,v1,v2,s1,s2,va,M,snr_db,nreal,ntrial,seed)
%
%   对每次实现生成 ntrial 条独立 IQ 记录，先平均其自相关 R0/R1，
%   再用式(3)(4)估计平均速度与谱宽（与 Yu et al. 2009 第3b节一致）。
%
%   输出
%     out.vhat, out.shat : [1 x nreal] 估计值
%     out.R0, out.R1     : 平均后的自相关
%     out.vtrue, out.svtrue, out.sk, out.kt : 理论值

if nargin < 9  || isempty(M),      M = 64;    end
if nargin < 10 || isempty(snr_db), snr_db = 30; end
if nargin < 11 || isempty(nreal),  nreal = 200; end
if nargin < 12 || isempty(ntrial), ntrial = 8;  end
if nargin < 13, seed = []; end
if ~isempty(seed), rng(seed,'twister'); end

v = (-va:2*va/M:va-2*va/M).';
S = y09_model_spectrum(S1,S2,v1,v2,s1,s2,v);
dv = 2*va/M;

iq = ysim_iq_from_psd(S*dv, M, snr_db, nreal*ntrial);
[R0all, R1all] = ysim_acf_lag01(iq);
R0all = reshape(R0all, ntrial, nreal);
R1all = reshape(R1all, ntrial, nreal);
R0 = mean(R0all,1);
R1 = mean(R1all,1);

vhat = -(va/pi)*angle(R1);
shat = (sqrt(2)*va/pi)*sqrt(max(log(max(R0,eps)./max(abs(R1),eps)),0));

th = y09_dual_moments(S1,S2,v1,v2,s1,s2);
out.vhat = vhat; out.shat = shat;
out.R0 = R0; out.R1 = R1;
out.vtrue = th.vr; out.svtrue = th.sv; out.sk = th.sk; out.kt = th.kt;
out.bias_v = vhat - th.vr;
out.bias_sv = shat - th.sv;
out.Pn = sum(S)*dv/10^(snr_db/10);
end
