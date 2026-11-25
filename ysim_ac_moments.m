function out = ysim_ac_moments(iq, va, noisePower)
%YSIM_AC_MOMENTS  自协方差(脉冲对)法估计平均速度与谱宽
%
%   out = ysim_ac_moments(iq, va, noisePower)
%
%   iq         : [M x N] 复 IQ
%   va         : 最大不模糊速度 [m/s]
%   noisePower : 每样本噪声功率（可选，默认 0）。若给出，则先从 R0 中扣除。
%
%   公式（Doviak & Zrnic 1993；Yu et al. 2009 式(3)(4)）：
%       v_hat  = -(va/pi) * angle(R1)
%       sv_hat = (sqrt(2)*va/pi) * sqrt( ln(R0/|R1|) )
%
%   对高斯谱无偏；对双高斯等非高斯谱存在偏差，正是 Yu et al. (2009)
%   要量化的对象。

if nargin < 3 || isempty(noisePower), noisePower = 0; end
[R0, R1] = ysim_acf_lag01(iq);
R0 = max(R0 - noisePower, eps);
A  = max(abs(R1), eps);
out.R0 = R0;
out.R1 = R1;
out.v  = -(va/pi) * angle(R1);
out.sv = (sqrt(2)*va/pi) * sqrt(max(log(R0./A), 0));
out.snr = 10*log10(max(R0 - noisePower,eps)./max(noisePower,eps));
end
