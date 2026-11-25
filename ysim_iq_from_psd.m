function [iq, Ps, Pn] = ysim_iq_from_psd(S, M, snr_db, nreal, seed)
%YSIM_IQ_FROM_PSD  由目标多普勒功率谱生成复 IQ 时间序列 (Zrnic 1975 谱成形法)
%
%   [iq, Ps, Pn] = ysim_iq_from_psd(S, M, snr_db, nreal, seed)
%
%   S      : [M x 1] 目标功率谱（线性功率，一个 Nyquist 周期，速度 bin
%            -va:dv:va-dv）。总信号功率 Ps=sum(S)。
%   M      : 每实现的脉冲数；默认 numel(S)
%   snr_db : 加性白噪声信噪比 10log10(Ps/Pn)；Inf 表示不加噪声
%   nreal  : 实现数（列数），默认 1
%   seed   : 随机种子（可选）
%
%   iq     : [M x nreal] 复 IQ
%
%   原理：频域独立复高斯系数 X_k ~ CN(0,S_k)，iq = ifft(X)*sqrt(M)，
%   则 E|iq_n|^2 = sum(S)，周期图期望为 S。与 Zrnic (1975) 的谱成形法一致。

if nargin < 2 || isempty(M),      M = numel(S);     end
if nargin < 3 || isempty(snr_db), snr_db = Inf;     end
if nargin < 4 || isempty(nreal),  nreal = 1;        end
if nargin >= 5 && ~isempty(seed)
    rng(seed, 'twister');
end

S = S(:);
if numel(S) ~= M
    error('ysim_iq_from_psd: S must have M elements.');
end
S = max(real(S), 0);

% S 按速度轴 -va:dv:va-dv 排列；FFT 的频率顺序是 0,+...,-...，
% 因此先 ifftshift 再作为频域系数，否则会产生一个固定速度偏移。
Sf = ifftshift(S);
X = sqrt(Sf/2) .* (randn(M, nreal) + 1i*randn(M, nreal));
% 取共轭使正速度对应负多普勒相位（与脉冲对公式 v=-(va/pi)angle(R1) 一致）
% ifft 自带 1/M，乘 M 后 E|iq|^2 = sum(S)（总功率守恒）。
iq = conj(ifft(X, [], 1)) * M;

Ps = sum(S);
if isfinite(snr_db)
    Pn = Ps / 10^(snr_db/10);
    iq = iq + sqrt(Pn/2) .* (randn(M, nreal) + 1i*randn(M, nreal));
else
    Pn = 0;
end
end
