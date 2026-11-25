function [S, v] = ysim_periodogram(iq, va, nfft, winName)
%YSIM_PERIODOGRAM  加窗周期图（速度轴已 fftshift）
%
%   [S, v] = ysim_periodogram(iq, va, nfft, winName)
%
%   iq      : [M x 1] 复 IQ
%   va      : 最大不模糊速度 [m/s]
%   nfft    : FFT 点数，默认 M（补零不改变总功率）
%   winName : 'hann'（默认，von Hann 窗）、'rect'、'hamming'
%
%   S : [nfft x 1] 每速度 bin 的功率，sum(S) ≈ sum(|iq|^2)
%   v : [nfft x 1] 速度轴 [m/s]，-va:dv:va-dv，dv=2va/nfft

iq = iq(:);
M = numel(iq);
if nargin < 2 || isempty(va),      va = 0.5;          end
if nargin < 3 || isempty(nfft),    nfft = M;          end
if nargin < 4 || isempty(winName), winName = 'hann'; end

switch lower(winName)
    case 'hann'
        w = 0.5 - 0.5*cos(2*pi*(0:M-1)'/max(M-1,1));
    case 'hamming'
        w = 0.54 - 0.46*cos(2*pi*(0:M-1)'/max(M-1,1));
    case 'rect'
        w = ones(M,1);
    otherwise
        error('ysim_periodogram: unknown window "%s".', winName);
end

X = fft(iq.*w, nfft);
% 归一化使白噪声情形 sum(S)=sum(|iq|^2)
S = abs(X).^2 / max(sum(w.^2), eps) * (M/nfft);
S = fftshift(S);
faxis = (-nfft/2 : nfft/2-1).';
% 本项目的 IQ 生成约定：正速度对应负多普勒频率（与脉冲对公式
% v=-(va/pi)angle(R1) 一致）。因此周期图频率轴取负号后再排序，
% 使返回的 v 与 y09_model_spectrum 中定义的正速度方向一致。
v = -faxis * (2*va/nfft);
[v, idx] = sort(v);
S = S(idx);
end
