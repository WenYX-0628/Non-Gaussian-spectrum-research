function feat = ysim_spectrum_features(S, v, pribOpt)
%YSIM_SPECTRUM_FEATURES  由多普勒谱计算谱宽、谱平坦度和 PRIB
%
%   feat = ysim_spectrum_features(S, v, pribOpt)
%
%   S : 功率谱（线性功率，任意单位）
%   v : 与 S 对应的径向速度轴 [m/s]，长度 M
%   pribOpt : 可选结构体
%       .a      PRIB 斜率，默认 0.9
%       .nfft   PRIB FFT 点数，默认 256
%       .method 'absfft'（默认）或 'direct'
%       .floor_db 计算 dB 谱时的功率下限（相对 S 峰值），默认 -80 dB
%
%   输出 feat 字段
%       .P       总功率 sum(S)
%       .vm      一阶谱矩（平均径向速度）
%       .sv      二阶中心矩平方根（谱宽）
%       .x       dB 谱 x=10log10(S/Smax)
%       .sigma_s 谱平坦度：dB 谱的标准差（论文 Fig.3/Fig.7 的 sigma_s）
%       .PRIB    PRIB 相位 [deg]
%       .snr_db  按峰值与中位数底噪估计的 SNR（仅参考）
%
%   注意：论文把谱平坦度 sigma_s 定义为 x(k) 的“统计方差”，但 Fig.3/Fig.7
%   的数值范围（约 5~15，且白噪声谱约 5.6）表明实际使用的是标准差；
%   本函数返回标准差，并在 README 中说明。

if nargin < 3 || isempty(pribOpt), pribOpt = struct(); end
if ~isfield(pribOpt,'a'),        pribOpt.a = 0.9;          end
if ~isfield(pribOpt,'nfft'),     pribOpt.nfft = 256;       end
if ~isfield(pribOpt,'method'),   pribOpt.method = 'absfft'; end
if ~isfield(pribOpt,'floor_db'), pribOpt.floor_db = -80;   end

S = S(:); v = v(:);
P = sum(S);
if P <= 0 || ~isfinite(P)
    feat = struct('P',P,'vm',NaN,'sv',NaN,'x',[],'sigma_s',NaN, ...
                  'PRIB',NaN,'snr_db',NaN);
    return;
end

vm = sum(v.*S)/P;
sv = sqrt(max(sum((v-vm).^2.*S)/P, 0));

Smax = max(S);
x = 10*log10(max(S/Smax, 10^(pribOpt.floor_db/10)));

feat.P       = P;
feat.vm      = vm;
feat.sv      = sv;
feat.x       = x;
feat.sigma_s = std(x, 1);          % 总体标准差（M 个 bin）
feat.PRIB    = ysim_prib(x, pribOpt.a, pribOpt.nfft, pribOpt.method);

% 粗 SNR：峰值 bin 与最低 20% bin 中位数之比
q = sort(S);
noise = median(q(1:max(1,round(0.2*numel(q)))));
feat.snr_db = 10*log10(max(Smax,eps)/max(noise,eps));
end
