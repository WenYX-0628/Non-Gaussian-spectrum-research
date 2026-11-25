function [Pdeg, I] = ysim_prib(x, a, nfft, method)
%YSIM_PRIB  径向积分双谱相位 (Phase of Radially Integrated Bispectrum, PRIB)
%
%   P = ysim_prib(x, a, nfft, method)
%
%   对应 Yu et al. (2007, JTECH) 式(3):
%       I(a) = ∫_0^{1/(1+a)} B(f1, a f1) df1
%       P    = arg(I)
%   其中 B(f1,f2) = X(f1) X(f2) X*(f1+f2), X 为输入一维序列 x(k) 的 DFT。
%   积分路径 f2 = a*f1，a=0.9（论文取值）。双谱在离散频点上按
%   Chandran & Elgar (1993) 的方式线性插值（插值双谱而非插值傅里叶系数）。
%
%   输入
%     x      : 一维序列（论文中为 dB 多普勒谱 x(k)=10log10 S(vk)）
%     a      : 积分路径斜率，默认 0.9
%     nfft   : FFT 点数（论文补零到 256），默认 256
%     method : 'direct'  : 直接对 x 做双谱（偶对称序列相位退化为 0/180°）
%              'absfft'  : 论文第3b节给出的对称性规避方法——先取 |X(f)|，
%                          再对其前一半作双谱（默认）
%
%   输出
%     Pdeg : 相位 [deg]，范围 (-180,180]
%     I    : 复积分值
%
%   备注：论文没有给出可直接复现的代码；对完全对称的无噪图案，直接法
%   的双谱为实数，P 退化为 0/180°，因此这里默认采用论文明确描述的
%   |X(f)| 替代序列方法。实际含噪谱两种方法给出的统计分布不同。

if nargin < 2 || isempty(a),      a = 0.9;        end
if nargin < 3 || isempty(nfft),   nfft = 256;     end
if nargin < 4 || isempty(method), method = 'absfft'; end

x = x(:);
x = x(isfinite(x));
if isempty(x)
    Pdeg = NaN; I = NaN; return;
end

switch lower(method)
    case 'direct'
        y = zeros(nfft, 1);
        y(1:min(numel(x), nfft)) = x(1:min(numel(x), nfft));
        X = fft(y);
    case 'absfft'
        X0 = fft(x, nfft);
        A  = abs(X0(1:nfft/2));
        A  = A - mean(A);
        y  = zeros(nfft, 1);
        y(1:numel(A)) = A;
        X  = fft(y);
    otherwise
        error('ysim_prib: unknown method "%s".', method);
end

% 去直流与 Nyquist（Chandran & Elgar 代码中的处理，避免插值引入直流依赖）
X(1) = 0;
X(nfft/2 + 1) = 0;

% 积分上限：f1 + f2 <= 0.5 (归一化到采样频率)，f2=a f1
kmax = floor(0.5*nfft/(1 + a)) - 1;
if kmax < 1
    Pdeg = NaN; I = NaN; return;
end
k   = (1:kmax).';
f2  = a*k;
lo  = floor(f2);
hi  = lo + 1;
w1  = f2 - lo;          % 对 ceil 的权重
w2  = 1 - w1;           % 对 floor 的权重

Bfloor = X(k+1) .* X(lo+1) .* conj(X(k+lo+1));
Bceil  = X(k+1) .* X(hi+1) .* conj(X(k+hi+1));
I = sum(w1.*Bceil + w2.*Bfloor);

Pdeg = atan2d(imag(I), real(I));
end
