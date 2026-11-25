function [R0, R1] = ysim_acf_lag01(iq)
%YSIM_ACF_LAG01  复 IQ 在滞后 0、1 的自相关（脉冲对法）
%
%   [R0, R1] = ysim_acf_lag01(iq)
%
%   R0 = mean(|x|^2)
%   R1 = mean(x(n+1) * conj(x(n)))
%
%   对 [M x N] 输入按列计算，R0/R1 为 [1 x N]。

iq = iq(:,:);
R0 = mean(abs(iq).^2, 1);
% 采用周期（循环）滞后 1：对有限样本可避免丢掉最后一个脉冲对造成的
% 系统性谱宽低估；与理论自相关 R(1) 的定义一致。
R1 = mean(iq([2:end 1],:) .* conj(iq), 1);
end
