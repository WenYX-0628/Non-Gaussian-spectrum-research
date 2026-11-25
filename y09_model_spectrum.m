function S = y09_model_spectrum(S1,S2,v1,v2,s1,s2,v)
%Y09_MODEL_SPECTRUM 双高斯多普勒谱模型 (Yu et al. 2009, 式(1))
%
%   S = y09_model_spectrum(S1,S2,v1,v2,s1,s2,v)
%
%   Si 为线性功率（mW），vi [m/s]，si [m/s]。
%   Si 的单位与 S 相同（功率/(m/s)）。

v = v(:);
g1 = exp(-(v-v1).^2/(2*s1^2)) / (sqrt(2*pi)*s1);
g2 = exp(-(v-v2).^2/(2*s2^2)) / (sqrt(2*pi)*s2);
S = S1*g1 + S2*g2;
end
