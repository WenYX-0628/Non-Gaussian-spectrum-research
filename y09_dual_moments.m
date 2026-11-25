function out = y09_dual_moments(S1,S2,v1,v2,s1,s2)
%Y09_DUAL_MOMENTS 双高斯谱的理论谱矩 (Yu et al. 2009, 附录A 式A1-A4)
%
%   out = y09_dual_moments(S1,S2,v1,v2,s1,s2)
%
%   Si : 第 i 个高斯分量的信号功率（线性，任意单位；dBm 需先换算）
%   vi : 平均径向速度 [m/s]
%   si : 谱宽（标准差）[m/s]
%
%   输出：Ptot, vr, sv, mu3, mu4, sk, kt
%   其中 sk=mu3/sv^3（高斯=0），kt=mu4/sv^4（高斯=3）。

S  = S1 + S2;
P1 = S1/S;
P2 = S2/S;
d  = v2 - v1;

vr = P1*v1 + P2*v2;
sv2 = P1*s1^2 + P2*s2^2 + P1*P2*d^2;
sv  = sqrt(max(sv2,0));

mu3 = 3*P1*P2*d*(s1^2 - s2^2) + P1*P2*d^3*(P2 - P1);
mu4 = 3*(P1*s1^4 + P2*s2^4) ...
    + 6*P1*P2*d^2*(P2*s1^2 + P1*s2^2) ...
    + P1*P2*d^4*(P1^3 + P2^3);

out = struct();
out.Ptot = S; out.P1 = P1; out.P2 = P2;
out.vr = vr; out.sv = sv;
out.mu3 = mu3; out.mu4 = mu4;
out.sk = mu3/max(sv^3,eps);
out.kt = mu4/max(sv^4,eps);
end
