function out = y09_theory_ac(S1,S2,v1,v2,s1,s2,va)
%Y09_THEORY_AC 双高斯谱的自协方差法理论估计与偏差
%
%   out = y09_theory_ac(S1,S2,v1,v2,s1,s2,va)
%
%   实现 Yu et al. (2009) 式(2)-(4)：
%       R(1) = sum_i S_i exp[-0.5(pi s_i/va)^2] exp(-j pi v_i/va)
%       v_hat = -(va/pi) angle(R(1))
%       s_hat = (sqrt(2) va/pi) sqrt(ln(S/|R(1)|))
%   并与附录A的理论矩比较，给出偏差。

if nargin < 7 || isempty(va), va = 25; end
S = S1 + S2;
A1 = S1*exp(-0.5*(pi*s1/va)^2);
A2 = S2*exp(-0.5*(pi*s2/va)^2);
phi1 = -pi*v1/va;
phi2 = -pi*v2/va;
R1 = A1*exp(1i*phi1) + A2*exp(1i*phi2);

vhat = -(va/pi)*angle(R1);
shat = (sqrt(2)*va/pi)*sqrt(max(log(S/max(abs(R1),eps)),0));

th = y09_dual_moments(S1,S2,v1,v2,s1,s2);
out = th;
out.R1 = R1; out.vhat = vhat; out.shat = shat;
out.bias_v  = vhat - th.vr;
out.bias_sv = shat - th.sv;
end
