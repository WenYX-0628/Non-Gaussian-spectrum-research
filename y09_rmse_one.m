function [eac, edg, fit] = y09_rmse_one(S1,S2,v1,v2,s1,s2,va,M,snr_db,seed,nStart)
%Y09_RMSE_ONE 单次实现的 AC 重构与双高斯拟合 RMSE
%
%   返回
%     eac : 自协方差法（高斯重构）RMSE [dB]
%     edg : 双高斯拟合 RMSE [dB]
%     fit : 拟合结构体

if nargin < 11 || isempty(nStart), nStart = 8; end
[Sp,v,info] = y09_sim_spectrum(S1,S2,v1,v2,s1,s2,va,M,snr_db,1,8,seed);
yobs = 10*log10(max(Sp,eps));

% 自协方差法重构
Pnpsd = info.Pn/(2*va);
Pac = max(sum(Sp)*info.dv - info.Pn, 1e-9);
Sac = Pac*exp(-(v-info.vac).^2/(2*max(info.sac,0.2)^2))/(sqrt(2*pi)*max(info.sac,0.2)) + Pnpsd;
eac = sqrt(mean((yobs - 10*log10(max(Sac,eps))).^2));

% 双高斯拟合
fit = y09_fit_dual(Sp,v,va,struct('nStart',nStart));
edg = sqrt(fit.mse);
end
