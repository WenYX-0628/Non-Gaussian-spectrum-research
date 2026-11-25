function verify_yu_sim()
%VERIFY_YU_SIM  对两篇论文的关键仿真结果做数值校验并打印
%
%   校验项：
%     1) 2007：25 km 龙卷谱的边缘速度、谱宽、噪底
%     2) 2009：Expt II/III 在 dv=-20 的 AC 速度偏差（理论 vs 仿真）
%     3) 2009：Expt II/dv=-14 的双高斯拟合参数回收
%     4) 2009：Fig.3 的 RMSE 均值（若已运行 yu2009_fig3）

here = fileparts(mfilename('fullpath'));
addpath(here);
rng(42,'twister');

fprintf('\n========== Yu et al. (2007) 校验 ==========\n');
va = 65; M = 128;
o = struct('va',va,'M',M,'beamwidth',1,'sigma_r',87.5,'gridN',301, ...
    'vel_jitter',0.5,'smooth_bins',1.0,'r0',25000, ...
    'vortices',struct('x',0,'y',25000,'rt',200,'vtm',50,'vrm',0), ...
    'refl',struct('type','doughnut','Wz',60,'r0z',220,'zm',1));
[S,v] = y07_sim_volume(o);
S = S/max(S);
thr = 0.01*max(S);
v1 = v(find(S>thr,1)); v2 = v(find(S>thr,1,'last'));
vm = sum(v.*S)/sum(S);
sv = sqrt(sum((v-vm).^2.*S)/sum(S));
fprintf('25 km 谱: 边缘速度 %.1f / %.1f m/s, 谱宽 %.1f m/s, 均值 %.2f m/s\n',v1,v2,sv,vm);
fprintf('  预期边缘约 ±50 m/s (V_Tm=50)，谱宽约 35–40 m/s\n');
assert(abs(v1+50)<8 && abs(v2-50)<8, '25 km 谱边缘偏离预期');
assert(sv>25 && sv<45, '25 km 谱宽偏离预期');

fprintf('\n========== Yu et al. (2009) 校验 ==========\n');
S1=10; v1=10; s1=2; va=25;
fprintf('Expt II (S2=10 dBm, s2=6 m/s), dv=-20 m/s:\n');
thII = y09_theory_ac(S1,10,v1,-10,s1,6,va);
simII = y09_sim_bias(S1,10,v1,-10,s1,6,va,64,30,200,8,7);
fprintf('  理论: v=%.2f, sv=%.2f, sk=%.2f, kt=%.2f, bias_v=%.2f, bias_sv=%.2f\n', ...
    thII.vr,thII.sv,thII.sk,thII.kt,thII.bias_v,thII.bias_sv);
fprintf('  仿真: bias_v=%.2f±%.2f, bias_sv=%.2f±%.2f\n', ...
    mean(simII.bias_v),std(simII.bias_v),mean(simII.bias_sv),std(simII.bias_sv));
assert(abs(mean(simII.bias_v)-thII.bias_v)<0.6, 'Expt II 速度偏差与理论不符');

fprintf('Expt III (S2=1 dBm, s2=6 m/s), dv=-20 m/s:\n');
thIII = y09_theory_ac(S1,1,v1,-10,s1,6,va);
simIII = y09_sim_bias(S1,1,v1,-10,s1,6,va,64,30,200,8,8);
fprintf('  理论 bias_v=%.2f, 仿真 bias_v=%.2f\n',thIII.bias_v,mean(simIII.bias_v));
assert(abs(mean(simIII.bias_v)-thIII.bias_v)<0.6, 'Expt III 速度偏差与理论不符');

fprintf('双高斯拟合 (Expt II, dv=-14):\n');
[Sp,vv,info] = y09_sim_spectrum(10,10,10,-4,2,6,va,64,30,1,8,11);
fit = y09_fit_dual(Sp,vv,va,struct('nStart',12));
fprintf('  真值 S1=10, S2=10, v1=10, v2=-4, s1=2, s2=6, Pn=%.4f\n',info.Pn);
fprintf('  拟合 S1=%.2f, S2=%.2f, v1=%.2f, v2=%.2f, s1=%.2f, s2=%.2f, Pn=%.4f, RMSE=%.2f dB\n', ...
    fit.S1,fit.S2,fit.v1,fit.v2,fit.s1,fit.s2,fit.Pn,sqrt(fit.mse));
assert(abs(fit.v1-10)<1 && abs(fit.v2+4)<1.5, '拟合平均速度偏差过大');
assert(abs(fit.s1-2)<0.8 && abs(fit.s2-6)<1.5, '拟合谱宽偏差过大');

f3 = fullfile(here,'data','yu2009_fig3.mat');
if exist(f3,'file')
    D = load(f3);
    fprintf('Fig.3 RMSE 均值 (dB):\n');
    for e=1:4
        fprintf('  Expt %s: AC=%.2f, DG=%.2f\n', D.expt(e).name, ...
            mean(D.EAC{e}(:)), mean(D.EDG{e}(:)));
        assert(mean(D.EDG{e}(:))<3, 'DG RMSE 均值异常');
    end
end
fprintf('\n全部校验通过。\n');
end
