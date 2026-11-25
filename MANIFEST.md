# MANIFEST — Yu et al. (2007, 2009) 仿真复现包

- **版本**：v1.0
- **日期**：2026-11-25
- **运行环境**：MATLAB R2023a（建议 R2018b+）
- **必需工具箱**：Optimization Toolbox（`lsqnonlin`）、Signal Processing Toolbox（`hann`）、Statistics and Machine Learning Toolbox（`smoothdata`）
- **可选工具箱**：Parallel Computing Toolbox（仅加速 `yu2009_fig3`）
- **入口**：`run_all.m`；校验：`verify_yu_sim.m`

## 1. 代码清单

### 1.1 公共函数（6）

| 文件 | 作用 |
|---|---|
| `ysim_iq_from_psd.m` | 由目标 PSD 生成复 IQ（Zrnić 1975 谱成形法），含噪声与功率归一化 |
| `ysim_periodogram.m` | 加窗周期图与速度轴（Hann 默认） |
| `ysim_acf_lag01.m` | 滞后 0/1 周期自相关 |
| `ysim_ac_moments.m` | 脉冲对法平均速度/谱宽 |
| `ysim_prib.m` | 径向积分双谱相位 PRIB，含 Chandran–Elgar 插值与 `absfft` 对称性规避 |
| `ysim_spectrum_features.m` | 谱矩、dB 谱标准差 σ_s、PRIB |

### 1.2 Yu et al. (2007)（2 + 5）

| 文件 | 对应内容 |
|---|---|
| `y07_fields.m` | 式(A2)–(A5)：反射率场、Combined Rankine 涡旋、背景风/切变 |
| `y07_sim_volume.m` | 式(A1)：二维体积积分 → 多普勒谱，含波束/距离加权与速度混叠 |
| `yu2007_fig1.m` | Fig.1：距离、背景风、水平切变、反射率宽度 |
| `yu2007_fig2.m` | Fig.2：龙卷+中气旋 PPI 与 5×5 谱 |
| `yu2007_fig4.m` | Fig.4：高斯/矩形图案的 PRIB |
| `yu2007_fig5.m` | Fig.5：Td/Tu/Gd 谱 |
| `yu2007_fig678.m` | Fig.6–8：σ_vn、σ_s、PRIB 距离廓线 |

### 1.3 Yu et al. (2009)（7 + 4）

| 文件 | 对应内容 |
|---|---|
| `y09_model_spectrum.m` | 式(1)：双高斯谱 |
| `y09_dual_moments.m` | 附录 A：双高斯理论矩、偏度、峰度 |
| `y09_theory_ac.m` | 式(2)–(4)：理论自相关与 AC 偏差 |
| `y09_sim_bias.m` | 四组实验的 AC 偏差数值实验 |
| `y09_fit_dual.m` | 式(5)：dB 域多初值非线性拟合 |
| `y09_sim_spectrum.m` | 含噪周期图 + AC 矩 |
| `y09_rmse_one.m` | 单次实现的 AC/DG 重构 RMSE |
| `yu2009_fig1.m` | Fig.1：理论偏差与 200 次实现验证 |
| `yu2009_fig2.m` | Fig.2：观测/AC 重构/双高斯拟合谱 |
| `yu2009_fig3.m` | Fig.3：RMSE 概率密度 |
| `yu2009_fig8.m` | Fig.8：垂直切变双峰/宽裙谱 |

### 1.4 驱动与校验（2）

| 文件 | 作用 |
|---|---|
| `run_all.m` | 顺序运行全部图件；`run_all(true)` 为快速模式 |
| `verify_yu_sim.m` | 关键数值校验并断言 |

## 2. 图件映射

| 输出文件 | 论文图 | 内容 |
|---|---|---|
| `figs/yu2007_fig1.png` | 2007 Fig.1 | 距离/风/切变/反射率对谱形的影响 |
| `figs/yu2007_fig2a.png` | 2007 Fig.2 上 | 反射率与径向速度 PPI |
| `figs/yu2007_fig2b.png` | 2007 Fig.2 下 | 黑框内 5×5 含噪谱 |
| `figs/yu2007_fig4.png` | 2007 Fig.4 | PRIB vs 谱宽 |
| `figs/yu2007_fig5.png` | 2007 Fig.5 | 9 个体积的 Td/Tu/Gd |
| `figs/yu2007_fig6.png` | 2007 Fig.6 | 归一化谱宽距离廓线 |
| `figs/yu2007_fig7.png` | 2007 Fig.7 | 谱平坦度距离廓线 |
| `figs/yu2007_fig8.png` | 2007 Fig.8 | PRIB 距离廓线 |
| `figs/yu2009_fig1.png` | 2009 Fig.1 | 双高斯 AC 偏差 |
| `figs/yu2009_fig2.png` | 2009 Fig.2 | 谱与两种重构 |
| `figs/yu2009_fig3.png` | 2009 Fig.3 | RMSE 概率密度 |
| `figs/yu2009_fig8.png` | 2009 Fig.8 | 垂直切变谱 |

## 3. 数据文件

| 文件 | 内容 |
|---|---|
| `data/yu2007_fig1.mat` | 12 条谱曲线 |
| `data/yu2007_fig2.mat` | PPI 场、5×5 谱、AC 速度 |
| `data/yu2007_fig4.mat` | PRIB–σ_v 曲线 |
| `data/yu2007_fig5.mat` | 9×3 条谱 |
| `data/yu2007_fig678.mat` | 9 面板 × 20 r_n 的统计量 |
| `data/yu2009_fig1.mat` | 四组实验理论/仿真偏差 |
| `data/yu2009_fig2.mat` | 12 个示例的观测/拟合 |
| `data/yu2009_fig3.mat` | AC/DG RMSE 样本 |
| `data/yu2009_fig8.mat` | 四组垂直切变谱 |

`data/` 可删除；`run_all` 会重新生成。

## 4. 运行

```matlab
cd Non-Gaussian-spectrum-research
run_all          % 全部图件，约 2 分钟
run_all(true)    % 快速模式
verify_yu_sim    % 数值校验
```

## 5. 校验摘要

| 项目 | 结果 |
|---|---|
| 2007 25 km 谱边缘 | −52.3 / +51.3 m/s |
| 2007 25 km 谱宽 | 35.5 m/s |
| 2009 双高斯拟合 RMSE 均值 | 1.40 / 1.71 / 1.55 / 1.40 dB |
| 2009 Expt II 最大速度偏差 | 约 2.7–3.0 m/s |
| 2009 Expt II/Δv=−14 参数回收 | v₁≈9.7、v₂≈−4.5 m/s；σ₁≈2.2、σ₂≈5.8 m/s |

详细讨论见 `README.md` 第 3–5 节和 `docs/verification_report.md`。

## 6. 完整性

发布 zip 包内附 `SHA256SUMS.txt`，覆盖包内除自身外的所有文件；仓库中不重复提交校验和文件。校验方法：

```powershell
Get-FileHash -Algorithm SHA256 <file>
```

