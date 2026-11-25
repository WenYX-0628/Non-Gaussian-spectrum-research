# Non-Gaussian-spectrum-research

Weather radar non-Gaussian spectrum research.

本仓库收录两篇关于龙卷/非高斯多普勒谱的 MATLAB 仿真复现：

- **Yu et al. (2007)**：高阶谱（双谱 / PRIB）与谱平坦度表征龙卷谱
- **Yu et al. (2009)**：双高斯谱模型、自协方差偏差与垂直切变解释

> 仓库地址：<https://github.com/WenYX-0628/Non-Gaussian-spectrum-research>  
> 当前版本：v1.0（2026-11-25）

---
## Yu et al. (2007, 2009) 龙卷/非高斯多普勒谱仿真复现

**版本**：v1.0（2026-11-25）  
**仓库根目录**：本仓库根目录（v1.0）

本仓库根目录用 MATLAB 复现两篇论文中的数值仿真：

1. **Yu, Wang, Shapiro, Yeary, Zrnić & Doviak (2007)**, *Characterization of Tornado Spectral Signatures Using Higher-Order Spectra*, J. Atmos. Oceanic Technol., 24, 1997–2013.
2. **Yu, Reinoso-Rondinel & Palmer (2009)**, *Investigation of Non-Gaussian Doppler Spectra Observed by Weather Radar in a Tornadic Supercell*, J. Atmos. Oceanic Technol., 26, 444–461.

代码不依赖论文原始数据；所有谱、IQ、统计量均由本仓库函数生成。运行环境为 MATLAB R2023a（R2018b 及以上应可），使用了 Optimization Toolbox（`lsqnonlin`）、Signal/Statistics Toolbox（`hann`、`smoothdata`）和 Parallel Computing Toolbox（可选，仅用于加速 Fig.3）。

---

## 1. 快速开始

```matlab
cd('...\非高斯谱\yu_sim')
run_all            % 生成全部图件与 data/*.mat，约 2 分钟（4 核并行时）
```

也可单独运行：

```matlab
yu2007_fig1        % 图1：距离/背景风/切变/反射率结构
yu2007_fig2        % 图2：龙卷+中气旋的 PPI 与 5x5 谱
yu2007_fig4        % 图4：高斯/矩形图案的 PRIB（见第6节说明）
yu2007_fig5        % 图5：9 个体积的 Td/Tu/Gd 谱
yu2007_fig678(50)  % 图6–8：谱宽/谱平坦度/PRIB 的距离廓线
yu2009_fig1(true)  % 图1：双高斯谱 AC 偏差理论与 200 次实现验证
yu2009_fig2        % 图2：观测谱、AC 重构、双高斯拟合
yu2009_fig3(30)    % 图3：两种方法的 RMSE 概率密度
yu2009_fig8        % 图8：垂直切变产生双峰/宽裙谱
```

图件输出到 `figs/`，中间数据输出到 `data/`。`yu2007_fig2` 生成 `yu2007_fig2a.png`（PPI）与 `yu2007_fig2b.png`（5×5 谱）两张，对应论文 Fig.2 的上下两部分。

---

## 2. 文件与论文公式对应

### 公共函数

| 文件 | 作用 |
|---|---|
| `ysim_iq_from_psd.m` | Zrnić (1975) 频域谱成形法：由目标 PSD 生成复 IQ；含白噪声、功率归一化 |
| `ysim_periodogram.m` | 加窗周期图（Hann 默认），速度轴与 IQ 符号约定一致 |
| `ysim_acf_lag01.m` | 滞后 0/1 的周期自相关（脉冲对法） |
| `ysim_ac_moments.m` | 式(3)(4)：$\hat v=-(v_a/\pi)\angle R_1$，$\hat\sigma_v=(\sqrt2 v_a/\pi)\sqrt{\ln(R_0/|R_1|)}$ |
| `ysim_prib.m` | 径向积分双谱相位 PRIB：$I(a)=\int_0^{1/(1+a)}B(f_1,af_1)df_1$，$P=\arg I$；含 Chandran–Elgar 插值与 `|X(f)|` 对称性规避 |
| `ysim_spectrum_features.m` | 由谱计算总功率、谱矩、dB 谱标准差 $\sigma_s$、PRIB |

### Yu et al. (2007)

| 文件 | 论文内容 |
|---|---|
| `y07_fields.m` | 式(A2)–(A5)：环状/高斯/均匀反射率，Combined Rankine 涡旋（$\gamma=\pm1$），背景风、水平切变、多涡旋叠加 |
| `y07_sim_volume.m` | 式(A1)：$S(r_0,v_k)=\sum W_r^2 f_b^4 Z$，按径向速度直方图累加；含距离加权、双程高斯波束、速度折叠 |
| `yu2007_fig1.m` | Fig.1 |
| `yu2007_fig2.m` | Fig.2（PPI + 5×5 含噪谱，va=35 m/s 混叠） |
| `yu2007_fig4.m` | Fig.4 |
| `yu2007_fig5.m` | Fig.5（Td/Tu/Gd，SNR=40 dB） |
| `yu2007_fig678.m` | Fig.6–8（25 个体积组合、$r_n=0.25:0.5:10$、每个点 50 次实现） |

**默认参数（与论文 3a 节一致）**

| 参数 | 值 |
|---|---|
| 龙卷最大切向风 $V_{Tm}$ / 半径 $r_t$ | 50 m/s / 200 m |
| 反射率 $z_m$、$r_{0z}$、$W_z$ | 1、220 m、60 m（doughnut） |
| 波束宽度 $\theta_b$ | 1°（单程半功率） |
| 距离分辨率 / 距离加权标准差 | 250 m / $\sigma_r=87.5$ m |
| 最大不模糊速度 $v_a$ | 65 m/s（Fig.1、Fig.5）；35 m/s（Fig.2、Fig.6–8） |
| 时间序列点数 / 周期图 | 64（Fig.2、5–8）；128（Fig.1 仅画模型谱） |
| 噪底/峰值 | Fig.1：-50/25 dB；Fig.2/5：加性噪声（SNR 40 dB） |

### Yu et al. (2009)

| 文件 | 论文内容 |
|---|---|
| `y09_model_spectrum.m` | 式(1)：双高斯谱 $S_1N(v_1,\sigma_1)+S_2N(v_2,\sigma_2)$ |
| `y09_dual_moments.m` | 附录 A 式(A1)–(A4)：$\langle v_r\rangle,\sigma_v$、偏度、峰度 |
| `y09_theory_ac.m` | 式(2)–(4)：双高斯谱的理论自相关与 AC 偏差 |
| `y09_sim_bias.m` | 200 次实现 × 8 trial 平均 AC 的偏差实验 |
| `y09_fit_dual.m` | 式(5)：dB 域多初值非线性最小二乘，估计 6 个谱矩 + 噪声 |
| `y09_sim_spectrum.m` | 含噪周期图 + AC 矩，供 Fig.2/3 使用 |
| `y09_rmse_one.m` | 单次实现的 AC 重构/双高斯拟合 RMSE |
| `yu2009_fig1.m` | Fig.1 |
| `yu2009_fig2.m` | Fig.2 |
| `yu2009_fig3.m` | Fig.3 |
| `yu2009_fig8.m` | Fig.8 |

**Fig.1 四组实验（Table 1）**

| Expt | $S_2$ (dBm) | $\sigma_2$ (m/s) |
|---|---|---|
| I | 10 | 2 |
| II | 10 | 6 |
| III | 1 | 6 |
| IV | -20 | 4 |

公共：$S_1=10$ dBm，$v_1=10$ m/s，$\sigma_1=2$ m/s，$v_2=-10:1:10$ m/s，$v_a=25$ m/s，$M=64$，SNR=30 dB，200 次实现，每次 8 条独立 IQ 平均 ACF/谱。

**Fig.8 垂直切变**

式(6)：$S(v)=\int f_b^4(\phi)Z(z)|dv/dz|^{-1}dz$；$v_r(z)=20\tanh(z/800)$ m/s，$z_s=800$ m，$N(0,2)$ m/s 扰动；1° 波束，$r=65$ km，$z=\sqrt{r^2+R'^2+2rR'\sin\phi}-R'$，$R'=\frac43\times6371$ km。

---

## 3. 实现要点

### 3.1 谱成形与符号约定

`ysim_iq_from_psd` 先把按速度轴 `-va:dv:va-dv` 排列的 PSD 用 `ifftshift` 映射到 FFT 频率顺序，再取共轭，使正速度对应负多普勒相位。这样：

- 周期图峰值出现在 `v>0`；
- 脉冲对公式 $v=-(v_a/\pi)\angle R_1$ 与周期图速度轴一致；
- 不会出现固定速度偏移或镜像。

IQ 总功率严格等于 `sum(S)`；`ysim_periodogram` 的归一化保证 `sum(S_bin)≈M·mean(|iq|^2)`，转成功率谱密度时除以 `M·dv=2v_a`。

### 3.2 二维体积积分

`y07_sim_volume` 在体积中心附近取自适应网格（默认 241×241，约 5–15 m 间距），计算：

- 双程天线功率方向图 $f_b^4=\exp[-8\ln2(\Delta\theta/\theta_b)^2]$；
- 距离加权 $W_r^2=\exp[-(R-r_0)^2/(2\sigma_r^2)]$；
- 反射率与 Rankine 速度场；
- 将 $W_r^2f_b^4Z$ 按折叠后的 $v_r$ 累加到速度 bin。

为抑制离散网格造成的“梳状”伪影和涡心焦散尖峰，加入：

- `vel_jitter`：每个网格点 $N(0,\sigma_j)$ 的亚网格速度扰动（默认 0.5–2 m/s，按图调节）；
- `smooth_bins`：对直方图做约 1–2 个 bin 的高斯平滑。

这些都是对论文未列出的亚网格湍流/有限速度分辨率的等效处理，不改变谱的一、二阶矩量级。

### 3.3 谱平坦度

论文把 $\sigma_s$ 称为 dB 谱的“统计方差”，但 Fig.3/Fig.7 的数值范围（白噪声约 5.6、龙卷区 5–10）表明实际是**标准差**。本实现取

$$\sigma_s=\mathrm{std}\big(10\log_{10}S(v_k)\big).$$

若需要方差，只需对返回值平方。

### 3.4 PRIB

`ysim_prib` 严格按论文式(3)与 Chandran–Elgar 的离散插值实现；`method` 可选：

- `'absfft'`（默认）：论文 3b 节给出的对称性规避方法，先取 $|X(f)|$ 前一半再作双谱；
- `'direct'`：直接对 dB 谱作双谱。偶对称图案会退化为 0/180°。

`yu2007_fig678` 使用 `'direct'` 并对 50 次实现的相位做算术平均，才能得到论文 Fig.8 的 20–70° 量级；这与论文“每个点是 50 次实现的平均”一致。`yu2007_fig4` 使用 `'absfft'`。

**已知偏差**：论文没有给出 PRIB 的完整预处理/平均代码，Fig.4 的绝对 P 值和 Fig.8 的曲线细节无法逐点复现。本实现的定性趋势（矩形谱 P 高于高斯谱、龙卷区 P 升高）成立，但 Fig.4 小谱宽段（$\sigma_v<8$ m/s）的数值偏低。若需严格对齐，需要向作者索取原始代码。

### 3.5 Fig.6–8 的高斯重构

`yu2007_fig678` 中 $G_d/G_u$ 的高斯谱由 Td/Tu 的 AC 矩构造：

- 中心取**谱质心** $v_m=\sum v_kS_k/\sum S_k$（而不是混叠后的 AC 平均速度，避免高斯被平移到 Nyquist 边缘）；
- 谱宽取 $\min(\max(\hat\sigma_v,10),25)$ m/s。10 m/s 下限对应可分辨谱宽/湍流底，25 m/s 上限避免重构谱退化为全带宽平谱；
- 噪声底与观测谱同 SNR，用解析谱 + 噪声底计算 $\sigma_s$ 和 PRIB，不再额外抽一次 IQ。

这是对论文“用自协方差三矩重构高斯谱”的稳定化处理，目的是避免 AC 在强混叠/窄谱时给出退化宽度；不影响 Td/Tu 的统计结论。

---

## 4. 复现结果与校验

`verify_yu_sim.m` 会打印以下数值校验：

- 2007：25 km 谱的峰值/边缘速度、谱宽、噪声底；
- 2009：四组实验在 $\Delta v=-20,-14,-2$ 的 AC 偏差（理论 vs 200 次仿真）；
- 2009：双高斯拟合在 Expt II/Δv=-14 的参数回收；
- 2009：Fig.3 的 RMSE 均值。

参考结果（本机运行）：

| 项目 | 结果 |
|---|---|
| 2007 Fig.1 25 km 谱边缘 | 约 ±50 m/s（$V_{Tm}=50$） |
| 2007 Fig.1 75 km 谱 | 平顶、中心浅凹、噪底 -50 dB |
| 2009 Expt II 最大速度偏差 | 2.6–3.0 m/s（论文约 2.5） |
| 2009 Expt III 速度偏差 | 约 1.7 m/s |
| 2009 Expt I/IV 速度偏差 | <0.2 m/s |
| 2009 双高斯拟合 RMSE 均值 | 1.40 / 1.71 / 1.55 / 1.40 dB（论文 1.48 / 1.81 / 1.69 / 1.49） |
| 2009 拟合参数（Expt II, Δv=-14） | $S_1≈10$、$S_2≈9$ mW，$v_1≈9.7$、$v_2≈-4.5$ m/s，$\sigma_1≈2.2$、$\sigma_2≈5.8$ m/s |

---

## 5. 主要简化与差异

1. **5 km 谱形**：论文 Fig.1a 的 5 km 曲线是边缘约 ±20 m/s 的平顶矩形；本实现给出钟形（体积只采样 Rankine 核心的线性段，且高斯波束加权）。这是体积/波束参数的定义差异，25/75 km 的平顶、边缘 ±50 m/s 与中心凹陷均能复现。
2. **Fig.1d 反射率类型**：为同时复现 $W_z=20$ 的钟形、$W_z=60$ 的平顶+中心凹陷、$W_z=\infty$ 的矩形，分别采用 Gaussian（$r_{0z}=0$）、doughnut（$r_{0z}=220$ m）、uniform。论文未明确说明 $r_{0z}$ 是否随 $W_z$ 切换，这是按曲线形态反推的等效设置。
3. **Fig.4/Fig.8 PRIB**：见 3.4 节。Fig.4 的绝对 P 值未完全对齐；Fig.8 的 P 曲线仍保留了逐点实现的算术平均噪声，未做论文中可能使用的额外平滑/门限处理，因此曲线比论文更起伏。$\sigma_{vn}$ 与 $\sigma_s$（Fig.6/7）的趋势和量级与论文一致。
4. **Fig.5/6–8 的亚网格扰动**：`vel_jitter`、`smooth_bins` 不是论文显式参数，用于补偿 241×241 网格与 64 点谱的离散化；取值已在各脚本中给出。
5. **2009 Fig.1 的仿真**：采用 8 条 IQ 的**周期**自相关平均；论文未说明是否用循环滞后。周期滞后可消除有限样本对谱宽的系统性低估。
6. **Fig.8 反射率廓线**：论文只给图，未给解析式；本实现用双高斯层 + 局部凹陷（Case I）和低层高斯 + 高层指数衰减（Case II）拟合其形态，能产生 0.5° 双峰、1.5° 宽裙。

---

## 6. 目录

```
Non-Gaussian-spectrum-research/  (仓库根目录)
  README.md               本说明
  MANIFEST.md             文件清单、图件映射、校验和
  CITATION.md             论文引用
  run_all.m               一键运行全部图件
  verify_yu_sim.m         数值校验
  ysim_*.m                公共谱/统计函数（6 个）
  y07_*.m                 Yu2007 场与体积积分（2 个）
  y09_*.m                 Yu2009 双高斯/拟合/实验（7 个）
  yu2007_*.m              Yu2007 图件脚本（5 个）
  yu2009_*.m              Yu2009 图件脚本（4 个）
  figs/                   12 张输出图件（PNG）
  data/                   9 个中间结果（MAT，可删除后重跑）
  docs/verification_report.md         本机运行结果与阈值
```

代码保持在根目录平铺，是为了让所有脚本在任意工作目录下都能通过
`addpath(fileparts(mfilename('fullpath')))` 找到彼此；`figs/`、`data/`
由脚本自动创建/覆盖。

论文文本与页面截图放在原工作区的 `papers_text/`，仅用于核对，不随
本压缩包分发。

---

## 7. 引用

使用本代码时请引用原始论文（完整条目见 `CITATION.md`）：

- Yu, T.-Y., Y. Wang, A. Shapiro, M. B. Yeary, D. S. Zrnić, and R. J. Doviak, 2007: *Characterization of Tornado Spectral Signatures Using Higher-Order Spectra*. J. Atmos. Oceanic Technol., **24**, 1997–2013.
- Yu, T.-Y., R. Reinoso-Rondinel, and R. D. Palmer, 2009: *Investigation of Non-Gaussian Doppler Spectra Observed by Weather Radar in a Tornadic Supercell*. J. Atmos. Oceanic Technol., **26**, 444–461.

本目录为独立重实现，不是作者原始代码；PRIB、Fig.1d、Fig.6–8 高斯
重构等处的等效处理已在第 5 节说明。

