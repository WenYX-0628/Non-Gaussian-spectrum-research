# 引用与说明

本代码是对以下两篇论文中数值仿真的独立重实现：

1. **Yu, T.-Y., Y. Wang, A. Shapiro, M. B. Yeary, D. S. Zrnić, and R. J. Doviak, 2007**:
   *Characterization of Tornado Spectral Signatures Using Higher-Order Spectra.*
   Journal of Atmospheric and Oceanic Technology, **24**(12), 1997–2013.
   DOI: [10.1175/2007JTECHA934.1](https://doi.org/10.1175/2007JTECHA934.1)

2. **Yu, T.-Y., R. Reinoso-Rondinel, and R. D. Palmer, 2009**:
   *Investigation of Non-Gaussian Doppler Spectra Observed by Weather Radar in a Tornadic Supercell.*
   Journal of Atmospheric and Oceanic Technology, **26**(3), 444–461.
   DOI: [10.1175/2008JTECHA1124.1](https://doi.org/10.1175/2008JTECHA1124.1)

## 说明

- 本代码不是论文作者的原始代码，而是根据论文正文、附录和图形独立编写的。
- 论文未明确给出的细节（PRIB 预处理/平均、Fig.1d 的反射率类型切换、Fig.8 的反射率廓线解析式等）采用了等效或稳定化处理，具体差异见 `README.md` 第 5 节。
- 若用于学术发表，请引用上述原始论文，并在方法描述中注明“independent reimplementation”。

## 第三方参考

PRIB 的离散插值方式参考了：

- Chandran, V., and S. L. Elgar, 1993: *Pattern recognition using invariants defined from higher order spectra—One-dimensional inputs.* IEEE Trans. Signal Process., **41**(1), 205–212. DOI: [10.1109/TSP.1993.193139](https://doi.org/10.1109/TSP.1993.193139)

其作者公开的示例代码（`bispinvS.m`）仅用于核对插值公式，未直接包含在本包中。
