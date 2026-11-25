function [S, v, meta] = y07_sim_volume(opt)
%Y07_SIM_VOLUME  由二维反射率/风场积分得到雷达分辨率体积内的多普勒谱
%
%   [S, v, meta] = y07_sim_volume(opt)
%
%   实现 Yu et al. (2007, JTECH) 附录式(A1)-(A5)：
%       S(r0,vk) = sum_{i,j} Wr^2(i,j) fb^4(i,j) Z(i,j)
%       仅累加 vk-dv/2 <= vr(i,j) < vk+dv/2 的网格点
%
%   opt 字段（均有默认值）
%     .r0        体积中心距雷达距离 [m]              (25000)
%     .az0       波束方位角 [deg]（0=正北）          (0)
%     .va        最大不模糊速度 [m/s]                (65)
%     .M         速度 bin 数                         (128)
%     .beamwidth 单程半功率波束宽度 [deg]            (1)
%     .sigma_r   距离加权函数标准差 [m]              (87.5 = 0.35*250m)
%     .gridN     网格点数（每维）                    (241)
%     .vortices  涡旋数组（见 y07_fields）           (单个龙卷)
%     .u0/.v0    背景风 [m/s]                        (0/0)
%     .shear_dvdx 水平切变 [1/s]                     (0)
%     .refl      反射率结构（见 y07_fields）          (doughnut, Wz=60)
%     .returnFields 是否返回二维场                   (false)
%
%   输出
%     S    : [M x 1] 模型谱（未归一化功率）
%     v    : [M x 1] bin 中心速度 [m/s]
%     meta : 网格/场/权重等信息

if nargin < 1 || isempty(opt), opt = struct(); end
def = struct('r0',25000,'az0',0,'va',65,'M',128,'beamwidth',1, ...
             'sigma_r',87.5,'gridN',241,'u0',0,'v0',0,'shear_dvdx',0, ...
             'returnFields',false);
fn = fieldnames(def);
for i = 1:numel(fn)
    if ~isfield(opt,fn{i}) || isempty(opt.(fn{i})), opt.(fn{i}) = def.(fn{i}); end
end
if ~isfield(opt,'vortices') || isempty(opt.vortices)
    opt.vortices = struct('x',0,'y',opt.r0,'rt',200,'vtm',50,'vrm',0);
end
if ~isfield(opt,'refl') || isempty(opt.refl)
    opt.refl = struct('type','doughnut','Wz',60,'r0z',220,'zm',1);
end

r0   = opt.r0;
az0  = deg2rad(opt.az0);
thb  = deg2rad(opt.beamwidth);
cx   = r0*sin(az0);
cy   = r0*cos(az0);

% 网格范围：覆盖波束（横向）与距离加权（纵向），同时覆盖涡旋核心
rmax = 200;
for k = 1:numel(opt.vortices), rmax = max(rmax, opt.vortices(k).rt); end
half_x = max(2.4*rmax, 2.2*r0*thb);
half_y = max(2.6*rmax, 4.0*opt.sigma_r);

N = opt.gridN;
x = cx + linspace(-half_x, half_x, N);
y = cy + linspace(-half_y, half_y, N);
[X, Y] = meshgrid(x, y);

cfg = opt;
cfg.shear_x0 = cx;
[Z, U, V] = y07_fields(X, Y, cfg);

% scene 反射率以 dBZ 给出，积分时换算成线性功率
if isfield(opt.refl,'Z_is_dBZ') && opt.refl.Z_is_dBZ
    Z = 10.^(Z/10);
end

% 径向速度：正=远离雷达
R  = sqrt(X.^2 + Y.^2);
R  = max(R, 1);
vr = (U.*X + V.*Y)./R;

% 小尺度速度扰动（湍流/亚网格运动）：打散 Rankine 涡心与等速线在
% 直方图中的焦散尖峰，使模拟谱接近论文中光滑的连续谱。
if ~isfield(opt,'vel_jitter') || isempty(opt.vel_jitter)
    opt.vel_jitter = 0.0;
end
if opt.vel_jitter > 0
    vr = vr + opt.vel_jitter*randn(size(vr));
end

% 速度折叠（论文附录：在谱域模拟速度混叠）
va = opt.va;
vr = mod(vr + va, 2*va) - va;

% 双程天线方向图与距离加权
dth = atan2(X, Y) - az0;
dth = mod(dth + pi, 2*pi) - pi;
fb4 = exp(-8*log(2)*(dth/thb).^2);          % 双程功率方向图
wr2 = exp(-((R - r0).^2)/(2*opt.sigma_r^2)); % 距离加权（功率）

W = wr2 .* fb4 .* max(Z,0);

% 直方图累加到速度 bin
M  = opt.M;
edges = linspace(-va, va, M+1);
dv = edges(2)-edges(1);
v  = edges(1:end-1).' + dv/2;
idx = floor((vr(:) + va)/dv) + 1;
idx = min(max(idx,1),M);
S = accumarray(idx, W(:), [M 1]);

% 网格离散化会把连续的速度分布采样成“梳状”；用约 1 个 bin 的高斯核
% 平滑，等效于有限速度分辨率/湍流展宽，也避免空 bin。
if ~isfield(opt,'smooth_bins') || isempty(opt.smooth_bins)
    opt.smooth_bins = 1.0;
end
if opt.smooth_bins > 0
    sb = opt.smooth_bins;
    nk = max(3, 2*ceil(3*sb)+1);
    kk = (-(nk-1)/2:(nk-1)/2)';
    g  = exp(-0.5*(kk/sb).^2); g = g/sum(g);
    S  = conv(S, g, 'same');
end

meta = struct();
meta.r0 = r0; meta.az0 = opt.az0; meta.va = va; meta.M = M;
meta.dv = dv; meta.cx = cx; meta.cy = cy;
meta.half_x = half_x; meta.half_y = half_y;
meta.vortices = opt.vortices; meta.refl = opt.refl;
meta.SumW = sum(W(:));
if opt.returnFields
    meta.X = X; meta.Y = Y; meta.Z = Z; meta.U = U; meta.V = V; meta.W = W;
end
end
