function [Z, U, V, rv] = y07_fields(X, Y, cfg)
%Y07_FIELDS  Yu et al. (2007) 二维反射率场与水平风场
%
%   [Z, U, V, rv] = y07_fields(X, Y, cfg)
%
%   X,Y : 笛卡尔网格（m），雷达位于原点，x 向东、y 向北
%   cfg : 结构体
%       .vortices  : 涡旋结构体数组，字段 x,y,rt,vtm,vrm
%                    （Combined Rankine：r<=rt 时 V=Vm*(r/rt)，r>rt 时
%                     V=Vm*(rt/r)，对应论文式(A4)(A5)，gamma=+1/-1）
%       .u0,.v0    : 背景风 [m/s]
%       .shear_dvdx: 水平切变 dv/dx [1/s]（论文图1c）
%       .shear_x0  : 切变参考点 x0（通常取雷达体积中心）
%       .refl      : 反射率结构体
%           .type   : 'doughnut' | 'uniform' | 'gaussian' | 'scene'
%           .Wz     : 反射率宽度 [m]（式(A2)）
%           .r0z    : 最大反射率半径 [m]（式(A2)）
%           .zm     : 最大反射率（线性功率/相对值）
%           .center : [x y] 龙卷中心（默认第一个涡旋中心）
%           .bg     : 均匀背景反射率（默认 0）
%           .scene  : 'scene' 类型的参数（龙卷+中气旋）
%
%   输出 Z 为线性反射率（相对值），U/V 为水平风分量 [m/s]。
%   rv 为第一个涡旋（龙卷）中心的距离场，便于绘图。

Z = zeros(size(X));
U = zeros(size(X));
V = zeros(size(X));

if isfield(cfg,'u0') && ~isempty(cfg.u0), U(:) = cfg.u0; end
if isfield(cfg,'v0') && ~isempty(cfg.v0), V(:) = cfg.v0; end

% ---- 水平切变：V = V + (dv/dx)(x-x0) ----
if isfield(cfg,'shear_dvdx') && abs(cfg.shear_dvdx) > 0
    x0 = 0;
    if isfield(cfg,'shear_x0') && ~isempty(cfg.shear_x0), x0 = cfg.shear_x0; end
    V = V + cfg.shear_dvdx * (X - x0);
end

% ---- 叠加一个或多个 Combined Rankine 涡旋 ----
rv = inf(size(X));
if isfield(cfg,'vortices') && ~isempty(cfg.vortices)
    for k = 1:numel(cfg.vortices)
        vo = cfg.vortices(k);
        dx = X - vo.x;
        dy = Y - vo.y;
        r  = sqrt(dx.^2 + dy.^2);
        r  = max(r, eps);
        gam = ones(size(r));
        gam(r > vo.rt) = -1;                 % 式(A4)(A5) 的 gamma
        prof = (r/vo.rt).^gam;               % r/rt (内) 或 rt/r (外)
        vt = vo.vtm * prof;                  % 切向速度
        vr = 0;
        if isfield(vo,'vrm') && ~isempty(vo.vrm), vr = vo.vrm * prof; end
        % 逆时针切向 + 径向外流
        U = U - vt .* dy ./ r + vr .* dx ./ r;
        V = V + vt .* dx ./ r + vr .* dy ./ r;
        if k == 1, rv = r; end
    end
end

% ---- 反射率 ----
rf = struct('type','uniform','Wz',1e9,'r0z',0,'zm',1,'bg',0);
if isfield(cfg,'refl') && ~isempty(cfg.refl)
    fn = fieldnames(cfg.refl);
    for i = 1:numel(fn), rf.(fn{i}) = cfg.refl.(fn{i}); end
end
if isinf(rv)
    if isfield(rf,'center') && ~isempty(rf.center)
        rv = sqrt((X-rf.center(1)).^2 + (Y-rf.center(2)).^2);
    else
        rv = sqrt(X.^2 + Y.^2);
    end
end

switch lower(rf.type)
    case 'uniform'
        Z = rf.zm * ones(size(X));
    case 'doughnut'
        Z = rf.zm * exp(-((rv - rf.r0z)./rf.Wz).^2);
    case 'gaussian'
        Z = rf.zm * exp(-(rv./rf.Wz).^2);
    case 'scene'
        % 龙卷（doughnut）+ 中气旋（Gaussian）；用于 Fig.2
        sc = rf.scene;
        rvT = sqrt((X - sc.tornado(1)).^2 + (Y - sc.tornado(2)).^2);
        rvM = sqrt((X - sc.meso(1)).^2    + (Y - sc.meso(2)).^2);
        Zt = sc.zt * exp(-((rvT - sc.tornado_r0z)./sc.tornado_Wz).^2);
        Zm = sc.zm * exp(-(rvM./sc.meso_Wz).^2);
        Z  = max(sc.bg, max(Zt, Zm));
        % 龙卷中心低反射率孔（图2a 的暗洞）
        Z  = Z .* (1 - sc.holeDepth*exp(-(rvT./sc.holeRadius).^2));
        rv = rvT;
    otherwise
        error('y07_fields: unknown reflectivity type "%s".', rf.type);
end

Z = Z + rf.bg * ones(size(X));   % 允许额外均匀背景
if strcmpi(rf.type,'scene'), Z = Z - rf.bg; end
end
