function yu2009_fig8()
%YU2009_FIG8  Yu et al. (2009) Fig.8：垂直切变产生双峰/宽裙谱的一维波束积分
%
%   实现式(6)：
%       S(v) = ∫ fb^4(phi) Z(z) |dv/dz|^{-1} dz
%   速度廓线 v_r(z)=v_s tanh(z/z_s)，v_s=20 m/s，z_s=800 m；
%   1° 双程高斯波束，r=65 km，地球曲率高度公式；
%   两组理想反射率廓线（Case I 高层仍有回波，Case II 高层迅速衰减）；
%   每个高度叠加 N(0,2 m/s) 速度扰动。

here = fileparts(mfilename('fullpath'));
addpath(here);
figdir = fullfile(here,'figs'); datadir = fullfile(here,'data');
if ~exist(figdir,'dir'), mkdir(figdir); end
if ~exist(datadir,'dir'), mkdir(datadir); end
rng(2009,'twister');

r = 65000; Rp = 4/3*6371e3; thb = 1; vs = 20; zs = 800;
va = 30; dv = 0.5;
edges = -va:dv:va; ctr = edges(1:end-1)+dv/2;

zfun = @(phi) sqrt(r^2 + Rp^2 + 2*r*Rp*sind(phi)) - Rp;
Z1 = @(z) -5 + 24*exp(-((z-350)/170).^2) + 22*exp(-((z-1200)/330).^2) ...
          - 8*exp(-((z-720)/120).^2);
Z2 = @(z) -5 + 24*exp(-((z-350)/190).^2) - 13*max(z-800,0)/1000;

% ---------------- 谱积分 ----------------
Sspec = cell(2,2);
Smax = 0;
p0List = [0.5 1.5];
for c = 1:2
    for ie = 1:2
        p0 = p0List(ie);
        phi = linspace(p0-2.2, p0+2.2, 2201);
        z = zfun(phi);
        if c==1, Z = Z1(z); else, Z = Z2(z); end
        fb = exp(-8*log(2)*((phi-p0)/thb).^2);
        W = fb .* 10.^(Z/10);
        S = zeros(numel(ctr),1);
        for k = 1:numel(phi)
            vv = vs*tanh(z(k)/zs) + 2*randn(1,21);
            idx = floor((vv+va)/dv)+1;
            idx = min(max(idx,1),numel(ctr));
            S = S + accumarray(idx(:), W(k)/21, [numel(ctr) 1]);
        end
        % 轻微平滑，压制离散直方图锯齿
        g = exp(-0.5*((-3:3)'/1.2).^2); g=g/sum(g);
        S = conv(S,g,'same');
        Smax = max(Smax, max(S));
        Sspec{c,ie} = S/max(S);   % 各谱独立归一化（论文各曲线的绝对功率未标定）
    end
end
save(fullfile(datadir,'yu2009_fig8.mat'),'Sspec','ctr','r','Rp','vs','zs');

% ---------------- 绘图 ----------------
fig = figure('Color','w','Position',[70 60 1100 620]);
tiledlayout(2,6,'Padding','compact','TileSpacing','compact');

% (top-left) velocity profile
nexttile([1 2]); hold on; box on; grid on;
zpl = linspace(0,3000,301);
plot(vs*tanh(zpl/zs), zpl/1000,'k','LineWidth',1.4);
ylim([0 3.2]); xlim([0 22]);
xlabel('V_r (m s^{-1})'); ylabel('Height (km)');
title('Velocity profile');

% (top-middle) reflectivity profiles
nexttile([1 2]); hold on; box on; grid on;
plot(Z1(zpl), zpl/1000,'k-','LineWidth',1.4);
plot(Z2(zpl), zpl/1000,'k--','LineWidth',1.4);
ylim([0 3.2]); xlim([-25 20]);
xlabel('Model Z (dBZ)'); ylabel('Height (km)');
title('Reflectivity profiles');
legend({'Case I','Case II'},'Location','southwest','FontSize',7); legend boxoff;

% (top-right) effective weighting
nexttile([1 2]); hold on; box on; grid on;
styles = {'-','--',':','-.'};
colw = lines(4);
kk=0;
for c=1:2
    for ie=1:2
        kk=kk+1; p0=p0List(ie);
        phi=linspace(p0-2.2,p0+2.2,2201); z=zfun(phi);
        if c==1, Z=Z1(z); else, Z=Z2(z); end
        fb=exp(-8*log(2)*((phi-p0)/thb).^2);
        W=fb.*10.^(Z/10);
        plot(W, z/1000, styles{kk}, 'Color', colw(kk,:), 'LineWidth',1.2);
    end
end
ylim([0 3.2]);
xlabel('W (linear)'); ylabel('Height (km)');
title('Effective weighting');
legend({'I,0.5^\circ','I,1.5^\circ','II,0.5^\circ','II,1.5^\circ'}, ...
       'Location','northeast','FontSize',6); legend boxoff;

% (bottom) spectra
for c=1:2
    nexttile([1 3]); hold on; box on; grid on;
    for ie=1:2
        S = Sspec{c,ie};
        Sdb = 10*log10(max(S,10^(-55/10))) - 20;
        if ie==1
            plot(ctr,Sdb,'k-','LineWidth',1.5);
        else
            plot(ctr,Sdb,'k-.','LineWidth',1.4);
        end
    end
    xlim([-30 30]); ylim([-60 5]);
    xlabel('V_r (m s^{-1})'); ylabel('PSD (dB)');
    title(sprintf('Case %s', repmat('I',1,c)));
    legend({'S(v) @ 0.5^\circ','S(v) @ 1.5^\circ'},'Location','northwest','FontSize',7);
    legend boxoff;
end
sgtitle('Yu et al. (2009) Fig. 8 reproduction','FontWeight','bold');
saveas(fig,fullfile(figdir,'yu2009_fig8.png'));
end
