function yu2007_fig2()
%YU2007_FIG2  Yu et al. (2007) Fig.2：龙卷+中气旋的反射率/速度场与 5x5 谱
%
%   输出
%     figs/yu2007_fig2a.png  模型反射率与虚拟 WSR-88D 径向速度
%     figs/yu2007_fig2b.png  黑框内 5 个方位 x 5 个距离门的含噪多普勒谱
%
%   场景：龙卷位于 (0,25 km)，中气旋中心在东北方向 1 km；VTm=50/25 m/s，
%   rt=200/2000 m，VRm=5 m/s；背景风 (u0,v0)=(5,10) m/s，dv/dx=0.01 s^-1；
%   va=35 m/s（发生速度混叠），64 点周期图。

here = fileparts(mfilename('fullpath'));
addpath(here);
figdir = fullfile(here,'figs'); datadir = fullfile(here,'data');
if ~exist(figdir,'dir'), mkdir(figdir); end
if ~exist(datadir,'dir'), mkdir(datadir); end
rng(2007,'twister');

va = 35;
vort = [struct('x',0,'y',25000,'rt',200,'vtm',50,'vrm',5), ...
        struct('x',1000*sin(pi/4),'y',25000+1000*cos(pi/4),'rt',2000,'vtm',25,'vrm',5)];
scene = struct('tornado',[0 25000],'meso',[1000*sin(pi/4),25000+1000*cos(pi/4)], ...
    'zt',45,'tornado_r0z',220,'tornado_Wz',60, ...
    'zm',45,'meso_Wz',2000,'bg',5,'holeDepth',0.98,'holeRadius',180);
cfg = struct('vortices',vort,'u0',5,'v0',10,'shear_dvdx',0.01, ...
             'refl',struct('type','scene','scene',scene));

% ---------------- PPI ----------------
xg = linspace(-1000,3500,260); yg = linspace(23000,30000,280);
[X,Y] = meshgrid(xg,yg);
[Z,U,V] = y07_fields(X,Y,cfg);
R = sqrt(X.^2+Y.^2); R = max(R,1);
vr = (U.*X + V.*Y)./R;
vr = mod(vr + va, 2*va) - va;      % 速度混叠

fig = figure('Color','w','Position',[80 80 980 460]);
tiledlayout(1,2,'Padding','compact','TileSpacing','compact');
nexttile;
imagesc(xg/1000,yg/1000,Z); axis xy equal tight; clim([-10 50]);
colormap(gca,parula); cb=colorbar; cb.Label.String='Reflectivity (dBZ)';
xlabel('x-distance (km)'); ylabel('y-distance (km)');
title('(a) Model reflectivity and horizontal wind');
hold on;
% 黑色取样框
plot([-0.5 0.5 0.5 -0.5 -0.5],[24.5 24.5 25.5 25.5 24.5],'k-','LineWidth',1.4);
quiver(X(1:20:end,1:20:end)/1000,Y(1:20:end,1:20:end)/1000, ...
       U(1:20:end,1:20:end),V(1:20:end,1:20:end),0.8,'w');

nexttile;
imagesc(xg/1000,yg/1000,vr); axis xy equal tight; clim([-va va]);
colormap(gca,jet); cb=colorbar; cb.Label.String='v_r (m s^{-1})';
xlabel('x-distance (km)'); ylabel('y-distance (km)');
title('(b) Radial velocity observed by virtual WSR-88D');
hold on;
plot([-0.5 0.5 0.5 -0.5 -0.5],[24.5 24.5 25.5 25.5 24.5],'k-','LineWidth',1.4);
sgtitle('Yu et al. (2007) Fig. 2a,b reproduction','FontWeight','bold');
saveas(fig,fullfile(figdir,'yu2007_fig2a.png'));

% ---------------- 5x5 spectra ----------------
azList = [-2.2 -1.2 -0.2 0.8 1.8];
rList  = [24.5 24.75 25.0 25.25 25.5]*1000;
M = 64;
base = struct('va',va,'M',M,'beamwidth',1,'sigma_r',87.5,'gridN',241, ...
    'vel_jitter',0.5,'smooth_bins',1.0,'u0',5,'v0',10,'shear_dvdx',0.01, ...
    'vortices',vort, ...
    'refl',struct('type','scene','scene',scene,'Z_is_dBZ',true));

fig2 = figure('Color','w','Position',[60 40 980 920]);
tiledlayout(5,5,'Padding','tight','TileSpacing','tight');
Sshow = cell(5,5); vac = zeros(5,5); vax = [];
for ir = 1:5
    for ia = 1:5
        o = base; o.r0 = rList(ir); o.az0 = azList(ia);
        [S,v] = y07_sim_volume(o);
        % 峰值 25 dB、噪底 -55 dB 的模型谱
        Sn = S/max(S)*10^(2.5) + 10^(-5.5);
        iq = ysim_iq_from_psd(Sn, M, Inf, 1, 200000+100*ir+ia);
        [Sb, vb] = ysim_periodogram(iq, va, M, 'hann');
        Sshow{ir,ia} = Sb/max(Sb);
        am = ysim_ac_moments(iq, va);
        vac(ir,ia) = am.v;
        vax = vb;

        nexttile((ir-1)*5+ia); hold on; box on;
        plot(vb, 10*log10(max(Sshow{ir,ia},10^(-75/10))) + 25, 'k-', 'LineWidth',0.8);
        yl = [-60 40]; xlim([-va va]); ylim(yl);
        plot([vac(ir,ia) vac(ir,ia)],yl,'k:','LineWidth',0.5);
        if vac(ir,ia) >= 0
            plot(vac(ir,ia), 28, '^k', 'MarkerFaceColor','k','MarkerSize',3);
        else
            plot(vac(ir,ia), 28, 'vk', 'MarkerFaceColor','k','MarkerSize',3);
        end
        set(gca,'XTick',-30:30:30,'YTick',-60:30:30,'FontSize',6);
        if ir==1, title(sprintf('\\phi_0=%.1f^o',azList(ia)),'FontSize',8); end
        if ia==1, ylabel(sprintf('%.2f km',rList(ir)/1000),'FontSize',7); end
        if ir==5, xlabel('v_r','FontSize',7); end
        if ia>1, set(gca,'YTickLabel',[]); end
        if ir<5, set(gca,'XTickLabel',[]); end
    end
end
sgtitle('Yu et al. (2007) Fig. 2 (bottom) reproduction','FontWeight','bold');
saveas(fig2,fullfile(figdir,'yu2007_fig2b.png'));

save(fullfile(datadir,'yu2007_fig2.mat'),'Sshow','vax','vac','azList','rList','Z','U','V','xg','yg');
end
