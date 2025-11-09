%% 
clc
clear
close all

%% Set working folder
restoredefaultpath;
main_folder = pwd;
addpath([main_folder,'\functions'])

%% Load workspace
load("workspace\p2d_workspace_LIB4.mat")
n1 = mesh.n1;
n2 = mesh.n2;
n3 = mesh.n3;
n4 = length(MatC.Lcn);
ord = mesh.ord;

%% Initial solid concentration
V0 = median([data.V_V(1) - 0.1, min(SOCTab.OCV), max(SOCTab.OCV)]);
[~, index_unique] = unique(SOCTab.OCV);
xn_init = interp1(SOCTab.OCV(index_unique), SOCTab.xn(index_unique), V0);
xp_init = interp1(SOCTab.OCV(index_unique), SOCTab.xp(index_unique), V0);

%% Input signal
dt = 0.5;
t_end = sim_ref.t_s(end);

t0 = data.t_s;
I0 = data.C_rate*Const.I_1C;
index_CV = find(data.Control_mode==2);
I0(index_CV) = I0(index_CV(1)-1);
t_eoc = t0(index_CV(end));

t_sim = zeros(1);
I_data = zeros(1);

ii = 1;
tt = t0(ii);
t_sim(ii) = tt;
I_data(ii) = I0(ii);

while tt<t0(end)
    ii = ii + 1;
    tt = tt + dt;
    t_sim(ii,1) = tt;
    index_cut = find(t0<=tt,1,"last");
    I_data(ii,1) = I0(index_cut);
end
T_data = repmat(298, size(t_sim));

I_sim = timeseries(I_data, t_sim);
T_sim = timeseries(T_data, t_sim);
V_max = 4.35;

%% Run simulation
Result = sim("Li_ion_cell_p2d_with_control_R2021b.slx");

%% Postprocessing and plotting
figure 
ax =gca;
line(Result.tout/60, Result.Vex.data, 'parent',ax,'linewidth',2,'color','r')
line(data.t_s/60,data.V_V, 'parent',ax,'linestyle',':','linewidth',4,'color','b');
xlabel("Time [mim]")
ylabel("Voltage [V]")
legend("Model", "Data")
title("Figure 5(a) of Ref 1")

figure 
ax =gca;
line(Result.tout/60, Result.I_control.data, 'parent',ax,'linewidth',2,'color','r')
line(data.t_s/60,data.C_rate*Const.I_1C, 'parent',ax,'linestyle',':','linewidth',4,'color','b');
xlabel("Time [mim]")
ylabel("Current [A]")
legend("Model", "Data")
title("Figure 5(b) of Ref 1")

xn_2nd = linspace(0, mesh.Ln, ord*n1+1);
xs_2nd = linspace(mesh.Ln, mesh.Ln+mesh.Ls, ord*n2+1);
xp_2nd = linspace(mesh.Ln+mesh.Ls, mesh.Ln+mesh.Ls+mesh.Lp, ord*n3+1);
xcoord_2nd = unique([xn_2nd, xs_2nd, xp_2nd])/1e-6;

t_plot = 1200;
figure 
ax =gca;
cl_plot = Result.cl.data(find(Result.tout<=t_plot,1,"last"),:);
line(xcoord_2nd, cl_plot, 'parent',ax,'linewidth',1,'Marker','+','color','r')
line([mesh.Ln,mesh.Ln]/1e-6, [min(cl_plot), max(cl_plot)],'parent',ax,'linewidth',0.5, 'LineStyle','--', 'color','k')
line([mesh.Ln+mesh.Ls,mesh.Ln+mesh.Ls]/1e-6, [min(cl_plot), max(cl_plot)],'parent',ax,'linewidth',0.5, 'LineStyle','--', 'color','k')
xlabel("x [um]")
ylabel("Electrolyte concentration [mol/m^3]")
title(['Electrolyte concentration at t = ', num2str(t_plot), 's'])

t_plot = 60;
figure 
ax =gca;
phil_plot = Result.phil.data(find(Result.tout<=t_plot,1,"last"),:);
line(xcoord_2nd, phil_plot, 'parent',ax,'linewidth',1,'Marker','+','color','r')
line([mesh.Ln,mesh.Ln]/1e-6, [min(phil_plot), max(phil_plot)],'parent',ax,'linewidth',0.5, 'LineStyle','--', 'color','k')
line([mesh.Ln+mesh.Ls,mesh.Ln+mesh.Ls]/1e-6, [min(phil_plot), max(phil_plot)],'parent',ax,'linewidth',0.5, 'LineStyle','--', 'color','k')
xlabel("x [um]")
ylabel("Electrolyte potential [V]")
title(['Electrolyte potential at t = ', num2str(t_plot), 's'])

xn_1st = linspace(0, mesh.Ln, n1+1);
xs_1st = linspace(mesh.Ln, mesh.Ln+mesh.Ls, n2+1);
xp_1st = linspace(mesh.Ln+mesh.Ls, mesh.Ln+mesh.Ls+mesh.Lp, n3+1);
xcoord_1st = unique([xn_1st, xs_1st, xp_1st])/1e-6;

t_plot = 1200;
js_orig = Result.js.data(find(Result.tout<=t_plot,1,"last"),:);
js_plot = nan(size(xcoord_1st));
js_plot(1:n1+1) = [js_orig(1), js_orig(1:n1-1)*0.5 + js_orig(2:n1)*0.5, js_orig(n1)];
js_plot(n1+n2+1:n1+n2+n3+1) =[js_orig(n1+1), js_orig(n1+1:n1+n3-1)*0.5 + js_orig(n1+2:n1+n3)*0.5, js_orig(n1+n3)];

if n1>=3
    js_plot(1) = 3/2*js_orig(1) - 1/2*js_orig(2);
    js_plot(n1+1) = 3/2*js_orig(n1) - 1/2*js_orig(n1-1);
end

if n3>=3
    js_plot(n1+n2+1) = 3/2*js_orig(n1+1) - 1/2*js_orig(n1+2);
    js_plot(n1+n2+n3+1) = 3/2*js_orig(n1+n3) - 1/2*js_orig(n1+n3-1);
end

figure 
ax =gca;
line(xcoord_1st, js_plot, 'parent',ax,'linewidth',0.5,'marker','o','color','b')
line([mesh.Ln,mesh.Ln]/1e-6, [min(js_plot([1:n1+1,n1+n2+1:n1+n2+n3+1])), max(js_plot([1:n1+1,n1+n2+1:n1+n2+n3+1]))],'parent',ax,'linewidth',0.5, 'LineStyle','--', 'color','k')
line([mesh.Ln+mesh.Ls,mesh.Ln+mesh.Ls]/1e-6, [min(js_plot([1:n1+1,n1+n2+1:n1+n2+n3+1])), max(js_plot([1:n1+1,n1+n2+1:n1+n2+n3+1]))],'parent',ax,'linewidth',0.5, 'LineStyle','--', 'color','k')
xlabel("x [um]")
ylabel("Reaction current [A/m^2]")
title(['Reaction current at t = ', num2str(t_plot), 's'])


