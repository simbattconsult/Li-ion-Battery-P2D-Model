clc
clear 
close all

%% Set working folder
restoredefaultpath;
main_folder = pwd;
addpath([main_folder,'\functions'])

%% Load worksapce
load("workspace\p2d_workspace_LIB4.mat")

%% Run simulation
[ Result ] = P2D_run( mesh, Const, MatC, SOCTab, data  );

%% Plot results
figure 
ax =gca;
line(Result(:,1)/60, Result(:,2), 'parent',ax,'linewidth',2,'color','r')
line(data.t_s/60,data.V_V, 'parent',ax,'linestyle',':','linewidth',4,'color','b');
xlabel("Time [mim]")
ylabel("Voltage [V]")
legend("Model", "Data")

figure 
ax =gca;
line(Result(:,1)/60, Result(:,3)*Const.Ac, 'parent',ax,'linewidth',2,'color','r')
line(data.t_s/60,data.C_rate*Const.I_1C, 'parent',ax,'linestyle',':','linewidth',4,'color','b');
xlabel("Time [mim]")
ylabel("Current [A]")
legend("Model", "Data")

