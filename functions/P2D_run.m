function [ Result ] = P2D_run( mesh, Const, MatC, SOCTab, data )
%%%
n1 = mesh.n1;
n2 = mesh.n2;
n3 = mesh.n3;
ord = mesh.ord;

V0 = data.V_V(1) - 0.1;
us = zeros(length(MatC.Lcn),n1+n3);
us(end,:) = [ ones(n1,1)*Interp_lin( V0, SOCTab.OCV, SOCTab.xn); ones(n3,1)*Interp_lin( V0, SOCTab.OCV, SOCTab.xp ) ];
xs = transpose(sum(us,1));
xl = ones((n1+n2+n3)*ord+1,1);
ul = MatC.Vcl\xl;
T = 298;
II =1;
if data.Control_mode(II)==1
    Crate = data.C_rate(II);
    [ y ] = BVsolve_CC( mesh, Const, MatC, xl , xs, T ,Crate );
elseif data.Control_mode(II)==2
    Vapp = data.V_V(II);
    [ y ] = BVsolve_VC( mesh, Const, MatC, xl , xs, T ,Vapp );
elseif data.Control_mode(II)==3
    Prate = data.C_rate(II);
    [ y  ] = BVsolve_PC( mesh, Const, MatC, xl , xs, T ,Prate );
elseif data(II,4)==4
    Rrate = data.C_rate(II);
    [ y  ] = BVsolve_RC( mesh, Const, MatC, xl , xs, T ,Rrate );
else
end
Vex = y(end-1)+ y(end)*Const.Rex;
js = y(1:n1+n3);
iapp = y(end);
Result = [data.t_s(1),Vex,iapp];

while II < size(data.t_s,1) && Vex > 2.8    
    II=II+1;
    dt = data.t_s(II,1) - data.t_s(II-1,1);
    [ ul, us ] = SVupdate( ul, us, js, dt, T, mesh, Const, MatC );
    xs = transpose(sum(us,1));
    xl = MatC.Vcl*ul;
    if data.Control_mode(II)==1
        Crate = data.C_rate(II);
        [ y ] = BVsolve_CC( mesh, Const, MatC, xl , xs, T ,Crate );
    elseif data.Control_mode(II)==2
        Vapp = data.V_V(II);
        [ y ] = BVsolve_VC( mesh, Const, MatC, xl , xs, T ,Vapp );
    elseif data.Control_mode(II)==3
        Prate = data.C_rate(II);
        [ y ] = BVsolve_PC( mesh, Const, MatC, xl , xs, T ,Prate );
    elseif data.Control_mode(II)==4
        Rrate = data.C_rate(II);
        [ y ] = BVsolve_RC( mesh, Const, MatC, xl , xs, T ,Rrate );
    else
    end
    Vex = y(end-1) + y(end)*Const.Rex;    
    js = y(1:n1+n3);   
    iapp = y(end);
    Result(II,1) = Result(II-1,1) + dt;
    Result(II,2) = Vex;
    Result(II,3) = iapp;
    Result(II,4) = 0;
end
end


%%
function [ y_test ] = Interp_lin( x_test, x, y )
y_test=zeros(size(x_test));
    for ii=1:length(y_test)
        if (x_test(ii)-x(end))*(x_test(ii)-x(1))>=0
            if abs(x_test(ii)-x(end))>abs(x_test(ii)-x(1))
                y_test(ii)=y(1);
            else
                y_test(ii)=y(end);
            end
        else
            nlim=[1,size(x,1)];
            iter=0;
            while ((nlim(2)-nlim(1))>1)&&(iter<size(x,1))
                iter=iter+1;
                idx=floor((nlim(1)+nlim(2))/2);
                if (x_test(ii)-x(idx))*(x_test(ii)-x(1))>=0
                    nlim(1)=idx;
                else
                    nlim(2)=idx;
                end
            end
            y_test(ii)=(y(nlim(2))-y(nlim(1)))/(x(nlim(2))-x(nlim(1)))*(x_test(ii)-x(nlim(1)))+y(nlim(1));
        end
    end
end

%%
function [ y ] = BVsolve_CC( mesh, Const, MatC, xl , xs, T ,Crate)
%%% Basic settings
n1 = mesh.n1;
n2 = mesh.n2;
n3 = mesh.n3;

%%% Variables
iapp = Crate*Const.I_1C/Const.Ac;
xlm = MatC.AVG1*xl;
Kl = MatC.EXT*kappa(xlm*Const.ce0,T);
Kd = Kl./(MatC.EXT*xlm)*2*Const.R*T/Const.Fra*(1+Const.dlnfdce)*(1-Const.tplus);
UU = [ Un(xs(1:n1));Up(xs(n1+1:n1+n3));];
j0 = [ Const.ienref*sqrt(xs(1:n1).*(1-xs(1:n1)).*xlm(1:n1)) ; Const.iepref*sqrt(xs(n1+1:n1+n3).*(1-xs(n1+1:n1+n3)).*xlm(n1+n2+1:n1+n2+n3)); ];
Rf = [ ones(n1,1)*Const.Rfilmn ; ones(n3,1)*Const.Rfilmp ];

%%% Charge conservation equations
%( BA*diag(kl)*B )*phil = BF*js + (BA*diag(kd)*B)*cl ;
MM1 = MatC.BA*diag(Kl)*MatC.B;
MM1(1,:) = 0;
MM1(1,1) = 1;
BF = MatC.BF;
BF(1,:) = 0;
MM2 = (MatC.BA*diag(Kd)*MatC.B);
MM2(1,:) = 0;
% eta = detady + feta
dphildy = sparse([MatC.AVG2*(MM1\BF), ones(n1+n3,1), zeros(n1+n3,1)]);
feta = - MatC.AVG2*(MM1\MM2)*xl - UU;
dphisldy = MatC.dphisdy - dphildy;
detady = [dphisldy,zeros(n1+n3,1)] - [diag(Rf), zeros(n1+n3,3)];

%%% Inilization
Jac = [ [ diag(Rf + Const.R*T./(j0*Const.Fra)) , zeros(n1+n3,3) ] - [dphisldy,zeros(n1+n3,1)]; [ MatC.intn,MatC.intp, 0, 0, 0] ; [ MatC.intn,-MatC.intp, 0, 0, 2];[zeros(1,n1+n3+2),1] ];
Res = [ feta; 0; 0; iapp];
y = Jac\Res;
nr = 1;
iter = 1;

%%% Iteration
while iter<100 && nr>1e-8
    iter =iter+1;
    eta = detady*y + feta;
    Res = [j0.*( exp(0.5*Const.Fra/Const.R/T*eta) - exp(-0.5*Const.Fra/Const.R/T*eta)) - y(1:n1+n3);Jac(n1+n3+1:end,:)*y - [0; 0; iapp]];
    Jac(1:n1+n3,:) = (0.5*Const.Fra/Const.R/T)*diag(j0.*( exp(0.5*Const.Fra/Const.R/T*eta) + exp(-0.5*Const.Fra/Const.R/T*eta)))*detady - [eye(n1+n3),zeros(n1+n3,3)];
    y = y - Jac\Res;
    nr = norm(Res);
end
end

%%
function [ y ] = BVsolve_VC( mesh, Const, MatC, xl , xs, T ,Vapp)
%%% Basic settings
n1 = mesh.n1;
n2 = mesh.n2;
n3 = mesh.n3;

%%% Variables
xlm = MatC.AVG1*xl;
Kl = MatC.EXT*kappa(xlm*Const.ce0,T);
Kd = Kl./(MatC.EXT*xlm)*2*Const.R*T/Const.Fra*(1+Const.dlnfdce)*(1-Const.tplus);
UU = [ Un(xs(1:n1));Up(xs(n1+1:n1+n3));];
j0 = [ Const.ienref*sqrt(xs(1:n1).*(1-xs(1:n1)).*xlm(1:n1)) ; Const.iepref*sqrt(xs(n1+1:n1+n3).*(1-xs(n1+1:n1+n3)).*xlm(n1+n2+1:n1+n2+n3)); ];
Rf = [ ones(n1,1)*Const.Rfilmn ; ones(n3,1)*Const.Rfilmp ];

%%% Charge conservation equations
%( BA*diag(kl)*B )*phil = BF*js + (BA*diag(kd)*B)*cl ;
MM1 = MatC.BA*diag(Kl)*MatC.B;
MM1(1,:) = 0;
MM1(1,1) = 1;
BF = MatC.BF;
BF(1,:) = 0;
MM2 = (MatC.BA*diag(Kd)*MatC.B);
MM2(1,:) = 0;
% eta = detady + feta
dphildy = sparse([MatC.AVG2*(MM1\BF), ones(n1+n3,1), zeros(n1+n3,1)]);
feta = - MatC.AVG2*(MM1\MM2)*xl - UU;
dphisldy = MatC.dphisdy - dphildy;
detady = [dphisldy,zeros(n1+n3,1)] - [diag(Rf), zeros(n1+n3,3)];

%%% Inilization
Jac = [ [ diag(Rf + Const.R*T./(j0*Const.Fra)) , zeros(n1+n3,3) ] - [dphisldy,zeros(n1+n3,1)]; [ MatC.intn,MatC.intp, 0, 0, 0] ; [ MatC.intn,-MatC.intp, 0, 0, 2];[zeros(1,n1+n3+1),1,Const.Rex] ];
Res = [ feta; 0; 0; Vapp];
y = Jac\Res;
nr = 1;
iter = 1;

%%% Iteration
while iter<100 && nr>1e-8
    iter =iter+1;
    eta = detady*y + feta;
    Res = [j0.*( exp(0.5*Const.Fra/Const.R/T*eta) - exp(-0.5*Const.Fra/Const.R/T*eta)) - y(1:n1+n3);Jac(n1+n3+1:end,:)*y - [0; 0; Vapp]];
    Jac(1:n1+n3,:) = (0.5*Const.Fra/Const.R/T)*diag(j0.*( exp(0.5*Const.Fra/Const.R/T*eta) + exp(-0.5*Const.Fra/Const.R/T*eta)))*detady - [eye(n1+n3),zeros(n1+n3,3)];
    y = y - Jac\Res;
    nr = norm(Res);
end
end

%%
function [ y ] = BVsolve_RC( mesh, Const, MatC, xl , xs, T ,Rrate)
%%% Basic settings
n1 = mesh.n1;
n2 = mesh.n2;
n3 = mesh.n3;

%%% Variables
Rload = 3.8/(Rrate*Const.I_1C/Const.Ac);
xlm = MatC.AVG1*xl;
Kl = MatC.EXT*kappa(xlm*Const.ce0,T);
Kd = Kl./(MatC.EXT*xlm)*2*Const.R*T/Const.Fra*(1+Const.dlnfdce)*(1-Const.tplus);
UU = [ Un(xs(1:n1));Up(xs(n1+1:n1+n3));];
j0 = [ Const.ienref*sqrt(xs(1:n1).*(1-xs(1:n1)).*xlm(1:n1)) ; Const.iepref*sqrt(xs(n1+1:n1+n3).*(1-xs(n1+1:n1+n3)).*xlm(n1+n2+1:n1+n2+n3)); ];
Rf = [ ones(n1,1)*Const.Rfilmn ; ones(n3,1)*Const.Rfilmp ];

%%% Charge conservation equations
%( BA*diag(kl)*B )*phil = BF*js + (BA*diag(kd)*B)*cl ;
MM1 = MatC.BA*diag(Kl)*MatC.B;
MM1(1,:) = 0;
MM1(1,1) = 1;
BF = MatC.BF;
BF(1,:) = 0;
MM2 = (MatC.BA*diag(Kd)*MatC.B);
MM2(1,:) = 0;
% eta = detady + feta
dphildy = sparse([MatC.AVG2*(MM1\BF), ones(n1+n3,1), zeros(n1+n3,1)]);
feta = - MatC.AVG2*(MM1\MM2)*xl - UU;
dphisldy = MatC.dphisdy - dphildy;
detady = [dphisldy,zeros(n1+n3,1)] - [diag(Rf), zeros(n1+n3,3)];

%%% Inilization
Jac = [ [ diag(Rf + Const.R*T./(j0*Const.Fra)) , zeros(n1+n3,3) ] - [dphisldy,zeros(n1+n3,1)]; [ MatC.intn,MatC.intp, 0, 0, 0] ; [ MatC.intn,-MatC.intp, 0, 0, 2];[zeros(1,n1+n3+1),1,Const.Rex + Rload] ];
Res = [ feta; 0; 0; 0];
y = Jac\Res;
nr = 1;
iter = 1;

%%% Iteration
while iter<100 && nr>1e-8
    iter =iter+1;
    eta = detady*y + feta;
    Res = [j0.*( exp(0.5*Const.Fra/Const.R/T*eta) - exp(-0.5*Const.Fra/Const.R/T*eta)) - y(1:n1+n3);Jac(n1+n3+1:end,:)*y - [0; 0; 0]];
    Jac(1:n1+n3,:) = (0.5*Const.Fra/Const.R/T)*diag(j0.*( exp(0.5*Const.Fra/Const.R/T*eta) + exp(-0.5*Const.Fra/Const.R/T*eta)))*detady - [eye(n1+n3),zeros(n1+n3,3)];
    y = y - Jac\Res;
    nr = norm(Res);
end
end

%%
function [ y , iapp] = BVsolve_PC( mesh, Const, MatC, xl , xs, T ,Prate)
%%% Basic settings
n1 = mesh.n1;
n2 = mesh.n2;
n3 = mesh.n3;

%%% Variables
Papp = Prate*Const.Pnom;
xlm = MatC.AVG1*xl;
Kl = MatC.EXT*kappa(xlm*Const.ce0,T);
Kd = Kl./(MatC.EXT*xlm)*2*Const.R*T/Const.Fra*(1+Const.dlnfdce)*(1-Const.tplus);
UU = [ Un(xs(1:n1));Up(xs(n1+1:n1+n3));];
j0 = [ Const.ienref*sqrt(xs(1:n1).*(1-xs(1:n1)).*xlm(1:n1)) ; Const.iepref*sqrt(xs(n1+1:n1+n3).*(1-xs(n1+1:n1+n3)).*xlm(n1+n2+1:n1+n2+n3)); ];
Rf = [ ones(n1,1)*Const.Rfilmn ; ones(n3,1)*Const.Rfilmp ];

%%% Charge conservation equations
%( BA*diag(kl)*B )*phil = BF*js + (BA*diag(kd)*B)*cl ;
MM1 = MatC.BA*diag(Kl)*MatC.B;
MM1(1,:) = 0;
MM1(1,1) = 1;
BF = MatC.BF;
BF(1,:) = 0;
MM2 = (MatC.BA*diag(Kd)*MatC.B);
MM2(1,:) = 0;
% eta = detady + feta
dphildy = sparse([MatC.AVG2*(MM1\BF), ones(n1+n3,1), zeros(n1+n3,1)]);
feta = - MatC.AVG2*(MM1\MM2)*xl - UU;
dphisldy = MatC.dphisdy - dphildy;
detady = [dphisldy,zeros(n1+n3,1)] - [diag(Rf), zeros(n1+n3,3)];

%%% Inilization
Jac = [ [ diag(Rf + Const.R*T./(j0*Const.Fra)) , zeros(n1+n3,3) ] - [dphisldy,zeros(n1+n3,1)]; [ MatC.intn,MatC.intp, 0, 0, 0] ; [ MatC.intn,-MatC.intp, 0, 0, 2]; ];
Res = [ feta; 0; 0;];
bb = Jac(:,1:n1+n3+2)\Res;
kk = -Jac(:,1:n1+n3+2)\Jac(:,end);
ub = bb(end);
Zb = kk(end);
upp = ub^2+4*(Zb+Const.Rex)*Papp;
iapp = (-ub + sqrt(max(upp,0)))/(2*(Zb+Const.Rex));
y = [bb + kk*iapp;iapp];
nr = 1;
iter = 1;

%%% Iteration
while iter<100 && nr>1e-8
    iter =iter+1;
    eta = detady*y + feta;
    Res = [j0.*( exp(0.5*Const.Fra/Const.R/T*eta) - exp(-0.5*Const.Fra/Const.R/T*eta)) - y(1:n1+n3);Jac(n1+n3+1:n1+n3+2,:)*y; (y(n1+n3+2)+y(end)*Const.Rex)*y(end)-Papp ];
    Jac(1:n1+n3,:) = (0.5*Const.Fra/Const.R/T)*diag(j0.*( exp(0.5*Const.Fra/Const.R/T*eta) + exp(-0.5*Const.Fra/Const.R/T*eta)))*detady - [eye(n1+n3),zeros(n1+n3,3)];
    Jac(n1+n3+3,:) =0;
    Jac(end,n1+n3+2) = y(end);
    Jac(end,end) = y(n1+n3+2)+2*y(end)*Const.Rex;
    y = y - Jac\Res;
    nr = norm(Res);
end
end

%%
function [ ul, us ] = SVupdate( ul, us, js, dt, T, mesh, Const, MatC )
%%% Basic settings
n1 = mesh.n1;
n3 = mesh.n3;
%%% Electrolyte diffusion
fcl = MatC.Fcl*js;
for ii = 1:size(MatC.Lcl,1)
    lambda = MatC.Lcl(ii)*exp(-Const.Eadl/Const.R*(1/T-1/Const.Tref));
    if ii ~= MatC.zidx
        ul(ii) = ul(ii)*exp(-lambda*dt) + (1 - exp(-lambda*dt))/lambda*fcl(ii);
    else
        ul(ii) = ul(ii) + fcl(ii)*dt;
    end
end
%%% Anode diffusion
lambda = MatC.Lcn*exp(-Const.Eadn/Const.R*(1/T-1/Const.Tref));
for jj = 1:n1    
    for ii = 1:size(lambda,2)-1
        us(ii,jj) = us(ii,jj)*exp(-lambda(ii)*dt) + (1 - exp(-lambda(ii)*dt))/lambda(ii)*MatC.Gcn(ii)*js(jj);
    end 
    us(ii+1,jj) = us(ii+1,jj) + dt*MatC.Gcn(ii+1)*js(jj);
end
%%% Cathode diffusion
lambda = MatC.Lcp*exp(-Const.Eadp/Const.R*(1/T-1/Const.Tref));
for jj = n1+1:n1+n3    
    for ii = 1:size(lambda,2)-1
        us(ii,jj) = us(ii,jj)*exp(-lambda(ii)*dt) + (1 - exp(-lambda(ii)*dt))/lambda(ii)*MatC.Gcp(ii)*js(jj);
    end 
    us(ii+1,jj) = us(ii+1,jj) + dt*MatC.Gcp(ii+1)*js(jj);
end
end

%%
function F=Un(xn)%Open circuit potential of negative electrode [V]
b=[0.164,0.48,90,0.01,0.025,300,0.136,0.032,100,0.161,0.032,50,0.22,0.032,80,0.505,0.1,90,0.965];
Un_const=b(1);
Un_exp1=b(2)*exp(-b(3)*(xn-b(4)));
Un_atan1=-b(5)/pi*atan(b(6)*(xn-b(7)));
Un_atan2=-b(8)/pi*atan(b(9)*(xn-b(10)));
Un_atan3=-b(11)/pi*atan(b(12)*(xn-b(13)));
Un_atan4=-b(14)/pi*atan(b(15)*(xn-b(16)));
Un_exp2=-b(17)*exp(b(18)*(xn-b(19)));
F=Un_const+Un_exp1+Un_atan1+Un_atan2+Un_atan3+Un_atan4+Un_exp2;
end

%%
function F=Up(xp)%Open circuit potential of positive electrode [V]
b=[3.67287671914126,0.0392134821245406,0.999470977852800,0.419937084206002,12.8799319870686,-47.2135901049774,64.6252944547241,-40.4701303771546,10.1693381074356;];
b_Nernst=b(1:3);
OCP_charge_Nernst=b_Nernst(1)+b_Nernst(2)*log(b_Nernst(3)-xp);
b_poly=b(4:end);
order=length(b_poly);
dOCP_charge_fit=0;
for i=1:order
    dOCP_charge_fit=dOCP_charge_fit+b_poly(order+1-i)*xp.^(i);
end
F=OCP_charge_Nernst+dOCP_charge_fit;
end

%%
function F=kappa(ce,T)%Electrolyte conductivity (S/m)
F=(3.45*exp(-798/T)*(ce/1000).^3-48.5*exp(-1080/T)*(ce/1000).^2+244*exp(-1440/T)*(ce/1000));
end

