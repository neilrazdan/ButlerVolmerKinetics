clear all; close all; clc
tic
%% Rate Constants
k10 = 1; k10rev = 10^1; 
k20 = 3; k20rev = 3*10^-1;

% Equilibrium Constants
K10 = k10/k10rev;
% Activities
aOH = 1; aH2O = 1; aH2 = 1;
% Symmetry Factors
B = 0.5;
% Faraday
f = 96500/8.314/298;
F = 96500;
RT = 8.314*298;
%% Arrays
eta = linspace(-0.03,0.03,50);
E = eta;
%% Rate Expression
r77 = K10.*k20.*exp(f.*eta).*aOH^2.*exp(B.*f.*eta)./(1 + K10.*exp(f.*eta).*aOH);
theta0 = 1./(1 + K10.*exp(f.*eta).*aOH);
%% DoRC, steady state parameters
delta = 10^-4;
eps = 10^-20;
tspan = [0 10^20];
%% Arrays for results
S = zeros(length(E),4)'; r = zeros(1,length(E)); A = r; O = r; sigma_app = r; TRC_A = r; r_trc = r;
%% Solvers and Nested Loops
for w = 1:length(E)
    disp(w)
    for m = 1:1
        k0 = [k10 k10rev k20 k20rev];
        for j = 1:length(k0)
            K = zeros(2,length(k0));
            K(1,:) = k0; K(2,:) = k0;
            K(2,j) = k0(j)*(1+delta);
            for i = 1:2
                %% Parameters
                k1 = K(i,1).*exp(E(w).*B.*F./RT);
                k1rev = K(i,2).*exp(E(w).*(B-1).*F./RT);
                k2 = K(i,3).*exp(E(w).*B.*F./RT);
                k2rev = K(i,4).*exp(E(w).*(B-1).*F./RT);
                
                %% Sites
                a = @(x) x(1); 
                o = @(x) 1 - a(x);
                
                % Step 1
                r1a = @(x) k1.*aOH.*o(x).*aH2 - k1rev.*a(x);
                
                % Step 2
                r2a = @(x) -k2.*a(x).*aOH + k2rev.*o(x).*aH2O;
                Fun = @(t,x) [r1a(x) + r2a(x)];
                
                C0 = [eps];
                
                [time,C] = ode23s(Fun,tspan,C0);
                t = time;
                a = C(:,1);
                o = 1 - a;
                
                em = round(length(a)*0.20);
                z1(i,j) = (k1rev.*mean(a((end-em):end)))./(k1.*aOH.*aH2.*mean(o((end-em):end)));
                z2(i,j) = (k2rev.*mean(o((end-em):end)).*aH2O)./(k2.*aOH.*mean(a((end-em):end)));
                Z_eff(i,j) = z1(i,j).*z2(i,j);
                R(i,j) = k2.*mean(a((end-em):end)).*aOH - k2rev.*mean(o((end-em):end)).*aH2O;
                R_fwd(i,j) = R(i,j)/(1-Z_eff(i,j));
                R_rev(i,j) = Z_eff(i,j).*R(i,j)/(1-Z_eff(i,j));
                
            end
            S(j,w) = (R(2,j) - R(1,j))./(delta*R(1,j));
            S_fwd(j,w) = (R_fwd(2,j) - R_fwd(1,j))./(delta*R_fwd(1,j));
            S_rev(j,w) = (R_rev(2,j) - R_rev(1,j))./(delta*R_rev(1,j));
            r(m,w) = R(1,1);
            A(m,w) = mean(a((end-em):end));
            O(m,w) = mean(o((end-em):end));
            Z1(m,w) = z1(1,1);
            Z2(m,w) = z2(1,1);
        end
    end
end
%% Extract results
XRC1 = S(1,:) + S(2,:); XRC2 = S(3,:) + S(4,:);
XRC1_fwd = S_fwd(1,:) + S_fwd(2,:); XRC2_fwd = S_fwd(3,:) + S_fwd(4,:);
XRC1_rev = S_rev(1,:) + S_rev(2,:); XRC2_rev = S_rev(3,:) + S_rev(4,:);
zeff = Z1.*Z2;

% Transfer coefficient
alpha = S(1,:).*B + S(2,:).*(B - 1) + S(3,:).*B + S(4,:).*(B - 1);
tafel = 1000*log(10)/f./alpha;

% Requisite ratio for QE step 1
ratio = k10rev.*exp(E.*(B-1).*F./RT)./(k10.*exp(E.*B.*F./RT));

% Unidirectional transfer coefficient
alpha_a = alpha - 2*zeff./(1-zeff);
alpha_c = alpha - 2./(1-zeff);

% Rxn orders
PhiH2fwd = S_fwd(1,:);
PhiOHfwd = S_fwd(1,:) + S_fwd(3,:);

%% Rates
figure
set(gca,'fontweight','bold','fontsize',11,'box','on');
set(gca,'linewidth', 2,'fontsize',12,'fontname','Arial')
hold on
plot(E*1000,alpha_a,'b-','linewidth',2);
plot(E*1000,5*A,'k--','linewidth',2);
plot(E*1000,S_fwd(1,:),'r.','linewidth',2);
xlim([-30 30])
ylim([0 2])
legend('\alpha_{a} = ... per Eq. ?','5 \times \theta_{H^{*}}','\Psi_{H_{2}}')
legend boxoff
%% Rates
figure
set(gca,'fontweight','bold','fontsize',11,'box','on');
set(gca,'linewidth', 2,'fontsize',12,'fontname','Arial')
hold on
plot(E*1000,XRC1,'k--','linewidth',2.5);
plot(E*1000,2*PhiH2fwd + A - alpha_a - B,'k-','linewidth',1);
xlim([-30 30])
ylim([0 0.5])
legend('X_{RC,1}','Eq. ?')
legend boxoff

%%
figure
set(gca,'fontweight','bold','fontsize',11,'box','on');
set(gca,'linewidth', 2,'fontsize',12,'fontname','Arial')
hold on
plot(E*1000,XRC1,'b--','linewidth',1);
plot(E*1000,XRC2,'r--','linewidth',1);

%% Coverages
figure
set(gca,'fontweight','bold','fontsize',11,'box','on');
set(gca,'linewidth', 2,'fontsize',12,'fontname','Arial')
hold on
plot(E,A,'b--','linewidth',2.5);
plot(E,S_fwd(2,:) + S_fwd(3,:),'b-','linewidth',1);
