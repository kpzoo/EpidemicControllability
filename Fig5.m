% Controllability margins under observation noise
clearvars; close all; clc;

% Assumptions and notes
% - assess several margins and measures of stability
% - consider combinations of reporting rates and delays

% Figure defaults
set(groot, 'defaultAxesTickLabelInterpreter', 'latex', 'defaultLegendInterpreter', 'latex');
fnt = 24; grey1 = 0.5*ones(1, 3); grey2 = 0.8*ones(1, 3);
set(0, 'defaultTextInterpreter', 'latex', 'defaultAxesFontSize', fnt);

% Complex s and times for any simulations
s = tf('s'); dt = 0.02; t = 0:dt:100; lent = length(t);
% Step and impulse input across times
ustep = 10*ones(1, lent); uimp = [10*ones(1, 50), zeros(1, lent-50)];

%% Statistics for fixed mean generation times

% Generation types [det exp gam gam bimod] with mean g0
g0 = 6.5; GTtype = [1 3 2 2 4]; lg = length(GTtype);
% Generation time shapes and scales (1 extra shape/scale for bimodal)
gshapes = [0 1 3 8 9 30]; gscales = [0 g0 g0/3 g0/8 1/3 1/3];

% Compute all distributions and Laplace transforms
W = s*ones(1, lg); w = zeros(lg, lent); wstat = cell(1, lg);
for ii = 1:lg
    % Assign parameters
    GT.mean = g0; GT.scale = gscales(ii); GT.shape = gshapes(ii);
    if ii == lg
        % Bimodal distribution case
        GT.scale = gscales(lg:lg+1); GT.shape = gshapes(lg:lg+1);
    end
    % Generation time properties in t and s domain
    [w(ii, :), W(ii), wstat{ii}] = generationLaplace(GT, GTtype(ii), s, t);
end


%% Control effort required if only under-reporting

% R0 range to examine
Rset = 1:0.1:5; l = length(Rset);
% Control effort without noise
Kideal = 1./Rset; ctrlType = 2;

% Control required to get critical epidemic
K = cell(1, l); rho = K;
for ii = 1:l
    % Define range of rho can plot
    if Rset(ii) >= 1
        rho{ii} = linspace(1-1/Rset(ii), 1, 50);
    else
        rho{ii} = linspace(0.1, 1, 50);
    end
    % Obtain control effort on cases required for criticality
    K{ii} = (Kideal(ii) - 1)./rho{ii} + 1;
end

%% Control effort required if only reporting delays

% Various delay magnitudes
tau = linspace(0, 5, 50); ltau = length(tau);

% Define the delay as exponential and gamma
He = s*ones(1, ltau); Hg = He;
for ii = 1:ltau
    He(ii) = 1/(tau(ii)*s + 1);
    Hg(ii) = 1/(tau(ii)*s/4 + 1)^(4);
end

% TFs for every delay choice
L = s*ones(lg, ltau); G = L; Lg = L; Gg = G;
p = cell(lg, ltau); pg = p; 
% Margins and crossover frequencies
gmarg = zeros(lg, ltau); dmarg = gmarg; gmargg = gmarg; 
dmargg = gmarg; wpc = gmarg; wpcg = gmarg;

% Pick an R and run over W and tau
for jj = 1:lg
    % Pick a W(s) and an R
    Wtau = W(jj); R = 2;

    % Compute TFs and margins
    for ii = 1:ltau
        % Open and closed loop TFs for exp
        L(jj, ii) = -R*Wtau*(1 + He(ii)*(0.5/R-1));
        G(jj, ii) = 1/(1 + L(jj, ii));

        % Open and closed loop TFs for gamma
        Lg(jj, ii) = -R*Wtau*(1 + Hg(ii)*(0.5/R-1));
        Gg(jj, ii) = 1/(1 + Lg(jj, ii));

        % Margins (minimum) for exp
        marg = allmargin(L(jj, ii));
        if isempty(marg.DelayMargin)
            marg.DelayMargin = inf;
        end
        if isempty(marg.GainMargin)
            marg.GainMargin = inf;
        end
        dmarg(jj, ii) = min(marg.DelayMargin);
        [gmarg(jj, ii), idwpc] = min(marg.GainMargin);
        wpc(jj, ii) = marg.GMFrequency(idwpc);

        % Margins (minimum) for gamma
        marg = allmargin(Lg(jj, ii));
        if isempty(marg.DelayMargin)
            marg.DelayMargin = inf;
        end
        if isempty(marg.GainMargin)
            marg.GainMargin = inf;
        end
        dmargg(jj, ii) = min(marg.DelayMargin);
        [gmargg(jj, ii), idwpc] = min(marg.GainMargin);
        wpcg(jj, ii) = marg.GMFrequency(idwpc);

        % Poles of controlled system
        p{jj, ii} = esort(pole(G(jj, ii)));
        pg{jj, ii} = esort(pole(Gg(jj, ii)));
    end
end

% For a mean delay of 3 examine step response
idtau = find(tau <= 3, 1, 'last'); tau0 = tau(idtau);
isteptau = zeros(lg, lent); isteptaug = isteptau; istep0 = isteptau;
for ii = 2:lg
    % Linear systems simulations
    isteptau(ii, :) = lsim(G(ii, idtau), ustep, t);
    isteptaug(ii, :) = lsim(Gg(ii, idtau), ustep, t);
    % Also simulation case with no delay
    istep0(ii, :) = lsim(G(ii, 1), ustep, t);
end

%% Across range of R0 obtain margins under specific delay and rho

% R0 range to examine
Re = 1:0.05:3; l = length(Re);

% Margins for direct control type
Le = s*ones(lg, l); Ge = Le; Gh = Le; GCe = Le;
% Define observation noise as fraction and delay
rhoe = 0.8; He = 1/(3*s+1);

% Poles and margins of systems
marge = cell(lg, l); pe = marge; pmaxe = marge;
gmarge = zeros(lg, l); dmarge = gmarge; 
% Compare when no He effect
margh = cell(lg, l); ph = margh; pmaxh = margh;
gmargh = zeros(lg, l); dmargh = gmargh; 

% Obtain TFs and their properties
for ii = 1:l
    for jj = 1:lg
        % Main function constructing open and closed loops
        [~, pe{jj, ii}, Le(jj, ii), Ge(jj, ii), marge{jj, ii}, ...
            pmaxe{jj, ii}] = getOLCLcontrolNoise(Re(ii), W(jj), 1/(3*Re(ii)), rhoe, He, ctrlType);
        % Store minimum gain and delay margins
        gmarge(jj, ii) = min(marge{jj, ii}.g); 
        dmarge(jj, ii) = min(marge{jj, ii}.d); 

        % Repeat but no delay now
        [~, ph{jj, ii}, ~, Gh(jj, ii), margh{jj, ii}, pmaxh{jj, ii}] = ...
            getOLCLcontrolNoise(Re(ii), W(jj), 1/(3*Re(ii)), rhoe, 1, ctrlType);
        % Store minimum gain and delay margins
        gmargh(jj, ii) = min(margh{jj, ii}.g);
        dmargh(jj, ii) = min(margh{jj, ii}.d);

    end
end

% Dominant poles and check worst perturbations
pmaxe = cell2mat(pmaxe')'; 

% Obtain TF of cases versus imports
for ii = 1:l
    for jj = 1:lg
        GCe(jj, ii) = Ge(jj, ii)*rhoe*He;
    end
end

% Take an example and get step and impulse responses
ide = find(Rset <= 4, 1, 'last'); Ride = Re(ide);
istepe = zeros(lg, lent); isteph = istepe; cstepe = istepe;
for ii = 1:lg
    % Linear systems simulations
    istepe(ii, :) = lsim(Ge(ii, ide), ustep, t);
    cstepe(ii, :) = lsim(GCe(ii, ide), ustep, t);
    isteph(ii, :) = lsim(Gh(ii, ide), ustep, t);
end

% Check bounds on controllability as a polar plot
TFright = 1 - 1/(2*W(3));

%% Fig 4 analysis of controller with noise

figure('Renderer', 'painters', 'Position', [10 10 1000 1000]);
cols = {'g', 'b', grey1, grey2, 'r'};
% Critical control for a given reporting rate
subplot(2, 2, 1); hold on;
for ii = 3:l-1
    plot(rho{ii}, K{ii}, 'LineWidth', 2, 'Color', grey2);
end
plot(rho{2}, K{2}, 'LineWidth', 2, 'Color', 'b');
plot(rho{end}, K{end}, 'LineWidth', 2, 'Color', 'r');
plot(rho{1}, K{1}, 'k--', 'LineWidth', 2);
hold off; box off; grid off;
xlabel('$\rho$', 'FontSize', fnt); 
xlim([0.1, 1]); ylim([0 1.05]);
ylabel('$K^*$', 'FontSize', fnt);

% Margins with increasing delay
subplot(2, 2, 2); hold on;
for ii = 2:lg
    plot(tau, gmarg(ii, :), '--', 'Color', cols{ii}, 'LineWidth', 2);
    plot(tau, gmargg(ii, :), '-', 'Color', cols{ii}, 'LineWidth', 2);
end
plot(tau, 2*ones(size(tau)), 'k--', 'LineWidth', 2);
hold off; box off; grid off;
xlabel('$\tau$ (days)', 'FontSize', fnt); 
ylabel('$M_G$', 'FontSize', fnt); 
xlim([min(tau), max(tau)]); ylim([0.5 2.05]);

% Additional phase crossover
axes('Position',[0.615 0.61 0.15 0.15]);
hold on;
for ii = 2:lg
    plot(tau, wpc(ii, :), '--', 'Color', cols{ii}, 'LineWidth', 2);
    plot(tau, wpcg(ii, :), '-', 'Color', cols{ii}, 'LineWidth', 2);
end
hold off; box off; grid off;
h = gca; h(1).YColor = 'k';
title('$\omega_{PC}$', 'FontSize', fnt); 
xlabel('$\tau$', 'FontSize', fnt); 
xlim([min(tau), max(tau)]);

% Step response to delays in surveillance
subplot(2, 2, 3); 
hold on;
for ii = 2:lg
    plot(t, istep0(ii, :), '-.', 'Color', cols{ii}, 'LineWidth', 2);
    plot(t, isteptau(ii, :), '--', 'Color', cols{ii}, 'LineWidth', 2);
    plot(t, isteptaug(ii, :), '-', 'Color', cols{ii}, 'LineWidth', 2);
end
hold off; box off; grid off;
ylim([10 30]);
xlabel('$t$ (days)', 'FontSize', fnt); 
ylabel('$i(t)$', 'FontSize', fnt); 

% Additional delay margin
axes('Position',[0.287,0.137,0.15,0.15]); 
hold on;
for ii = 2:lg
    semilogy(tau, dmarg(ii, :), '--', 'Color', cols{ii}, 'LineWidth', 2);
    semilogy(tau, dmargg(ii, :), '-', 'Color', cols{ii}, 'LineWidth', 2);
end
hold off; box off; grid off;
xlabel('$\tau$', 'FontSize', fnt); 
title('$M_D$', 'FontSize', fnt); 
xlim([2, max(tau)]); ylim([0 30]);

% Step response to specific delay and rho
subplot(2, 2, 4); 
hold on;
for ii = 2:lg
    plot(t, istepe(ii, :), '-', 'Color', cols{ii}, 'LineWidth', 2);
    plot(t, cstepe(ii, :), '--', 'Color', cols{ii}, 'LineWidth', 2);
end
hold off; box off; grid off;
xlabel('$t$ (days)', 'FontSize', fnt); 
ylabel('$i(t)$', 'FontSize', fnt); 

% Additional delay margin
axes('Position',[0.587,0.137,0.15,0.15]); 
hold on;
for ii = 2:lg
    semilogy(Rset, dmarge(ii, :), '-', 'Color', cols{ii}, 'LineWidth', 2);
end
hold off; box off; grid off;
xlabel('$R$', 'FontSize', fnt); 
title('$M_D$', 'FontSize', fnt); 
xlim([min(Rset), max(Rset)]); ylim([0 30]);

