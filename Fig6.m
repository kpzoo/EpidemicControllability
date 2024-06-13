% Reproduces Figure 6
clearvars; close all; clc;

% Assumptions and notes
% - assess several margins and measures of stability
% - targeted control focuses on 1 of the 2 types
% - different types have different generation times
% - superspreading and fast-slow transmission

% Figure defaults
set(groot, 'defaultAxesTickLabelInterpreter', 'latex', 'defaultLegendInterpreter', 'latex');
fnt = 24; grey1 = 0.5*ones(1, 3); grey2 = 0.8*ones(1, 3);
set(0, 'defaultTextInterpreter', 'latex', 'defaultAxesFontSize', fnt);

% Complex s and times for any simulations
s = tf('s'); dt = 0.02; t = 0:dt:100; lent = length(t);
% Step and impulse input across times
ustep = 10*ones(1, lent); uimp = [10*ones(1, 50), zeros(1, lent-50)];

%% Component generation time distributions

% Mean generation times of both types
lg = 5; g1 = linspace(2, 4, lg); g2 = 8*ones(1, lg);
% Shape and scale parameters of types
g1pms = zeros(2, lg); g1pms(1, :) = 3*ones(1, lg);
g2pms = g1pms; g2pms(1, :) = 9*ones(1, lg);
g1pms(2, :) = g1./g1pms(1, :); g2pms(2, :) = g2./g2pms(1, :);

% Distributions and Laplace transforms for all components
W1 = s*ones(1, lg); w1 = zeros(lg, lent); W2 = W1; 
w2 = w1; wstat1 = cell(1, lg); wstat2 = wstat1;
for ii = 1:lg
    % Generation time distribution formula with null weight
    [wtemp, Wtemp, wsttemp] = generationGammaWeighted(g1pms(:, ii),...
        g2pms(:, ii), 0.5, s, t, lent);
    % Extract component distributions
    w1(ii, :) = wtemp(1, :); W1(ii) = Wtemp(1); wstat1{ii} = wsttemp{1};
    w2(ii, :) = wtemp(2, :); W2(ii) = Wtemp(2); wstat2{ii} = wsttemp{2};
end

%% Control effort if controlling superspreading (fixed W)

% Define a range of superspreading proportions
leps = 25; eps1 = linspace(0.05, 0.5, leps); eps2 = 1 - eps1;
% Define R2, range of R1 for superspreading and critical eps1
R2 = 1.1; R1 = 1.1:0.2:5.5; l1 = length(R1); ecrit = 1 - 1/R2;

% Critical constant targeted and non-selective controllers
K1crit = zeros(leps, l1); Kcrit = K1crit; 
% Checks of gain and delay margins
checkmarg = K1crit; checkmarg1 = K1crit;

% Compute expected critical controllers and check margins are as expected
for jj = 1:leps
    for ii = 1:l1
        % Controllers expected to obtain gain margin of 1
        Kcrit(jj, ii) = 1/(eps1(jj)*R1(ii) + eps2(jj)*R2);
        K1crit(jj, ii) = (1 - eps2(jj)*R2)/(eps1(jj)*R1(ii));

        % Check controllers for a number of generation times
        for gg = 1:lg
            % Loop TFs for targeted and non-selective control
            L = -(eps1(jj)*R1(ii)*W1(gg) + eps2(jj)*R2*W2(gg))*Kcrit(jj, ii);
            L1 = -(eps1(jj)*R1(ii)*W1(gg)*K1crit(jj, ii) + eps2(jj)*R2*W2(gg));

            % Margins (minimum) for non-selective
            marg = allmargin(L);
            if isempty(marg.GainMargin)
                marg.GainMargin = inf;
            end
            [gmarg, idwpc] = min(marg.GainMargin);
            wpc = marg.GMFrequency(idwpc);

            if max(abs(gmarg - 1)) < 10^-9 && wpc == 0
                checkmarg(jj, ii) = 1;
            else
                assignin('base', 'gmarg', gmarg);
                assignin('base', 'wpc', wpc);
                error('Margins not critical');
            end

            % Margins (minimum) for targeted
            marg1 = allmargin(L1);
            if isempty(marg1.GainMargin)
                marg1.GainMargin = inf;
            end
            [gmarg1, idwpc1] = min(marg1.GainMargin);
            wpc1 = marg1.GMFrequency(idwpc1);

            if max(abs(gmarg1-1)) < 10^-9 && wpc1 == 0
                checkmarg1(jj, ii) = 1;
            else
                assignin('base', 'gmarg1', gmarg1);
                assignin('base', 'wpc1', wpc1);
                error('Margins not critical');
            end
        end
    end
end


%% Control effort if fast and slow transmission and targetting either

% Assume eps*R is fixed across types to A and define controller form
A = 0.7; Kdel = (1/7)*(1+20*s)/((1+s)^2);
%Kdel = (1/7)*(1+20*s)/((1+2*s)^2);
% Ensure controller has correct integral
t = 0.0001:0.0001:200;  syms('s'); 
Kx = ilaplace((1/7)*(1+20*s)/((1+2*s)^2)); 
Kint = eval(Kx); int = trapz(t, Kint); 
disp(['Integral of k(t): ' num2str(int)]);

% Reset variable types
s = tf('s'); t = 0:dt:100;

% Find maximum magnitude of controller
[magtemp, ~, wtemp] = bode(Kdel); [magK, idw] = max(magtemp);
disp(['Max |K(s)| are: ' num2str(magK)]);

% Store TFs and poles for targeted controllers
L = s*ones(lg, 2); G = L; p = cell(lg, 2); 
% Margins and poles and crossover frequencies
gmarg = zeros(lg, 2); dmarg = gmarg; wpc = gmarg; 

% For every R1 and W compute with controllers
for ii = 1:lg
    % Compute TFs and margins for both targeted controllers
    L(ii, 1) = -(Kdel*W1(ii) + W2(ii))*A; G(ii, 1) = feedback(1, L(ii, 1));
    L(ii, 2) = -(W1(ii) + Kdel*W2(ii))*A; G(ii, 2) = feedback(1, L(ii, 2));

    % All margins of TFs
    marg1 = allmargin(L(ii, 1));
    if isempty(marg1.DelayMargin)
        marg1.DelayMargin = inf;
    end
    if isempty(marg1.GainMargin)
        marg1.GainMargin = inf;
    end
    marg2 = allmargin(L(ii, 2));
    if isempty(marg2.DelayMargin)
        marg2.DelayMargin = inf;
    end
    if isempty(marg1.GainMargin)
        marg2.GainMargin = inf;
    end

    % MG, MD and crossover frequencies
    dmarg(ii, 1) = min(marg1.DelayMargin);
    dmarg(ii, 2) = min(marg2.DelayMargin);
    [gmarg(ii, 1), idwpc1] = min(marg1.GainMargin);
    wpc(ii, 1) = marg1.GMFrequency(idwpc1);
    [gmarg(ii, 2), idwpc2] = min(marg2.GainMargin);
    wpc(ii, 2) = marg2.GMFrequency(idwpc2);

    % Poles from both control strategies
    p{ii, 1} = esort(pole(G(ii, 1))); p{ii, 2} = esort(pole(G(ii, 2)));
end


% Examine step responses for control types
istep1 = zeros(lg, lent); istep2 = istep1; 
for ii = 1:lg
    % Linear systems simulations
    istep1(ii, :) = lsim(G(ii, 1), ustep, t);
    istep2(ii, :) = lsim(G(ii, 2), ustep, t);
end

%% Fig 6 analysis of superspreading and variants

figure('Position', [10 10 1000 600]);
cols = {'g', 'b', grey1, grey2, 'r'};
% Critical control for given proportions of types
subplot(1, 2, 1); hold on;
plot(eps1, K1crit(:, 2:end-1), '-', 'LineWidth', 2, 'Color', grey2);
plot(eps1, Kcrit(:, 2:end-1), '--', 'LineWidth', 2, 'Color', grey2);
plot(eps1, K1crit(:, 1), '-', 'LineWidth', 2, 'Color', 'b');
plot(eps1, Kcrit(:, 1), '--', 'LineWidth', 2, 'Color', 'b');
plot(eps1, K1crit(:, end), '-', 'LineWidth', 2, 'Color', 'r');
plot(eps1, Kcrit(:, end), '--', 'LineWidth', 2, 'Color', 'r');
plot([ecrit ecrit], [0 0.25], 'k--', 'LineWidth', 2);
plot([0.2 0.2], [0 1], '--', 'LineWidth', 2, 'Color', grey1);
ylabel('$K^*, K_1 ^*  \,|\, M_G = 1$', 'FontSize', fnt);
hold off; box off; grid off; 
ylim([0 1.05]); xlim([0 0.5]);
xlabel('$\epsilon_1 | R_2 = 1.1$', 'FontSize', fnt);

% Dynamical effects of control of variants of different speed
subplot(1, 2, 2); hold on;
for ii = 1:lg
    plot(t, istep1(ii, :), '--', 'LineWidth', 2, 'Color', cols{ii});
    plot(t, istep2(ii, :), '-', 'LineWidth', 2, 'Color', cols{ii});
end
hold off; box off; grid off;
ylabel('$i(t) | \alpha = 0.7$', 'FontSize', fnt); 
xlabel('$t$ (days)', 'FontSize', fnt); 




