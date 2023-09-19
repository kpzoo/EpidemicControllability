% Margins from different generation time shapes
clearvars; close all; clc;

% Assumptions and notes
% - changing generation time structure and r or R
% - assess several margins and measures of stability
% - define margins by their min value if multiple

% Figure defaults
set(groot, 'defaultAxesTickLabelInterpreter', 'latex', 'defaultLegendInterpreter', 'latex');
fnt = 24; grey1 = 0.5*ones(1, 3); grey2 = 0.8*ones(1, 3);
set(0, 'defaultTextInterpreter', 'latex', 'defaultAxesFontSize', fnt);

% Complex s and times for any simulations
s = tf('s'); dt = 0.02; t = 0:dt:50; lent = length(t);
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

%% Across range of R0 at same g0 get margins with K = 1

% R0 range to examine
Rset = 0.2:0.1:5; l = length(Rset); 
% The last R0 <= 1
id1 = find(Rset <= 1, 1, 'last'); Rid1 = Rset(id1);

% Margins for direct control type
ctrlType = 1; L = s*ones(lg, l); G = L; 

% Zeros, poles and margins of systems
marg = cell(lg, l); z = marg; p = marg; pmax = marg;
% Reduced form poles (balanced truncation) 
pred = zeros(lg, l); idwarn = pred;
% Store margins from first case (including diskmargin)
gmarg = pred; dmarg = pred; disk = pred; dworst = marg;

% Obtain TFs and their properties
for ii = 1:l
    for jj = 1:lg
        % Main function constructing open and closed loops
        [z{jj, ii}, p{jj, ii}, L(jj, ii), G(jj, ii), marg{jj, ii}, ...
            pmax{jj, ii}] = getOLCLcontrolNoise(Rset(ii), W(jj), 1, 1, 1, ctrlType);
        % Store minimum gain and delay margins
        gmarg(jj, ii) = min(marg{jj, ii}.g); 
        dmarg(jj, ii) = min(marg{jj, ii}.d);
        disk(jj, ii) = marg{jj, ii}.disk(1);
        dworst{jj, ii} = marg{jj, ii}.worst;

        % Get poles of reduced functions (ensure no static gains)
        try
            if order(balred(G(jj, ii), 1)) > 0
                pred(jj, ii) = pole(balred(G(jj, ii), 1));
            else
                pred(jj, ii) = pole(balred(G(jj, ii), 2));
            end
        catch
            % Note issue and store in idwarn
            pred(jj, ii) = max(real(pole(G(jj, ii))));
            idwarn(jj, ii) = 1;
        end  
    end
end

% Dominant poles and check worst perturbations
pmax = cell2mat(pmax')'; 
try
    dworst = cell2mat(dworst);
catch
    % Not static gains in worst case
    disp('Worst case perturbation not a pure gain');
end

% Take an example and get step and impulse responses
id2 = find(Rset <= 0.5, 1, 'last'); Rid2 = Rset(id2);
istep = zeros(lg, lent); iimp = istep;
for ii = 1:lg
    % Linear systems simulations
    istep(ii, :) = lsim(G(ii, id2), ustep, t);
    iimp(ii, :) = lsim(G(ii, id2), uimp, t);
end

% For step response get expected steady state error
ess = zeros(1, lg); etheo = ess;
for ii = 1:lg
    % Practical steady state error
    ess(ii) = ustep(1)*(1 - evalfr(G(ii, id2), 0));
    % Theoretical computation from R
    switch(ctrlType)
        case 1 
            etheo(ii) = ustep(1)*(1 - 1/(1 - Rid2));
        case 2
            etheo(ii) = ustep(1)*((1 - Rid2)/(2 - Rid2));
    end
end

%% Fig 2 of margins with constant controllers

figure('Renderer', 'painters', 'Position', [10 10 1000 1000]);
cols = {'g', 'b', grey1, grey2, 'r'};
% Distributions of generation time
subplot(2, 2, 1); hold on;
for ii = 1:lg
    plot(t, w(ii, :), 'Color', cols{ii}, 'LineWidth', 2);
end
h = gca; plot([g0 g0], [0 h.YLim(2)], 'k--', 'LineWidth', 2);
hold off; box off; grid off;
xlabel('$t$ (days)', 'FontSize', fnt); xlim([0 20]);
ylabel('$w(t)$', 'FontSize', fnt); ylim([0 0.5]);
% Poles and growth rates (including truncation)
subplot(2, 2, 2); hold on;
for ii = 1:lg
    plot(Rset, pmax(ii, :), 'Color', cols{ii}, 'LineWidth', 2);
    plot(Rset, pred(ii, :), '--', 'Color', cols{ii}, 'LineWidth', 2);
end
plot([min(Rset) 1], [0 0], 'k--', 'LineWidth', 2);
h = gca; plot([1 1], [h.YLim(1) 0], 'k--', 'LineWidth', 2);
hold off; box off; grid off;
xlabel('$R$', 'FontSize', fnt); 
xlim([min(Rset) max(Rset)]);
ylabel('$p, p_1$', 'FontSize', fnt); 
% Margins of gain and 1/R
subplot(2, 2, 3); hold on;
for ii = 1:lg
    plot(Rset, gmarg(ii, :), 'Color', cols{ii}, 'LineWidth', 2);
end
plot(Rset, 1./Rset, 'k--', 'LineWidth', 2);
plot([1 1], [0 1], 'k--', 'LineWidth', 1);
plot([min(Rset) 1], [1 1], 'k--', 'LineWidth', 1);
hold off; box off; grid off;
xlabel('$R$', 'FontSize', fnt); 
xlim([min(Rset) max(Rset)]);
ylabel('$K^*$', 'FontSize', fnt);
% Additional diskmargin component
axes('Position',[0.27 0.27 0.15 0.15]);
hold on;
for ii = 1:lg
    plot(Rset(1:id1), disk(ii, 1:id1), '--', 'Color', cols{ii}, 'LineWidth', 2);
    plot(Rset(1:id1), dworst(ii, 1:id1), 'Color', cols{ii}, 'LineWidth', 2);
end
hold off; box off; grid off;
xlim([min(Rset) 1]);
xlabel('$R$', 'FontSize', fnt); 
ylabel('$D$', 'FontSize', fnt); 

% Step responses of closed loop TFs
subplot(2, 2, 4); hold on;
%plot(t, ustep, 'k--', 'LineWidth', 2);
for ii = 1:lg
    plot(t, istep(ii, :), 'Color', cols{ii}, 'LineWidth', 2);
end
hold off; box off; grid off;
xlim([min(t) max(t)]);
xlabel('$t$ (days)', 'FontSize', fnt); 
ylabel('$i(t)$', 'FontSize', fnt); 
% Additional impulse response component
axes('Position',[0.70 0.18 0.15 0.15]);
hold on;
for ii = 1:lg
    plot(t, iimp(ii, :), 'Color', cols{ii}, 'LineWidth', 2);
end
hold off; box off; grid off;
ylim([0 3]); xlim([min(t) max(t)]);
xlabel('$t$', 'FontSize', fnt); 
ylabel('$i(t)$', 'FontSize', fnt); 

