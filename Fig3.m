% Upset the margins by applying K(s) with dynamics
clearvars; close all; clc;

% Assumptions and notes
% - show impact of delay margins when K(s) has dynamics
% - examine margins for fixed r but varying RW

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


%% Lead compensations that destabilises the controller

% Define an R and magnitude of K(s)
R = 4; Kgain = 1/9; ctrlType = 1;
% Some controller choices
%K = [(1+s)/(1+8*s), (1+8*s)/(1+s), (1+8*s)/(s^2 + s + 1), 1 + 8/s, 1 + 8*s]; 
K = [(1+s)/(1+8*s), (1+8*s)/(1+s), (1+8*s)/(s^2 + s + 1), (1/((5/50)*s+1))^50, 1]; 

% Find maximum magnitude of controllers
lenk = length(K); magK = zeros(1, lenk); wmag = magK;
for ii = 1:lenk
    % Maximum magnitude from bode
    [magtemp, ~, wtemp] = bode(Kgain*K(ii)); 
    [magK(ii), idw] = max(magtemp);
    wmag(ii) = wtemp(idw);
end
disp(['Max |K(s)| are: ' num2str(magK)]);

% Initialise TFs and margins
L = s*ones(lg, lenk); G = L; marg = cell(lg, lenk);
% Store margins from first case (including diskmargin)
gmarg = zeros(lg, lenk); dmarg = gmarg; disk = gmarg;
% Zeros, poles and margins of systems
z = marg; p = marg; pmax = marg; dworst = marg; pred = gmarg;

% Obtain TFs and their properties
for jj = 1:lg
    for ii = 1:lenk
        % Main function constructing open and closed loops
        [z{jj, ii}, p{jj, ii}, L(jj, ii), G(jj, ii), marg{jj, ii}, ...
            pmax{jj, ii}] = getOLCLcontrolNoise(R, W(jj), Kgain*K(ii), 1, 1, ctrlType);
        % Store minimum gain and delay margins
        gmarg(jj, ii) = marg{jj, ii}.g; 
        dmarg(jj, ii) = marg{jj, ii}.d;
        % Disk margin and worst perturbation
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

% Step and impulse responses under third controller
istep = zeros(lg, lent); iimp = istep; istep0 = istep; iimp0 = istep;
for ii = 1:lg
    % Linear systems simulations
    istep(ii, :) = lsim(G(ii, 2), ustep, t);
    iimp(ii, :) = lsim(G(ii, 2), uimp, t);
    istep0(ii, :) = lsim(G(ii, 1), ustep, t);
    iimp0(ii, :) = lsim(G(ii, 1), uimp, t);
end


%% Delay to destabilise the system

% Delay in the loop TF
Ldel = L*exp(-3.5*s); Gdel = cell(lg, lenk); margdel = Gdel;
for jj = 1:lg
    for ii = 1:lenk
        % Closed loop TF and margins
        Gdel{jj, ii} = feedback(1, Ldel(jj, ii));
        margdel{jj, ii} = allmargin(Ldel(jj, ii));
    end
end

% Step and impulse responses under third controller
istepdel = zeros(lg, lent); iimpdel = istepdel; 
istepdel0 = istepdel; iimpdel0 = istepdel;
for ii = 1:lg
    % Linear systems simulations
    istepdel(ii, :) = lsim(Gdel{ii, 2}, ustep, t);
    iimpdel(ii, :) = lsim(Gdel{ii, 2}, uimp, t);
    istepdel0(ii, :) = lsim(Gdel{ii, 1}, ustep, t);
    iimpdel0(ii, :) = lsim(Gdel{ii, 1}, uimp, t);
end

% Inverse Laplace transform some controllers
syms('s'); Ks = Kgain*[(1+s)/(1+8*s), (1+8*s)/(1+s)];
Kt = ilaplace(Ks);

%% Fig 3 of destabilising controllers and dynamics

figure('Renderer', 'painters', 'Position', [10 10 1000 600]);
cols = {'g', 'b', grey1, grey2, 'r'};
subplot(1, 2, 1); hold on;
for ii = 2:lg
    plot(t, istep0(ii, :), '--', 'Color', cols{ii}, 'LineWidth', 2);
    plot(t, istep(ii, :), 'Color', cols{ii}, 'LineWidth', 2);
end
xlabel('$t$ (days)', 'FontSize', fnt); 
ylim([10 25]); ylabel('$i(t)$', 'FontSize', fnt);

subplot(1, 2, 2); hold on;
for ii = 2:lg
    plot(t, istepdel0(ii, :), '--', 'Color', cols{ii}, 'LineWidth', 2);
    plot(t, istepdel(ii, :), 'Color', cols{ii}, 'LineWidth', 2);
end
xlabel('$t$ (days)', 'FontSize', fnt); 
ylim([10 25]); ylabel('$i(t)$', 'FontSize', fnt);





