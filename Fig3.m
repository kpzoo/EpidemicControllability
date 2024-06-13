% Reproduces Figure 3 
clearvars; close all; clc;

% Assumptions and notes
% - show impact of delay margins when K(s) has dynamics
% - upset the margins by applying K(s) with dynamics

% Figure defaults
set(groot, 'defaultAxesTickLabelInterpreter', 'latex', 'defaultLegendInterpreter', 'latex');
fnt = 24; grey1 = 0.5*ones(1, 3); grey2 = 0.8*ones(1, 3);
set(0, 'defaultTextInterpreter', 'latex', 'defaultAxesFontSize', fnt);

% Complex s and times for any simulations
s = tf('s'); dt = 0.02; ts = 0:dt:100; lent = length(ts);
% Time for checking integrals of controllers
t = 0.0001:0.0001:200;
% Step and impulse input across times
ustep = 100*ones(1, lent); uimp = [100*ones(1, 50), zeros(1, lent-50)];

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
    [w(ii, :), W(ii), wstat{ii}] = generationLaplace(GT, GTtype(ii), s, ts);
end


%% Lead compensations that destabilises the controller

% Define an R and magnitude of K(s)
R = 4; Kgain = 1/8; ctrlType = 1; ctrlSet = 1;

% Controller sets to consider with inverse Laplace
switch(ctrlSet)
    case 1
        % Version of exp type with extra poles
        tau = 2; K = [(1+g0*s)/((1+tau*s)^4), (1+g0*s)/((1+(tau/4)*s)^4)]; 
        syms('s'); taudel = 3.5;
        Ks = ilaplace(Kgain*[(1+g0*s)/((1+tau*s)^4), (1+g0*s)/((1+(tau/4)*s)^4)]); 

        % Other option
        %tau = 8; K = [(1+s)/((1+tau*s)*(1+s)), (1+tau*s)/(1+s)^2];
        %Ks = ilaplace(Kgain*[(1+s)/(1+tau*s)/(1+s), (1+tau*s)/(1+s)^2]); 
    case 2
        % Modification of actual generation time distribution type 2
        tau = 3.5; K = [1, ((1 + (g0/3)*s)^3)/(1 + (tau/6)*s)^6];
        syms('s'); taudel = 2; 
        Ks = ilaplace(Kgain*[1, ((1 + (g0/3)*s)^3)/(1 + (tau/6)*s)^6]);
end
% Integrate contoller
Kint1 = eval(Ks(1)); Kint2 = eval(Ks(2)); s = tf('s'); 
int1 = trapz(t, Kint1); int2 = trapz(t, Kint2);
disp(['Integrals of k(t) are: ' num2str(int1) ' and ' num2str(int2)]);

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
z = marg; p = marg; pmax = marg; dworst = marg; 

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
    istep(ii, :) = lsim(G(ii, 2), ustep, ts);
    iimp(ii, :) = lsim(G(ii, 2), uimp, ts);
    istep0(ii, :) = lsim(G(ii, 1), ustep, ts);
    iimp0(ii, :) = lsim(G(ii, 1), uimp, ts);
end


%% Delay to destabilise the system

% Delay in the loop TF
Ldel = L*exp(-taudel*s); Gdel = cell(lg, lenk); margdel = Gdel;
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
    istepdel(ii, :) = lsim(Gdel{ii, 2}, ustep, ts);
    iimpdel(ii, :) = lsim(Gdel{ii, 2}, uimp, ts);
    istepdel0(ii, :) = lsim(Gdel{ii, 1}, ustep, ts);
    iimpdel0(ii, :) = lsim(Gdel{ii, 1}, uimp, ts);
end


%% Fig 3 of destabilising controllers and dynamics

figure('Position', [10 10 1000 600]);
cols = {'g', 'b', grey1, grey2, 'r'};
subplot(1, 2, 1); hold on;
for ii = 2:lg
    plot(ts, istep0(ii, :), '--', 'Color', cols{ii}, 'LineWidth', 2);
    plot(ts, istep(ii, :), 'Color', cols{ii}, 'LineWidth', 2);
end
xlabel('$t$ (days)', 'FontSize', fnt); 
ylabel('$i(t)$', 'FontSize', fnt); ylim([ustep(1) 2.5*ustep(1)]);

subplot(1, 2, 2); hold on;
for ii = 2:lg
    plot(ts, istepdel0(ii, :), '--', 'Color', cols{ii}, 'LineWidth', 2);
    plot(ts, istepdel(ii, :), 'Color', cols{ii}, 'LineWidth', 2);
end
xlabel('$t$ (days)', 'FontSize', fnt); 
ylabel('$i(t)$', 'FontSize', fnt); ylim([ustep(1) 2.5*ustep(1)]);
