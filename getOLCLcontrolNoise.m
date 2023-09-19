% Compute transfer functions and margins with a controller applied
function [z, p, L, G, marg, pmax] = getOLCLcontrolNoise(R0, W, K, rho, H, ctrlType)

% Assumptions and notes
% - given R0 and W(s) get open L(s) and closed loop G(s) TFs
% - also considers reporting rate rho and delay distrib H(s)
% - applies a controller K(s) or if K(s) = 1 gives raw margins

% Construct open and closed loop transfer functions
switch(ctrlType)
    case 1
        % Controlled epidemic model with true infections all seen
        L = -R0*W*K; G = feedback(1, L); 
    case 2
        % Control with prsymptomatic spread or surveillance noise
        L = -R0*W*(1+rho*H*(K - 1)); G = feedback(1, L); 
end

% Closed loop zeros and poles sorted by real part
p{1} = esort(pole(G)); z{1} = esort(zero(G)); 
% Maximum real poles from systems
try
    pmax(1) = p{1}(max(p{1} == real(p{1})));
catch 
    %warning('Max pole may not be real');
    pmax(1) = max(real(p{1})); 
end

% All margins from open loop transfer functions
d = allmargin(L); 
marg.d = min(d.DelayMargin(d.DelayMargin > 0));
marg.g = min(d.GainMargin(d.GainMargin > 0));
marg.ph = min(d.PhaseMargin(d.PhaseMargin > 0));

% Check for infinite margins
if isempty(marg.d)
    marg.d = inf;
end
if isempty(marg.g)
    marg.g = inf;
end
if isempty(marg.ph)
    marg.ph = inf;
end

% Also compute a disc margin and worst perturbation
disk = diskmargin(L); marg.disk = disk.DiskMargin;
marg.worst = disk.WorstPerturbation;
% If perturbation is static gain just save gain
[zdata, pdata] = tfdata(marg.worst);
if cell2mat(pdata) == 1
    marg.worst = cell2mat(zdata);
end

% Remove any unstable margin warnings
warning('off', 'Control:analysis:MarginUnstable');