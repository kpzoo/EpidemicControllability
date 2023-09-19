% Generation time distributions and Laplace transform
function [w, W, wstat] = generationLaplace(GT, GTtype, s, t)

% Assumptions and notes
% - one of 4 types of generation time distribution
% - generation time parameters for exp, det, gamma and bimodal
% - outputs PDF w(t) and Laplace transform W(s)

% Parametrise different generation time structures
switch(GTtype)
    case 1
        % Deterministic generation time
        %W = exp(-s*GT.mean); W = pade(W, 2);
        dshape = 80;
        W = (1 + (GT.mean/dshape)*s)^(-dshape);
        w = gampdf(t, dshape, GT.mean/dshape);
    case 2
        % Gamma generation time (scale=1 is exp)
        W = (1 + GT.scale*s)^(-GT.shape);
        w = gampdf(t, GT.shape, GT.scale);
    case 3
        % Exponential generation time
        W = (1 + GT.mean*s)^(-1);
        w = gampdf(t, 1, GT.mean);
    case 4
        % Bimodal distribution is sum of gamma
        W = 0.5*((1 + GT.scale(1)*s)^(-GT.shape(1))) + ...
            0.5*((1 + GT.scale(2)*s)^(-GT.shape(2)));
        w = 0.5*gampdf(t, GT.shape(1), GT.scale(1)) + ...
            0.5*gampdf(t, GT.shape(2), GT.scale(2));
end

% Statistics of every distribution
wstat.norm = trapz(t, w); wstat.mean = trapz(t, w.*t);
wstat.var = trapz(t, w.*(t.^2)) - wstat.mean^2;
% Order and zero freq of Laplace (and MGF)
wstat.normLap = evalfr(W, 0); wstat.order = order(W);