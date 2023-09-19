% Generation time distributions and Laplace transform for two types
function [w, W, wstat] = generationGammaWeighted(g1param, g2param, eps, s, t, lent)

% Assumptions and notes
% - weighted gamma distributions by eps and 1-eps
% - outputs PDF w(t) and Laplace transform W(s)
% - shape-scale formulation of gamma distributions

% Component gamma distributions and mixture
W = s*ones(1, 3); w = zeros(3, lent); wstat = cell(1, 3);

% Laplace transforms
W(1) = (1 + g1param(2)*s)^(-g1param(1));
W(2) = (1 + g2param(2)*s)^(-g2param(1));
W(3) = W(1)*eps + W(2)*(1 - eps);

% Time domain w(t) 
w(1, :) = gampdf(t, g1param(1), g1param(2));
w(2, :) = gampdf(t, g2param(1), g2param(2));
w(3, :) = eps*w(1, :) + (1 - eps)*w(2, :);

% Statistics of all distributions
for ii = 1:3
    % Statistics of every distribution
    wstat{ii}.norm = trapz(t, w(ii, :)); wstat{ii}.mean = trapz(t, w(ii, :).*t);
    wstat{ii}.var = trapz(t, w(ii, :).*(t.^2)) - wstat{ii}.mean^2;
    % Order and zero freq of Laplace (and MGF)
    wstat{ii}.normLap = evalfr(W(ii), 0); wstat{ii}.order = order(W(ii));
end


