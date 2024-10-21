function c = tapas_nassar_delta_config
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Contains the configuration for the Approximate Delta learning model from Nassar et al. 2010.
%
% Nassar, M. R., Wilson, R. C., Heasly, B., & Gold, J. I. (2010). An approximately 
% Bayesian delta-rule model explains the dynamics of belief updating in a changing environment. 
% Journal of Neuroscience, 30(37), 12366-12378.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% hazExp  is estimated in 'logit-space' because it is bounded inside the unit
% interval.
%
% 'Logit-space' is a logistic sigmoid transformation of native space
% 
% tapas_logit(x) = ln(x/(1-x)); x = 1/(1+exp(-tapas_logit(x)))
%
% Any of the parameters can be fixed (i.e., set to a fixed value) by setting the variance of their
% prior to zero. To fix vhat_0 to 0.5 set the mean as well as the variance of the prior to zero.
% --------------------------------------------------------------------------------------------------
% Copyright (C) 2012-2013 Christoph Mathys, TNU, UZH & ETHZ
%
% This file is part of the HGF toolbox, which is released under the terms of the GNU General Public
% Licence (GPL), version 3. You can redistribute it and/or modify it under the terms of the GPL
% (either version 3 or, at your option, any later version). For further details, see the file
% COPYING or <http://www.gnu.org/licenses/>.

% Config structure
c = struct;

% Model name
c.model = 'nassar_delta';

% Initial v
c.vhat_0mu = 50;
c.vhat_0sa = 0;

% Initial log sd
c.logN_0mu = log(5);
c.logN_0sa = 0;
% Alpha
c.logithazExpmu = tapas_logit(0.3, 1);
c.logithazExpsa = 1;

% Gather prior settings in vectors
c.priormus = [
    c.vhat_0mu,...
    c.logN_0mu,...
    c.logithazExpmu,...
         ];

c.priorsas = [
    c.vhat_0sa,...
    c.logN_0sa,...
    c.logithazExpsa,...
         ];

% Model function handle
c.prc_fun = @tapas_nassar_delta;

% Handle to function that transforms perceptual parameters to their native space
% from the space they are estimated in
c.transp_prc_fun = @tapas_nassar_delta_transp;

return;
