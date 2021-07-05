function [bfp, LLVec, completedruns] = find_ML_parameters(data,model,runlist,runmax)
% FIND_ML_PARAMETERS estimates the parameters that maximize the likelihood
% of the data given the parameters.
%
% ================= INPUT VARIABLES ===============
% DATA: struct, participant data
% MODEL: vector of length 4
%     model(1): encoding scheme. 1: V, 2: F
%     model(2): inference scheme: 1: V, 2: F, 3: L, 4: S
%     model(3): decision rule: 1: optimal, 2: max
%     model(4): decision noise: 0: none, 1: local, 2: global 
%       (only 2 is presented in Yoo, Acerbi, & Ma 2021)
% RUNLIST: vector of indices of which run to complete
% RUNMAX: scalar integer, number of total runs
%
% ============== OUTPUT VARIABLES ================
% BFP: maximum-likelihood parameter estimate
% LLVEC: associated likelihood
% COMPLETEDRUNS: indices of completed runs
% 
% aspen h yoo, 2018
% aspen.yoo@nyu.edu

if nargin < 4 || isempty(runmax); runmax = 50; end

% get model fitting settings
condition = data.pres2stimuli;
[logflag,LB,UB,PLB,PUB] = getFittingSettings(model, condition);

% BADS options
options.UncertaintyHandling = 'on';

% Generate set of starting point with a Latin hypercube design
rng(0); % Same set for all
nvars = numel(PLB);
x0_list = lhs(runmax,nvars,PLB,PUB,[],1e3);

filename = sprintf('fits/%s/subj%s_%s_model%d%d%d%d.mat',condition,data.subjid,condition,model(1),model(2),model(3),model(4));

% ibs bads settings
dMat = data.Delta;
rels = unique(data.rel);
blah = data.rel;
for irel = 1:length(rels)
    blah(blah == rels(irel)) = irel;
end
dMat = [dMat blah];

for iter = 1:numel(runlist)
    fprintf('iteration number: %d \n',runlist(iter))
    
    % Fix random seed based on iteration (for reproducibility)
    rng(runlist(iter));
    
    % ibs bad settings
    options_ibs = ibslike('defaults');
    options_ibs.Vectorized = 'on';
    options_ibs.MaxIter = 10000;
    
    x0 = x0_list(runlist(iter),:);
    x0 = x0+rand(size(x0))*0.1-0.05; % adding a tiny bit of noise for jobs that are not converging
    x0 = min(max(x0,LB+1e-3),UB-1e-3); % make sure it is within bounds
    
    % this is for ibs bads
    fun = @(x,y) fun_LL(x,y,model,condition,logflag);
    [xbest,LL,~,~] = bads(@(x) ibslike(fun,x,data.resp,dMat,options_ibs),x0,LB,UB,PLB,PUB,[],options);
    
    xbest(logflag) = exp(xbest(logflag)); % getting parameters back into natural units
    
    % it is necessary to reload the file at every iteration in case multiple processors are
    % saving the file at the same time
    if exist(filename,'file')
        load(filename,'bfp','LLVec','completedruns')
    else
        [bfp, LLVec, completedruns] = deal([]);
    end
    
    % update and save variables
    bfp = [bfp; xbest];
    LLVec = [LLVec; LL];
    completedruns = [completedruns; runlist(iter)];
    save(filename,'bfp','LLVec','completedruns')
end

end
