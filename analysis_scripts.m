%% ====================================================================
%                      COLLECT / FIT YOUR OWN DATA!
% =====================================================================


%% run experiment

run_experiment

% 1. you will be prompted to give subject initials, you can give any string
% 2. You will be prompted to give the Experiment ID. You should input 
%       R for Reliability (main experimental block) 
%       T for Threshold block
%       P for Practice block
% 3. You will be prompted on the condition. You should input
%       E for Ellipse 
%       L for Line
%
% Note that the low-reliability ellipse eccentricity is determined from the
% Threshold block, so that block is required before running the main 
% Reliability block. 


%% fit models to data

clear all

subjid = 'S02';             % subject id
condition = 'Ellipse';      % condition: 'Ellipse' or 'Line'
model = [2 2 3 2];          % [encoding inference decision rule decision noise]
% MODEL is a 1 x 4 vector containing indices corresponding to the 
% [encoding, inference, decision rule, and decision noise] settings desired
%
%  encoding: encoding noise. 1: VP, 2: FP
%  inference: inferred/assumed noise. 1: VP, 2: FP, 3: diss assumed for 
%       ellipse and line, 4: ignore completely
%  decision rule: decision rule. 1: optimal, 2: max
%  decision noise: decision noise: 0: none, 1: local, 2: global
%       (publication only used global decision noise)


% optimization settings
runmax = 20;                % how many times to run optimization (used to sample parameter space using lhs)
runlist = 1;                % which index (1:runmax) to run

load(sprintf('data/%s/%s_%s_simple.mat',condition,subjid,condition))
find_ML_parameters(data,model,runlist,runmax)


%% ====================================================================
%                   ANALYZE DATA, AS IN PUBLICATION
% =====================================================================

%% bayesian model selection: BMS. calculate exceedance probabilities
clear all

condition = 'Line';
load(sprintf('fits/%s/bfp_%s.mat',condition,condition));
nTrials = 2000;

% calculate AIC, AICc, and BIC
AICMat = 2*bsxfun(@plus,LLMat,nParamsVec');
AICcMat = bsxfun(@plus,AICMat,((2.*nParamsVec.*(nParamsVec+1))./(nTrials-nParamsVec-1))');
BICMat = 2*bsxfun(@plus,LLMat,nParamsVec'*log(nTrials));

Nsamp = 1e6;
do_plot = 1;
[alpha,exp_r,xp,pxp,bor] = spm_BMS(-0.5*AICcMat', Nsamp, do_plot)

% OUTPUT:
% alpha   - vector of model probabilities
% exp_r   - expectation of the posterior p(r|y)
% xp      - exceedance probabilities
% pxp     - protected exceedance probabilities
% bor     - Bayes Omnibus Risk (probability that model frequencies 
%           are equal)


%% model families: log marginal likelihoods 

clear all

condition = 'Line';
load(sprintf('fits/%s/bfp_%s.mat',condition,condition));
nTrials = 2000;
nSubjs = 13;

% calculate AIC, AICc, and BIC
AICMat = 2*bsxfun(@plus,LLMat,nParamsVec');
AICcMat = bsxfun(@plus,AICMat,((2.*nParamsVec.*(nParamsVec+1))./(nTrials-nParamsVec-1))');
BICMat = 2*bsxfun(@plus,LLMat,nParamsVec'*log(nTrials));

modelMat(:,4) = (modelMat(:,1) == modelMat(:,2))+1;

nCrits = 15;
critmax = 40;
critVec = linspace(0,critmax,nCrits);

legendlabels{1} = {'V','F'};
legendlabels{2} = {'V','F','L','S'};
legendlabels{3} = {'O','M'};
legendlabels{4} = {'inference mismatched', 'inference matched'};

% calculate log marginals for each factor
[L,M,CIs] = deal(cell(1,4));
for ifactor = 1:4
    
    condVec = unique(modelMat(:,ifactor));
    C = nchoosek(condVec,2);
    nComps = size(C,1);
    L{ifactor} = nan(nComps,nSubjs);
    M{ifactor} = nan(1,nComps);
    CIs{ifactor} = nan(2,nComps);
    for icomp = 1:nComps % each possible value for each factor
        idx1 = modelMat(:,ifactor) == C(icomp,1);
        idx2 = modelMat(:,ifactor) == C(icomp,2);

        level1 = -0.5.*AICcMat(idx1,:); 
        level2 = -0.5.*AICcMat(idx2,:);
        maxx = max(level1); % a somewhat arbitrary constant to help w computational stuff
        blah = log(sum(exp(bsxfun(@minus,level1,maxx)))) - log(sum(exp(bsxfun(@minus,level2,maxx))))...
            +log(size(level2,1)./size(level1,1));
        L{ifactor}(icomp,:) = blah;
        perms = sort(sum(blah(randi(nSubjs,1000,nSubjs)),2));
        CIs{ifactor}(:,icomp) = perms([25 975]);
        M{ifactor}(icomp) = sum(blah);
    end
    
end


%% ====================================================================
%                     CREATE FIGURES, AS IN PUBLICATION
% =====================================================================

%% FIGURE 2: PLOT DATA FOR BOTH EXPERIMENTS 
% this section produces a figure with four subplots:
% - Ellipse: p("change") as a function of magnitude of change, broken up 
%   by number of high-rel ellipses (nHigh)
% - Ellipse: hits and false alarms as a function of nHigh
% - Line: p("change") as a function of magnitude of change, broken up by 
%   nHigh
% - Line: hits and false alarms as a function of nHigh
% 
% (nHigh is the number of high-reliability ellipses in the first stimulus 
% presentation)

clear all
close all

load('plottingsettings.mat')
nBins = 7;
quantilebinning=1;

counter = 0;
for icond = 1:nConds
    condition = conditionVec{icond};
    
    [x_mean, pc_data] = deal(nan(5,nBins,nSubjs));
    [HRallVec,HRlowVec,HRhighVec,FARVec] = deal(nan(nSubjs,5));
    for isubj = 1:nSubjs;
        subjid = subjidVec{isubj};
        
        % load data
        load(sprintf('data/fitting_data/%s_%s_simple.mat',subjid,condition),'data')
        
        % get psychometric function binned data
        [x_mean(:,:,isubj), pc_data(:,:,isubj)] = plot_psychometric_fn(data,nBins,[],quantilebinning);
        
        % get hits/false alarms binned data
        [HRallVec(isubj,:),HRlowVec(isubj,:),HRhighVec(isubj,:),FARVec(isubj,:)] = ...
            plot_HR_FAR(data,[],0);
    end

    % get participant and model means: psychometric fn
    xrange = nanmean(x_mean,3);
    partM = nanmean(pc_data,3);
    partSEM = nanstd(pc_data,[],3)./sqrt(nSubjs-1);
    
    % get participant and model means: hits/false alarms
    m_HRall = mean(HRallVec);
    m_HRlow = mean(HRlowVec);
    m_HRhigh = mean(HRhighVec);
    m_FAR = mean(FARVec);
    sem_HRall = std(HRallVec)./sqrt(nSubjs);
    sem_HRlow = std(HRlowVec)./sqrt(nSubjs);
    sem_HRhigh = std(HRhighVec)./sqrt(nSubjs);
    sem_FAR = std(FARVec)./sqrt(nSubjs);
    
    figure(1);
    % PLOT PSYCHOMETRIC FUNCTIONS
    subplot(1,4,2*counter+1), hold on;
    xlim([-0.2 pi/2+0.2])
    ylim([0 1])
    for ii = 1:5;
        plot(xrange(ii,:),partM(ii,:),'Color',colorMat1(ii,:))
        errorb(xrange(ii,:),partM(ii,:),partSEM(ii,:),'color',colorMat1(ii,:))
    end
    set(gca,'XTick',[0:.25:1].*pi/2,'XTickLabel',{'0','','45','','90'})
    defaultplot
    xlabel('magnitude of change (deg)')
    ylabel('proportion report "change"')
    
    % PLOT HITS FALSE ALARMS
    subplot(1,4,2*counter+2), hold on;
    xlim([-0.5 4.5])
    ylim([0 1])
    plot(0:4,m_HRall,'Color',colorMat2(1,:))
    plot(0:4,m_HRlow,'Color',colorMat2(2,:))
    plot(0:4,m_HRhigh,'Color',colorMat2(3,:))
    plot(0:4,m_FAR,'Color',colorMat2(4,:))
    errorb(0:4,m_HRall,sem_HRall,'color',colorMat2(1,:))
    errorb(0:4,m_HRlow,sem_HRlow,'color',colorMat2(2,:))
    errorb(0:4,m_HRhigh,sem_HRhigh,'color',colorMat2(3,:))
    errorb(0:4,m_FAR,sem_FAR,'color',colorMat2(4,:))
    set(gca,'XTick',0:4,'XTickLabel',0:4,'YTickLabel','')
    defaultplot
    xlabel('number of high reliability ellipses')
    
    counter = counter+1;
    f2 = figure(2);
end

close(f2);

%% FIGURE 5: FACTORIAL MODEL FITS

clear all
condition = 'Line';

figure(2); clf;
modelcolMat = [1:4 9:11; ...
               5:8 12:14 ]';

load(sprintf('analysis/fits/bfp_%s.mat',condition));
load('plottingsettings.mat')

% prediction stuff
nSamples = 2;
nBins = 6;
quantilebinning=1;

% calculated AIC, AICc, and BIC
AICMat = 2*bsxfun(@plus,LLMat,nParamsVec');
AICcMat = bsxfun(@plus,AICMat,((2.*nParamsVec.*(nParamsVec+1))./(nTrials-nParamsVec-1))');
AICcMat = bsxfun(@minus,AICcMat,AICcMat(1,:));

% median
med_AICc = median(AICcMat,2);

% get 95 CI
nModels = size(modelMat,1);
CI_AICc= nan(2,nModels);
for imodel = 2:nModels
    AICcVec = AICcMat(imodel,:);
    
    blah = sort(median(AICcVec(randi(nSubjs,1000,nSubjs)),2));
    CI_AICc(:,imodel) = blah([25 975]);
end

nModels = size(modelcolMat,1);
for icol = 1:2
    imodelVec = modelcolMat(:,icol);
    
    for imodel = 1:nModels;%imodelVec%1:2;%nModels
        modelnum = imodelVec(imodel);
        model = modelMat(modelnum,:);
        
        % load bfp fits
        bfpmat = bfpMat{modelnum};
        nSubj = length(subjidVec);
        
        
        figure(1);
        [x_mean, pc_data, pc_pred] = deal(nan(5,nBins,nSubj));
        [HRallVec,HRlowVec,HRhighVec,FARVec,mod_HRallVec,mod_HRlowVec,mod_HRhighVec,mod_FARVec] = deal(nan(nSubj,5));
        for isubj = 1:nSubj
            subjid = subjidVec{isubj};
            bfp = bfpmat(isubj,:);
            
            % load data
            load(sprintf('data/fitting_data/%s_%s_simple.mat',subjid,condition),'data')
            
            % get predictions
            [LL,p_C_hat] = calculate_LL(bfp,data,model,[],nSamples);
            fprintf('subj %s: %5.2f \n',subjid,LL)
            
            % get psychometric function binned data
            [x_mean(:,:,isubj), pc_data(:,:,isubj), pc_pred(:,:,isubj)] = plot_psychometric_fn(data,nBins,p_C_hat,quantilebinning);
            
            % get hits/false alarms binned data
            [HRallVec(isubj,:),HRlowVec(isubj,:),HRhighVec(isubj,:),FARVec(isubj,:),...
                mod_HRallVec(isubj,:),mod_HRlowVec(isubj,:),mod_HRhighVec(isubj,:),mod_FARVec(isubj,:)] = ...
                plot_HR_FAR(data,p_C_hat,0);
            
        end
        
        % get participant and model means: psychometric fn
        xrange = nanmean(x_mean,3);
        partM = nanmean(pc_data,3);
        partSEM = nanstd(pc_data,[],3)./sqrt(nSubj-1);
        modelM = nanmean(pc_pred,3);
        modelSEM = nanstd(pc_pred,[],3)./sqrt(nSubj-1);
        
        % get participant and model means: hits/false alarms
        m_HRall = mean(HRallVec);
        m_HRlow = mean(HRlowVec);
        m_HRhigh = mean(HRhighVec);
        m_FAR = mean(FARVec);
        sem_HRall = std(HRallVec)./sqrt(nSubj);
        sem_HRlow = std(HRlowVec)./sqrt(nSubj);
        sem_HRhigh = std(HRhighVec)./sqrt(nSubj);
        sem_FAR = std(FARVec)./sqrt(nSubj);
        m_mod_HRall = mean(mod_HRallVec);
        m_mod_HRlow = mean(mod_HRlowVec);
        m_mod_HRhigh = mean(mod_HRhighVec);
        m_mod_FAR = mean(mod_FARVec);
        sem_mod_HRall = std(mod_HRallVec)./sqrt(nSubj);
        sem_mod_HRlow = std(mod_HRlowVec)./sqrt(nSubj);
        sem_mod_HRhigh = std(mod_HRhighVec)./sqrt(nSubj);
        sem_mod_FAR = std(mod_FARVec)./sqrt(nSubj);
        
        % plot
        figure(2);
        
        % PLOT MODEL FITS
        tight_subplot2(nModels,6,imodel,3*icol-2, gutter); hold on
        xlim([-0.2 pi/2+0.2])
        for ii = 1:5;
            plot_summaryfit(xrange(ii,:),partM(ii,:),partSEM(ii,:),modelM(ii,:),...
                modelSEM(ii,:),colorMat1(ii,:),colorMat1(ii,:))
        end
        set(gca,'Xlim',[-pi/16 9*pi/16],'Ylim',[0 1],...
            'XTick',[0:.25:1].*pi/2,'YTick',0:.2:1,...
            'XTickLabel','','YTickLabel','');
        if (imodel==1) && (icol==1); set(gca,'YTickLabel',{0,'','','','',1}); end
        ylabel(modelnamesVec{modelnum})
        if (imodel == nModels); 
            set(gca,'XTickLabel',{'0','','45','','90'}); 
            xlabel('change (deg)')
        end
        
        % PLOT HITS FALSE ALARMS
        tight_subplot2(nModels,6,imodel,3*icol-1, gutter); hold on
        xlim([-0.5 4.5])
        ylim([0 1])
        plot_summaryfit(0:4,m_HRall,sem_HRall,m_mod_HRall,sem_mod_HRall,colorMat2(2,:),colorMat2(2,:));
        plot_summaryfit(0:4,m_HRlow,sem_HRlow,m_mod_HRlow,sem_mod_HRlow,colorMat2(3,:),colorMat2(3,:));
        plot_summaryfit(0:4,m_HRhigh,sem_HRhigh,m_mod_HRhigh,sem_mod_HRhigh,colorMat2(1,:),colorMat2(1,:));
        plot_summaryfit(0:4,m_FAR,sem_FAR,m_mod_FAR,sem_mod_FAR,colorMat2(4,:),colorMat2(4,:));
        set(gca,'XTick',0:4,'YTick',0:.2:1,'XTickLabel','','YTickLabel','');% 0:4)
        if (imodel == nModels); 
            set(gca,'XTickLabel',0:4); 
            xlabel('N_{high}')
        end
        
        % \Delta AICc
        tight_subplot2(nModels,6,imodel,3*icol, gutter); hold on
        fill([0 0 nSubj+1 nSubj+1],...
            [CI_AICc(:,modelnum)' CI_AICc(2,modelnum) CI_AICc(1,modelnum)],0.85*ones(1,3),'EdgeColor','none');
        plot([0 nSubj+1], [med_AICc(modelnum) med_AICc(modelnum)], 'k-'); hold on
        bar(AICcMat(modelnum,:),'FaceColor',[234 191 51]./255,'EdgeColor','none','LineWidth',2)
        set(gca,'Xlim',[0 nSubj+1],'Ylim',[-100 1000],...
                'XTick',[],'XTickLabel',[],...
                'YTick', [-100 0:200:1000],'YTickLabel', '');
        if (imodel==1) && (icol==1); 
            set(gca,'YTickLabel',{'',0,'','','','',1000}); 
            ylabel('\Delta AICc')
        end
        if (imodel == nModels); xlabel('Participant'); end
        defaultplot
    end
end