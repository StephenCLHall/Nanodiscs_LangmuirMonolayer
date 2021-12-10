function CalcConfidenceIntervals(BayesFile,sigLevel)
%{
This script calculates mean and confidence intervals which are dependant on
fitted parameters, but not direcly fit.

    Input arguments:
    BayesFile: 'Filename' containing RasCAL MCMC Output
    sigLevel: 0<sigLevel<1 Sets the confidence interval to calculate

Input arguments can be defined either in the script body, or from the
terminal (if commented out in the script)
%}
%% LOAD THE DATA AND CHOOSE OPTIONS
    BayesFile = 'Bayes_PC_RAFT_SMALP_Nanodiscs.mat';
    output = load(BayesFile);
    output = output.output;
    [iterations,fitParams] = size(output.chain);
    close all;
    
%Set the Significance Level
    sigLevel = 0.95;
    
%% SELECT THE FITTED PARAMETERS

    Tail_Thickness_Initial	= output.chain(:,1);
    Headgroup_Thickness_Initial = output.chain(:,2); 
    Tail_Thickness_Final = output.chain(:,3);
    Headgroup_Thickness_Final = output.chain(:,4);
    Fraction_Lipid_Exchange = output.chain(:,5);
    VF_Polymer = output.chain(:,6);
    
%% CONSTANTS

    %Neutron Scattering lengths (b)
    % N.B. Units Ang
    b_PC_HG = 602.7e-6;
    b_dDM_tails = 5329.76e-6;
    b_hDM_tails = -291.64e-6;

    % Partial Specific Volumes
    % N.B. Units Ang^3
    Vol_PC_HG = 331; %https://doi.org/10.1016/j.chemphyslip.2006.04.002 for PC
    Vol_DM_tails_LE = 780;

    %SLDs
    % N.B. Units Ang^-2
    PC_HG_SLD = b_PC_HG/Vol_PC_HG;
    pSMA_SLD = 1.84e-6;

    dDM_tails_SLD_LE = b_dDM_tails / Vol_DM_tails_LE;
    hDM_tails_SLD_LE = b_hDM_tails / Vol_DM_tails_LE;

%% MODEL CALCULATIONS
%% INITIAL BILAYER AREA PER MOLECULE CALCULATIONS

% Calculation of Area Per Molecule
    Lipid_APM = (Vol_DM_tails_LE) ./ (Tail_Thickness_Initial);

%Calculate the apparent area per molecule of the lipid HG (including water)
    HG_APM_Apparent = Vol_PC_HG ./ Headgroup_Thickness_Initial;
    
%Calculate the Volume fraction of Heads in the layer given the APM of the
%headgroups must equal that of the tails;

    HG_Vol_Frac = HG_APM_Apparent ./ Lipid_APM;

%Calculate the Hydration of the headgroups from the volume fractions of the
%lipid heads present in the 'HG' layer. This enaures that the total number
%of headgroups = total number of tails.

    HG_Hydration = 100 - (100 .* HG_Vol_Frac);

%% CALCULATE BILAYER AREA PER MOLECULE TAKING INTO ACCOUNT POLYMER EMBEDDING

% Calculation of Area Per Lipid Molecule, accounting for the volume
% occupied by polymer chains in the tails
    VF_Lipid = 1 - VF_Polymer;

    Lipid_APM_Final = (Vol_DM_tails_LE ./ Tail_Thickness_Final) .* VF_Lipid;

%Calculate the apparent area per molecule of the lipid HG (including water)

    HG_APM_Apparent_Final = (Vol_PC_HG ./ Headgroup_Thickness_Final) .* VF_Lipid;

%Calculate the Volume fraction of Heads in the layer given the APM of the
%headgroups must equal that of the tails
    HG_Partial_Vol_Frac_After = HG_APM_Apparent_Final ./ Lipid_APM_Final;

%Calculate the Volume fraction of solvent taking into account the APM
%calculations and the volume fraction of solvent. 
%N.B. THIS IS NOT NESTED SO BE CAREFUL WITH PARAMETER LIMITS TO KEEP THIS 
%VARIABLE POSITIVE!

    HG_Vol_Frac_After = VF_Lipid .* HG_Partial_Vol_Frac_After;
    
    HG_Solvent_VF_After = VF_Lipid .* (1 - HG_Partial_Vol_Frac_After);
    
    HG_Hydration_After = 100 .* HG_Solvent_VF_After;

%% CALCULATE SLDs

%Calcualte the net SLD of headgroup and tail layers taking into account
%volume fractions of lipid, polymer (and solvent in the case of headgroups)

    Headgroup_SLD = PC_HG_SLD;

    HG_Pol_Combined_VF = VF_Polymer + HG_Vol_Frac_After;

    HG_SLD_After = ((VF_Polymer ./ HG_Pol_Combined_VF) .* pSMA_SLD)...
        + ((HG_Vol_Frac_After ./ HG_Pol_Combined_VF) .* Headgroup_SLD);
    
    hDM_Tail_SLD_After = ((1 - VF_Polymer) .* (((1 - Fraction_Lipid_Exchange) .* hDM_tails_SLD_LE)...
        + (Fraction_Lipid_Exchange .* dDM_tails_SLD_LE))) + (VF_Polymer .* pSMA_SLD);
    
    dDM_Tail_SLD_After = ((1 - VF_Polymer) .* (((1 - Fraction_Lipid_Exchange) .* dDM_tails_SLD_LE) ...
        + (Fraction_Lipid_Exchange .* hDM_tails_SLD_LE))) + (VF_Polymer .* pSMA_SLD);
    
    ddDM_Tail_SLD_After = ((1 - VF_Polymer) .* (dDM_tails_SLD_LE)) + (VF_Polymer .* pSMA_SLD);


%% LIST ALL PARAMETERS
    Ang = char(197);
    
    ParamNames = {...
            ['Initial Tail Thickness / ' Ang ''],...
            ['Initial Headgroup Thickness / ' Ang ''],...
            ['Final Tail Thickness / ' Ang ''],...
            ['Final Headgroup Thickness / ' Ang ''],...
            ['Volume Fraction Polymer'],...
            ['Fraction of Lipid Exchanged'],...
            ['Initial Lipid APM / ' Ang '^2'],...
            ['Intial Headgroup Hydration / %'],...
            ['Final Lipid APM/ ' Ang '^2'],...
            ['Final Headgroup Hydration / %'],...
            ['Final Headgroup SLD / *10^{-6} ' Ang '^{-2}'],...
            ['Final dTail SLD / *10^{-6} ' Ang '^{-2}'],...
            ['Final hTail SLD / *10^{-6} ' Ang '^{-2}'],...
            ['Final ddTail SLD / *10^{-6} ' Ang '^{-2}'],...
            };
    
    nParams = length(ParamNames);
    
    AllParams = zeros(iterations,nParams);
    
    AllParams(:,01) = Tail_Thickness_Initial;
    AllParams(:,02) = Headgroup_Thickness_Initial;
    AllParams(:,03) = Tail_Thickness_Final;
    AllParams(:,04) = Headgroup_Thickness_Final;
    AllParams(:,05) = VF_Polymer;
    AllParams(:,06) = Fraction_Lipid_Exchange;
    AllParams(:,07) = Lipid_APM;
    AllParams(:,08) = HG_Hydration;
    AllParams(:,09) = Lipid_APM_Final;
    AllParams(:,10) = HG_Hydration_After;
    AllParams(:,11) = HG_SLD_After;
    AllParams(:,12) = dDM_Tail_SLD_After;
    AllParams(:,13) = hDM_Tail_SLD_After;
    AllParams(:,14) = ddDM_Tail_SLD_After;
  

%% CALCULATE CONFIDENCE INTERVALS

    Names = {'Parameter', 'Mean', 'Relative_Low_Error', 'Relative_High_Error'};
    
    values = num2cell(zeros(nParams, length(Names)));

    for i = 1:nParams;
        
        thisParam = AllParams(:,i);
        % Histogram the parameter
        h = histogram(thisParam,length(thisParam),...
            'binlimits', [min(thisParam) max(thisParam)]...
            );
        % Calculate the cumulative probability distribution function and normalise
        norm_cpdf = cumsum(h.Values) / length(thisParam);
        [csu, iu, ju] = unique(norm_cpdf);

        % Get the mean
        Mean = mean(thisParam);
        % Get lower bound 
        lowRange = interp1(csu, h.BinEdges(iu), (1-sigLevel));

        % Get upper bound
        highRange = interp1(csu, h.BinEdges(iu), sigLevel);

        % Make table of results
        if i<nParams && Mean<1e-4
            values(i+1,:) = {ParamNames{i} Mean*1e6 (Mean*1e6)-(lowRange*1e6) (highRange*1e6)-(Mean*1e6)};
        elseif i<nParams;
            values(i+1,:) = {ParamNames{i} Mean Mean-lowRange highRange-Mean};
        elseif i>=nParams && Mean<1e-4
            values(i+1,:) = {ParamNames{i} Mean*1e6 (Mean*1e6)-(lowRange*1e6) (highRange*1e6)-(Mean*1e6)};         
        else
            values(i+1,:) = {ParamNames{i} Mean Mean-lowRange highRange-Mean};
        end
    end
    close all;

    Table = cell2table(values(2:end,:));
    Table.Properties.VariableNames = Names;
    
    disp(Table);
end



