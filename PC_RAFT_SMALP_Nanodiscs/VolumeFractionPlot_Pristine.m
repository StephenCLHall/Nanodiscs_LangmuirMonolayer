function VolumeFractionPlot_Pristine(BayesFile, SigLevel)
%{
This script uses the output of the RasCAL Bayesian analysis to plot the
best fit and confidence intervals of the model as a volume fraction profile
rather than traditional SLD profiles.

Here, the volume fraction profile of the pristine monolayer, BEFORE POLYMER
ADDITION, is calcualted.

    Input arguments:
    BayesFile: 'Filename' containing RasCAL MCMC Output
    sigLevel: 0<sigLevel<1 Sets the confidence interval to calculate

Input arguments can be defined either in the script body, or from the
terminal (if commented out in the script)

Note: This function has to process a large amount of data, so will take a
few seconds to generate the volume fraction plot
%}
%% 1. INITIAL DEFINITIONS
%{
This first section load the bayesian output file from rascal and sets up
the structures which will later be filled by the MCMC chains for the
parameters.
%}
%% Load the Data, set plotting range and options
    
    BayesFile = 'Bayes_PC_RAFT_SMALP_Nanodiscs.mat';
    SigLevel = 0.95;

    output = load(BayesFile);
    output = output.output;
    close all;


    TotalRange = -15:40;
    
    ShowXAxisLabels = 1;
    ShowYAxisLabels = 1;

    [iterations,nParams]=size(output.chain);
    
%% Setup relevant parameters and layers
%{  
List of relevant parameters to calculate volume fraction profile
This should NOT need to be modified
%}
    Parameters = {'Thickness',...
                  'Volume_Fraction_1',...
                  'Volume_Fraction_2',...
                  'LHS_Roughness',...
                  'RHS_Roughness',...
                  'Center'...
                  };

%{
List of Layers modelled in your analysis, including bulk phases
This WILL need to be modified for each analysis
Should be ordered BulkIn ---> BulkOut
%}
    Layers = {'Air',...
              'Tails',...
              'Heads',...
              'Solvent'...
              };

    nLayers = length(Layers);

    nParameters = length(Parameters);

    Bulk_Thick = 100;

%{ 
This loop generates a structure containing fields corresponding to
Parameters (as listed above). Within each Parameter field, a column is
inserted for each layer (as listed in Layers above), with some default
values inserted. 
This shouldn't need changing but defaults can be changed if needed
%}    

    for i = 1:nParameters
        for j = 1:nLayers
            if i==1 && j==1
                %Set Substrate_Thickness = Bulk_Thick
                AllLayers.(Parameters{i}) = ones(iterations,1)*Bulk_Thick;
            elseif i==3 && j==1
                %Sets Volume_Fraction_2 = 0 as default
                AllLayers.(Parameters{i}) = zeros(iterations,1);
            elseif j==1
                AllLayers.(Parameters{i}) = ones(iterations,1);
            elseif i==1 && j==nLayers
                %Sets Bulk Solvent Thickness = Bulk_Thick
                ThisLayer = ones(iterations,1)*Bulk_Thick;
                AllLayers.(Parameters{i}) = [AllLayers.(Parameters{i}), ThisLayer];
            elseif i==2 && j==nLayers || i==3
                %Sets all Volume_Fraction_2 of all leyers to 0 by default
                %and bulk solvent to 0 (required for solvation calculation
                %later)
                ThisLayer = zeros(iterations,1);
                AllLayers.(Parameters{i}) = [AllLayers.(Parameters{i}), ThisLayer];
            else
                %All other layer parameters set to 1 by default
                ThisLayer = ones(iterations,1);
                AllLayers.(Parameters{i}) = [AllLayers.(Parameters{i}), ThisLayer];
            end
        end
    end

%% 2. PARAMETER DEFINITIONS AND CONVERSIONS
%{  
This is a list of parameters, as listed in your custom model script.
After Bayesian analysis, the MCMC chain of each parameter will be
listed as a column in output.chain of the Bayesian output file.
The order will be the same as in your custom model script.

SLDs are not required for volume fraction representations, so comment
out any unrequired parameters to save memory and speed up the
calculation. Similarly, and unused parameters or variables can be
removed.

Any 'Hydration' Parameters need to be converted to a volume fraction
Should any calculations be done in the model script which are
dependent on fitted parameters which are then used to calculate the
SLD profile, include those calculations here as well.
%}

%% Fitting Parameters

    Tail_Thickness	= output.chain(:,1);
    Headgroup_Thickness = output.chain(:,2); 

%% Constants

    VF_PG = 0;
    Substrate_Roughness = 3.5;

    Vol_PC_HG = 331;
    Vol_DM_tails_LE = 780;
    Vol_PG_HG = 289; %PG: https://doi.org/10.1016/j.bbamem.2012.05.007
                     %PC: https://doi.org/10.1016/j.chemphyslip.2006.04.002

%% Model Calculations
%{    
Any 'Hydration' Parameters need to be converted to a volume fraction
Should any calculations be done in the model script which are
dependent on fitted parameters which are then used to calculate the
SLD profile, include those calculations here as well.
%}    
%INITIAL BILAYER AREA PER MOLECULE CALCULATIONS:

% Calculation of Area Per Molecules
    Lipid_APM = (Vol_DM_tails_LE)./(Tail_Thickness);

%Calculate the apparent area per molecule of the lipid HG (including water)
    HG_APM_Apparent = ((1-VF_PG)*(Vol_PC_HG./Headgroup_Thickness))+((VF_PG)*(Vol_PG_HG./Headgroup_Thickness));

%Calculate the Volume fraction of Heads in the layer given the APM of the
%headgroups must equal that of the tails;
    HG_Vol_Frac = HG_APM_Apparent./Lipid_APM;

%Calculate the Hydration of the headgroups from the volume fractions of the
%lipid heads present in the 'HG' layer. This enaures that the total number
%of headgroups = total number of tails.

    HG_Hydration = 1-HG_Vol_Frac;

    Tails_Hydration = 0;

%% 3. LAYERS AND COMPONENT PARAMETERS
%{    
This section populates the Parameter fields of the AllLayers structure
with the MCMC Chain corresponding to each layer.

Due to personal preferences in custom model writing, at this stage,
this needs to be entered manually for Thickness, Volume Fractions and
RHS_Roughness.
Other Parameters will be automatically populated later.
%}    
%% Populate the AllLayers structure with relevant parameter chains

% Thicknesses
    AllLayers.Thickness(1:iterations,2) = Tail_Thickness;
    AllLayers.Thickness(1:iterations,3) = Headgroup_Thickness;

% Volume Fraction (Component 1)
    AllLayers.Volume_Fraction_1(1:iterations,3) = HG_Vol_Frac;

% RHS_Roughness
    AllLayers.RHS_Roughness(1:iterations,1) = Substrate_Roughness;
    AllLayers.RHS_Roughness(1:iterations,2) = Substrate_Roughness;
    AllLayers.RHS_Roughness(1:iterations,3) = Substrate_Roughness;

%% Calculate LHS Roughness and Layer Centers

% Automatically Populates LHS Roughnesses based on the RHS Roughnesses 
% of the previous layer
    for i = 2:nLayers
        AllLayers.LHS_Roughness(1:iterations,(i)) = AllLayers.RHS_Roughness(1:iterations,(i-1));
    end


% Automatically populate Centers of each layer, starting from the bulk-in
% interface at 0
    for i = 1:(nLayers)
        if i == 1
            AllLayers.Center(:,i) = -(AllLayers.Thickness(1:iterations,i))/2;
        elseif i == 2
            AllLayers.Center(:,i) = AllLayers.Thickness(1:iterations,i)/2;
        else
            AllLayers.Center(1:iterations,i) = sum(AllLayers.Thickness(:,2:(i-1)),2)...
                + (AllLayers.Thickness(:,i)/2);
        end
    end

%% 4. CALCULATION OF COMPONENT VOLUME FRACTIONS
%{    
This section Calculates the volume fraction of each component.
This is done for each iteration within the MCMC chain, for each layer
and generates a structure, containing volume fractions of two
components in each layer as a function of distance

Uses asymconvstep, which produces a step function convoluted with 
differnt error functions on each side. (From RasCAL 2014):

asymconvstep (x,xw,xcen,s1,s2,h)
      x = vector of x values
     xw = Width of step function
   xcen = Centre point of step function
     s1 = Roughness parameter of left side
     s2 = Roughness parameter of right side
      h = Height of step function.

r = xcen + (xw/2);
l = xcen - (xw/2);

a = (x-l)./((2^0.5)*s1);
b = (x-r)./((2^0.5)*s2);

f = (h/2)*(erf(a)-erf(b));
%}
%% Calculate volume fractions of each component on a layer-by-layer basis

    for i = 1:nLayers
        AllComponents.(Layers{i}).Component_1 = zeros(iterations,length(TotalRange));
        AllComponents.(Layers{i}).Component_2 = zeros(iterations,length(TotalRange));
    end

    for i = 1:(nLayers-1)

        if any(AllLayers.Volume_Fraction_2(:,i)) == 1                        
            for k = 1:iterations
                AllComponents.(Layers{i}).Component_1(k,:) = asymconvstep(TotalRange,...
                                                         AllLayers.Thickness(k,i),...
                                                         AllLayers.Center(k,i),...
                                                         AllLayers.LHS_Roughness(k,i),...
                                                         AllLayers.RHS_Roughness(k,i),...
                                                         AllLayers.Volume_Fraction_1(k,i)...
                                                         );

                AllComponents.(Layers{i}).Component_2(k,:) = asymconvstep(TotalRange,...
                                                         AllLayers.Thickness(k,i),...
                                                         AllLayers.Center(k,i),...
                                                         AllLayers.LHS_Roughness(k,i),...
                                                         AllLayers.RHS_Roughness(k,i),...
                                                         AllLayers.Volume_Fraction_2(k,i)...
                                                         );

            end
        else

            for k = 1:iterations

                AllComponents.(Layers{i}).Component_1(k,:) = asymconvstep(TotalRange,...
                                                         AllLayers.Thickness(k,i),...
                                                         AllLayers.Center(k,i),...
                                                         AllLayers.LHS_Roughness(k,i),...
                                                         AllLayers.RHS_Roughness(k,i),...
                                                         AllLayers.Volume_Fraction_1(k,i)...
                                                         );
            end
        end
    end

%% Calculate the volume fraction of solvent throughout the interface
%{    
Up to now, we have negated the solvent contribution. In order to avoid
artefacts, this next section calculates the total volume fraction
accounted by non-solvent components as a function of distance such that
solvent must account for the rest.
%}

%Check if there is already a solvent field in the AllComponents
%structure from a previous run. If there is, this will remove that
%field before recalculating.


    if isfield(AllComponents,'Solvent') == 1
        AllComponents = rmfield(AllComponents,'Solvent');
        Non_Solvent_Components = fieldnames(AllComponents);
        AllComponents.Solvent.Component_1 = zeros(iterations,length(TotalRange));
        AllComponents.Solvent.Component_2 = zeros(iterations,length(TotalRange));
    else
        AllComponents.Solvent.Component_1 = zeros(iterations,length(TotalRange));
        AllComponents.Solvent.Component_2 = zeros(iterations,length(TotalRange));
    end

    Individual_Vol_Fracs = zeros(length(Non_Solvent_Components),length(TotalRange),iterations);

    for k = 1:iterations
        for i = 1:length(Non_Solvent_Components)

                Individual_Vol_Fracs(i,:,k) = AllComponents.(Non_Solvent_Components{i}).Component_1(k,:);...
                                      + AllComponents.(Non_Solvent_Components{i}).Component_2(k,:);
        end

        AllComponents.Solvent.Component_1(k,:) = 1 - (sum(Individual_Vol_Fracs(:,1:(length(TotalRange)),k)));

    end

%% Remove redundant fields
%{    
For Ease of plotting, now remove any fields for redundant components
In this case, Peripheral_Protein and Protein_Binding_Gap is now 
combined with component_2 from the Membrane layer, so these two fields
can be removed, and all component 2 fields are now redundant, so can be
removed
%}    
    ComponentList = fieldnames(AllComponents);
    for i = 1:length(ComponentList)
        AllComponents.(ComponentList{i}) = rmfield(AllComponents.(ComponentList{i}),'Component_2');
    end

%% 5. CALCULATION OF VOLUME FRACTION MEAN AND CONFIDENCE INTERVALS
%{    
Now the volume fractions for each interfacial component have been
determined, this section will calculate the mean and confidence
interval for the volume fraction of each component as a function of
distance from the interface.
%}

%% Select what components to calculate VFs for and setup stucture

%If you still have multiple components in each layer that you want to
%plot, then add 'Component_2' into this cell
    Components = {'Component_1'};

%TODO
%Make this more seamless!
    Layers = {'Air', 'Tails', 'Heads', 'Solvent'};

%Create Struct to fill with VF mean and CIs
    for i = 1:length(Layers)
        for c = 1:length(Components)
            VolumeFractions.(Layers{i}).(Components{c}) = zeros(length(TotalRange),3);
        end
    end

%% Calculate mean and confidence intervals from each MCMC iteration
    for i = 1:length(Layers)
        for j = 1:length(TotalRange)
            for c = 1:length(Components)

            ThisComponent = AllComponents.(Layers{i}).(Components{c})(:,j);
            % histogram the values of the component
            h = histogram(ThisComponent,...
                          iterations,...
                          'binlimits', [min(ThisComponent) max(ThisComponent)]);
            % calculate the cumulative probability distribution function and normalise
            norm_cpdf = cumsum(h.Values)/iterations;
            % Find the unique values of the histogram and the indices
            [csu, iu, ju] = unique(norm_cpdf);

            %Calculate the mean
            Mean = mean(ThisComponent);

            %Calculate the low uncertainty bound
            Low_Uncertainty = interp1(csu, h.BinEdges(iu), (1-SigLevel));

            %Calculate the high uncertainty bound
            High_Uncertainty = interp1(csu, h.BinEdges(iu), SigLevel);

            %Pupulate the VolumeFractons Structure with Mean, Low, High
            %columns, with each row corresponding to a different distance
            VolumeFractions.(Layers{i}).(Components{c})(j,:) = [Mean, Low_Uncertainty, High_Uncertainty];

            end
        end
    end

%% 6. PLOTTING
%{    
Now we have calculated everything to construct the component volume
fraction profiles, this next section will do the plotting and
formatting.
%}

%% Set up the figure
    Ang = char(197);
    clf;
    figure(1);
    hold on
    xlim([min(TotalRange) max(TotalRange)]);
    ylim([0 1.05]);

    set(gca,...
        'FontSize', 14,...
        'FontName', 'Arial',...
        'FontWeight', 'bold',...
        'LineWidth', 1.5,...
        'TickDir', 'both',...
        'Position',[0.13 0.19 0.75 0.75],...
        'Layer','top',...
        'Box', 'on'...
        );
    axis square
    grid off

%% Setup axes labels depending on if they are to be shown or not

    if ShowXAxisLabels == 1 && ShowYAxisLabels == 1
        ax1 = ylabel(['Component volume fraction']);
        ax2 = xlabel(['Distance from interface / ' Ang '']);
        set(ax1, 'Units');
        set(ax2, 'Units');
        
    elseif ShowXAxisLabels == 1
        ax2 = xlabel(['Distance from interface / ' Ang '']);
        set(ax2, 'Units');
        set(gca, 'yticklabel', '');
        
    elseif ShowYAxisLabels == 1
        ax1 = ylabel(['Component volume fraction']);
        set(ax1, 'Units');
        set(gca, 'xticklabel', '');
        
    else
        set(gca, 'xticklabel', '');
        set(gca, 'yticklabel', '');
    end
     

%% Choose colours for each component to be plotted

    blue = [0 0 255]./255;
    lightblue = [80 200 255]./255;
    orange = [255 120 0]./255;
    gold = [204 167 28]./255;
    yellow = [255 208 28]./255;
    red = [255 0 0]./255;
    green = [0 192 0]./255;
    darkgreen = [0 150 0]./255;
    lightgreen = [104 200 15]./255;
    purple = [165 0 255]./255;
    lightpurple = [231 98 255]./255;
    grey = [150 150 150]./255;
    black = [0 0 0]./255;
    ibm_blue = [100 143 255]./255;
    ibm_purple = [150 94 240]./255;
    ibm_pink = [220 38 127]./255;
    ibm_orange = [254 97 0]./255;
    ibm_gold = [255 176 0]./255;

    Components_To_Plot = {'Air', 'Tails', 'Heads', 'Solvent'};

    Colours = [black; ibm_orange; lightgreen; ibm_blue];
    
    Linestyle = {'-', '-', '-', '-'};

%% Setup a loop to do the plotting

     for i= 1:length(Components_To_Plot)
        MeanCurve = VolumeFractions.(Components_To_Plot{i}).Component_1(:,1);
        MinCurve = VolumeFractions.(Components_To_Plot{i}).Component_1(:,2);
        MaxCurve = VolumeFractions.(Components_To_Plot{i}).Component_1(:,3);
        hyErrorBar = [MaxCurve; flipud(MinCurve)];
        hxErrorBar = [TotalRange'; fliplr(TotalRange)'];

        if i == 1
            hFit = plot(TotalRange,MeanCurve,Linestyle{i},...
                'Linewidth', 2,...
                'color', Colours(i,:)...
                );
            set(hFit);
        else 
            hFit = plot(TotalRange,MeanCurve,Linestyle{i},...
                'Linewidth', 2,...
                'color', Colours(i,:)...
                );
            set(hFit);
            
            hPatch = patch(hxErrorBar, hyErrorBar, 1,...
                'FaceAlpha', 0.5,...
                'facecolor', Colours(i,:),...
                'linestyle', 'none'...
                );
           
           set(hPatch);

        end
        hold on
    end
end