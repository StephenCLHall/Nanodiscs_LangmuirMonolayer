function [output,sub_rough] = Model_MonolayerPolymer(params,bulk_in,bulk_out,contrast)
%{
RasCAL model to define a lipid monolayer before and after interaction with
deuterated SMA.

The model constrains headgroup hydration by area per molecule, allows the
monolayer structure to chage as a result of polymer interaction, and allows
separate embedding of the polymer into the headgroup and tail layers.
%}
%% PARAMETERS

Tail_Thickness = params (1);

Headgroup_Thickness = params(2);

VF_PG = params (3); % N.B. The volume fraction of PG is included as a 
                    % parameter, but should be FIXED AS A CONSTANT IN THE 
                    % GUI based on prior knowledge about the sample

VF_Polymer_Heads = params (4);

VF_Polymer_Tails = params (5);

Tail_Thickness_Final = params (6);

Headgroup_Thickness_Final = params (7);

%% FIXED VALUES

Substrate_Roughness = 3.5;

%Neutron Scattering lengths (b)
% N.B. Units Ang
b_PC_HG = 602.7e-6;
b_PG_HG = 717.18e-6;
b_dDM_tails = 5329.76e-6;
b_hDM_tails = -291.64e-6;

% Partial Specific Volumes
% N.B. Units Ang^3
Vol_PC_HG = 331;
Vol_DM_tails_LC = 663;
Vol_DM_tails_LE = 780;
Vol_PG_HG = 289; %See https://doi.org/10.1016/j.bbamem.2012.05.007 for PG
                 %https://doi.org/10.1016/j.chemphyslip.2006.04.002 for PC

%SLDs
% N.B. Units Ang^-2
PC_HG_SLD = b_PC_HG / Vol_PC_HG;
PG_HG_SLD = b_PG_HG / Vol_PG_HG;
pSMA_alt_SLD = 4.353e-6;
pS_SLD = 5.598e-6;
pSMA_SLD = 4.768e-6;

dDM_tails_SLD_LE = b_dDM_tails / Vol_DM_tails_LE;
dDM_tails_SLD_LC = b_dDM_tails / Vol_DM_tails_LC;
hDM_tails_SLD_LE = b_hDM_tails / Vol_DM_tails_LE;
hDM_tails_SLD_LC = b_hDM_tails / Vol_DM_tails_LC;

%% INITIAL BILAYER AREA PER MOLECULE CALCULATIONS

% Calculation of Area Per Molecule
Lipid_APM = (Vol_DM_tails_LC) / (Tail_Thickness);

%Calculate the apparent area per molecule of the lipid HG (including water)
HG_APM_Apparent = ((1 - VF_PG) * (Vol_PC_HG / Headgroup_Thickness))...
    + ((VF_PG) * (Vol_PG_HG / Headgroup_Thickness));

%Calculate the Volume fraction of Heads in the layer given the APM of the
%headgroups must equal that of the tails;

HG_Vol_Frac = HG_APM_Apparent / Lipid_APM;

%Calculate the Hydration of the headgroups from the volume fractions of the
%lipid heads present in the 'HG' layer. This enaures that the total number
%of headgroups = total number of tails.

HG_Hydration = 100 - (100 * HG_Vol_Frac);

Headgroup_SLD = (VF_PG * PG_HG_SLD) + ((1 - VF_PG) * PC_HG_SLD);

Tails_Hydration = 0;

%% CALCULATE BILAYER AREA PER MOLECULE TAKING INTO ACCOUNT POLYMER EMBEDDING

% Calculation of Area Per Lipid Molecule, accounting for the volume 
%occupied by polymer chains in the tails
Lipid_APM_After = (Vol_DM_tails_LC / Tail_Thickness_Final) * (1 - VF_Polymer_Tails);

%Calculate the volume fraction of Heads + Solvent
VF_Heads_Solvated = 1 - VF_Polymer_Heads;

%Calculate the Volume fraction of Heads in the layer given the APM of the
%headgroups must equal that of the tails;
HG_APM_Apparent_After = (((1 - VF_PG) * (Vol_PC_HG / Headgroup_Thickness_Final))...
    + ((VF_PG) * (Vol_PG_HG / Headgroup_Thickness_Final))) * VF_Heads_Solvated;

HG_Partial_Vol_Frac_After = HG_APM_Apparent_After / Lipid_APM_After;

%Calculate the Hydration of the headgroups from the volume fractions of the
%lipid heads present in the 'HG' layer. This enaures that the total number
%of headgroups = total number of tails.

HG_Vol_Frac_After = VF_Heads_Solvated * HG_Partial_Vol_Frac_After;
HG_Solvent_Vol_Frac_After = VF_Heads_Solvated * (1 - HG_Partial_Vol_Frac_After);

%% CALCULATE SLDs

%Calcualte the net SLD of headgroup and tail layers taking into account
%volume fractions of lipid, polymer (and solvent in the case of headgroups)

%NB. To make the maths work for calculating the SLDs, include the hydration
%in the SLD calculation for the headgroup. Make sure to set the hydration
%in HG_AFTER to 0.

HG_SLD_After = (HG_Solvent_Vol_Frac_After*bulk_out(contrast))+(VF_Polymer_Heads*pSMA_SLD)+(HG_Vol_Frac_After*Headgroup_SLD);
hDM_Tail_SLD_After = ((1-VF_Polymer_Tails)*hDM_tails_SLD_LC)+(VF_Polymer_Tails*pSMA_SLD);
dDM_Tail_SLD_After = ((1-VF_Polymer_Tails)*dDM_tails_SLD_LC)+(VF_Polymer_Tails*pSMA_SLD);

%% DEFINE LAYERS

%SYNTAX:
%LAYER = [Thickness SLD Roughness Hydration 1];

%--- INITIAL LAYERS (Before Nanodiscs)---

dTAILS = [Tail_Thickness dDM_tails_SLD_LC Substrate_Roughness Tails_Hydration 1];

hTAILS = [Tail_Thickness hDM_tails_SLD_LC Substrate_Roughness Tails_Hydration 1];

HG = [Headgroup_Thickness Headgroup_SLD Substrate_Roughness HG_Hydration 1];

%--- FINAL LAYERS (After Nanodiscs)---

dTAILS_AFTER = [Tail_Thickness_Final dDM_Tail_SLD_After Substrate_Roughness Tails_Hydration 1];

hTAILS_AFTER = [Tail_Thickness_Final hDM_Tail_SLD_After Substrate_Roughness Tails_Hydration 1];

HG_AFTER = [Headgroup_Thickness_Final HG_SLD_After Substrate_Roughness 0 1];

%% DEFINE LAYERS PRESENT IN EACH SAMPLE

%SYNTAX:
%output = [LAYER1; LAYER2; LayerN];

switch contrast
    case {1}
        output =  [hTAILS; HG];
    case {2}
        output =  [hTAILS_AFTER; HG_AFTER];
    case {3}
        output =  [dTAILS; HG];
    case {4}
        output =  [dTAILS_AFTER; HG_AFTER];
    case {5}
        output =  [dTAILS; HG];
    case {6}
        output =  [dTAILS_AFTER; HG_AFTER];
end
  
%% ROUGHNESS RESAMPLING

resamp = 1;

switch resamp
    case 1
        
        
        if size(output,2) == 5
            for i = 1:size(output,1);
                thisSLD = output(i,2);
                thisHydr = output(i,4)/100;
                thisSLD = (thisHydr*bulk_out(contrast)) + (1-thisHydr)* thisSLD;
                thisRough = output(i,3);
                newOutput(i,:) = [output(i,1) thisSLD thisRough];
            end
        end
        
        output = newOutput;
        
        
        
        [prof,subs_surface] = makeSLDProfile(bulk_in(contrast),bulk_out(contrast),Substrate_Roughness,output,size(output,1),1);
        [xn,yn] = resample_sld(prof,1);
        [xn,yn] = groupSamples(xn,yn,1e-9);
        rrr = ones(length(xn),1)*1e-9;
        output = [xn' yn' rrr];
        
end

sub_rough = Substrate_Roughness;

end

function [x, y] = groupSamples(x,y,tolerance)

 
debug = 0;

numberOfLayers = length(x);
count = 1;
newX = [];
newY = [];

thisLayerx = [x(1)];
thisLayery = [y(1)];
try
for i = 1:numberOfLayers-1;
    %debug
    diff = abs(y(i+1) - y(i));
    if debug; fprintf('Diff:- %5.5g  y(i):- %5.3g   pos:-  %5.2f \n',diff,y(i),i);end
    if diff <= tolerance ;
        thisLayerx = [thisLayerx x(i)];
        thisLayery = [thisLayery (y(i)+(diff/2))];
    else
        if debug; disp('grouping \n'); end
        newX = [newX sum(thisLayerx)];
        newY = [newY mean(thisLayery)];
        thisLayerx = [x(i)];
        thisLayery = [y(i)];
    end
end
catch
    [a,b] = lasterr;
    disp('debug');
end

    
newX = [newX sum(thisLayerx)];
newY = [newY mean(thisLayery)];  
    
x = newX;
y = newY;
end

%-------------------------------------------------------------------------

function [out,subs_surface] = makeSLDProfile(nbair,nbsub,ssub,layers,numberOfLayers,nrepeats);
 
layerThicks = sum(layers(:,1));
totalRange = layerThicks + 100;
x = 0:totalRange;
boxCen = 0;
boxWidth = 100;
boxRough = ssub;%layers(1,3);
subBox = asymconvstep(x,boxWidth,boxCen,boxRough,boxRough,nbair);
lastBoxEdge = boxCen + (boxWidth/2);
lastLayRough = ssub;
subs_surface = lastBoxEdge;
for i = 1:numberOfLayers;
    thisLayThick = layers(i,1);
    thisLaySLD = layers(i,2);
    thisLayRough = layers(i,3);
    %             if i<numberOfLayers
    %                 nextLayRough = layers(i+1,3);
    %             elseif (i == numberOfLayers)
    %                 nextLayRough = layers(i,3);
    %             end
    thisBoxCentre = lastBoxEdge + (thisLayThick/2);
    Lays(:,i) = asymconvstep(x,thisLayThick,thisBoxCentre,lastLayRough,thisLayRough,thisLaySLD);
    lastBoxEdge = thisBoxCentre + (thisLayThick/2);
    lastLayRough = thisLayRough;
end
thisLayThick = (x(end)-lastBoxEdge)*2;
thisLaySLD = nbsub;
thisLayRough = lastLayRough;
nextLayRough = lastLayRough;
thisBoxCentre = x(end);
Lays(:,(numberOfLayers+1)) = asymconvstep(x,thisLayThick,thisBoxCentre,thisLayRough,nextLayRough,thisLaySLD);

SLD = sum(Lays,2);

SLD = SLD(:) + subBox(:);

out = [x(:),SLD(:)];

 
end

 
function [x_new_tot,y_new_tot] = resample_sld(sld,y_min_step);

x = sld(:,1);
prof = sld(:,2);

 
thisx = x(1);
count = 1;
i = 2;
while i < length(x);
    this_step = prof(i) - prof(i-1);
    if abs(this_step) < y_min_step;
        %Keep original points
        x_new(count) = x(i);
        y_new(count) = prof(i);
        count = count + 1;
        i = i + 1;
    else
        if this_step > 0
            newsteps = prof(i-1):y_min_step:prof(i);
        else
            newsteps = prof(i-1):-y_min_step:prof(i);
        end
        new_xs = interp1(prof(i-1:i),x(i-1:i),newsteps);
        x_new = [x_new new_xs];
        y_new = [y_new newsteps];
        count = length(x_new);
        i = i + 1;
    end
end
%x_new_tot = x_new(:);
%y_new_tot = y_new(:);

%Make it into a 'histogram'
for i =1:length(x_new)-1;
    x_new_tot(i) = x_new(i + 1) - x_new(i);
    y_new_tot(i) = y_new(i);
end

 
%Remove any zero thickness layers...
good_layers = find(x_new_tot ~= 0);
x_new_tot = x_new_tot(good_layers);
y_new_tot = y_new_tot(good_layers);

 
end
