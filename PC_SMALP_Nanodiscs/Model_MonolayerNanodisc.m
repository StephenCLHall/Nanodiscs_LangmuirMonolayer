function [output,sub_rough] = Model_MonolayerNanodisc(params,bulk_in,bulk_out,contrast)
%{
RasCAL model to define a lipid monolayer before and after interaction with
pre-assembled nanodiscs.

The model constrains headgroup hydration by area per molecule, allows the
monolayer structure to chage as a result of polymer interaction, and allows
embedding of the polymer into the headgroup and tail layers. Due to the
lack of contrast between hydrogenous polymers and hydrogenous PC
headgroups, a single volume fration of polymer is fitted across the whole
monolayer.
%}
%% PARAMETERS

    Tail_Thickness_Initial = params (1);

    Headgroup_Thickness_Initial = params(2);

    Tail_Thickness_Final = params(3);

    Headgroup_Thickness_Final = params(4);

    Fraction_Lipid_Exchange = params (5);

    VF_Polymer = params (6);

%% FIXED VALUES

    Substrate_Roughness = 3.5;

%Neutron Scattering lengths (b)
% N.B. Units Ang
    b_PC_HG = 602.7e-6;
    b_dDM_tails = 5329.76e-6;
    b_hDM_tails = -291.64e-6;

% Partial Specific Volumes
% N.B. Units Ang^3

    Vol_PC_HG = 331; %https://doi.org/10.1016/j.chemphyslip.2006.04.002 for PC
    Vol_DM_tails_LC = 663;
    Vol_DM_tails_LE = 780;

%SLDs
% N.B. Units Ang^-2
    PC_HG_SLD = b_PC_HG/Vol_PC_HG;
    pSMA_SLD = 1.84e-6;

    dDM_tails_SLD_LE = b_dDM_tails / Vol_DM_tails_LE;
    dDM_tails_SLD_LC = b_dDM_tails / Vol_DM_tails_LC;
    hDM_tails_SLD_LE = b_hDM_tails / Vol_DM_tails_LE;
    hDM_tails_SLD_LC = b_hDM_tails / Vol_DM_tails_LC;

%% INITIAL BILAYER AREA PER MOLECULE CALCULATIONS

% Calculation of Area Per Molecules

    Lipid_APM = (Vol_DM_tails_LE) / (Tail_Thickness_Initial);

%Calculate the apparent area per molecule of the lipid HG (including water)
    HG_APM_Apparent = Vol_PC_HG/Headgroup_Thickness_Initial;

%Calculate the Volume fraction of Heads in the layer given the APM of the
%headgroups must equal that of the tails;

    HG_Vol_Frac = HG_APM_Apparent / Lipid_APM;

%Calculate the Hydration of the headgroups from the volume fractions of the
%lipid heads present in the 'HG' layer. This enaures that the total number
%of headgroups = total number of tails.

    HG_Hydration = 100 - (100 * HG_Vol_Frac);


    Tails_Hydration = 0;

%% CALCULATE BILAYER AREA PER MOLECULE TAKING INTO ACCOUNT POLYMER AND LIPID EXCHANGE

% Calculation of Area Per Lipid Molecule, accounting for the volume
% occupied by polymer chains in the tails
    VF_Lipid = 1 - VF_Polymer;

    Lipid_APM_Final = (Vol_DM_tails_LE / Tail_Thickness_Final) * VF_Lipid;

%Calculate the apparent area per molecule of the lipid HG (including water)

    HG_APM_Apparent_Final = (Vol_PC_HG / Headgroup_Thickness_Final) * VF_Lipid;

%Calculate the Volume fraction of Heads in the layer given the APM of the
%headgroups must equal that of the tails
    HG_Partial_Vol_Frac_After = HG_APM_Apparent_Final/Lipid_APM_Final;

%Calculate the Volume fraction of solvent taking into account the APM
%calculations and the volume fraction of solvent. 
%N.B. THIS IS NOT NESTED SO BE CAREFUL WITH PARAMETER LIMITS TO KEEP THIS 
%VARIABLE POSITIVE!

    HG_Vol_Frac_After = VF_Lipid*HG_Partial_Vol_Frac_After;
    HG_Solvent_VF_After = VF_Lipid*(1-HG_Partial_Vol_Frac_After);

% This is to print out a 'check' in the terminal that both the total volume
% fraction of all components in the headgroup layer is equal to 1, and that
% the headgroup solvation >= 0

    Check = VF_Polymer+HG_Solvent_VF_After+HG_Vol_Frac_After;

    if Check ==1 && HG_Solvent_VF_After >= 0
        disp({'Total Volume Fraction in headgroups =' Check ''});
        disp({'Headgroup Solvation = ' HG_Solvent_VF_After*100 ''});

    elseif Check == 1 && HG_Solvent_VF_After < 0
        disp({'Total Volume Fraction in headgroups =' Check ''});
        disp({'Headgroup Solvation = ' HG_Solvent_VF_After*100 '.'...
            'This value is negative. Adjust parameter limits for volume fraction of polymer.'});

    elseif Check ~= 1
        disp({'Total Volume Fraction in headgroups =' Check '. Check model!'});
    end

%% CALCULATE SLDs

%Calcualte the net SLD of headgroup and tail layers taking into account
%volume fractions of lipid, polymer (and solvent in the case of headgroups)

%NB. To make the maths work for calculating the SLDs, include the hydration
%in the SLD calculation for the headgroup. Make sure to set the hydration
%in HG_AFTER to 0.

    HG_SLD_After = (HG_Solvent_VF_After * bulk_out(contrast))...
        + (VF_Polymer * pSMA_SLD) + (HG_Vol_Frac_After * PC_HG_SLD);

    hDM_Tail_SLD_After = ((1 - VF_Polymer) * (((1 - Fraction_Lipid_Exchange) * hDM_tails_SLD_LE)...
        + (Fraction_Lipid_Exchange * dDM_tails_SLD_LE))) + (VF_Polymer * pSMA_SLD);
    
    dDM_Tail_SLD_After = ((1 - VF_Polymer) * (((1 - Fraction_Lipid_Exchange) * dDM_tails_SLD_LE) ...
        + (Fraction_Lipid_Exchange * hDM_tails_SLD_LE))) + (VF_Polymer * pSMA_SLD);
    
    ddDM_Tail_SLD_After = ((1 - VF_Polymer) * (dDM_tails_SLD_LE)) + (VF_Polymer * pSMA_SLD);


%% DEFINE LAYERS

%SYNTAX:
%LAYER = [Thickness SLD Roughness Hydration 1];

%--- INITIAL LAYERS (Before Nanodiscs)---

    dTAILS = [Tail_Thickness_Initial dDM_tails_SLD_LE Substrate_Roughness Tails_Hydration 1];

    hTAILS = [Tail_Thickness_Initial hDM_tails_SLD_LE Substrate_Roughness Tails_Hydration 1];

    HG = [Headgroup_Thickness_Initial PC_HG_SLD Substrate_Roughness HG_Hydration 1];
    
%--- FINAL LAYERS (After Nanodiscs)---

    ddTAILS_AFTER = [Tail_Thickness_Final ddDM_Tail_SLD_After Substrate_Roughness Tails_Hydration 1];

    HG_AFTER = [Headgroup_Thickness_Final HG_SLD_After Substrate_Roughness 0 1];

    dTAILS_AFTER = [Tail_Thickness_Final dDM_Tail_SLD_After Substrate_Roughness Tails_Hydration 1];

    hTAILS_AFTER = [Tail_Thickness_Final hDM_Tail_SLD_After Substrate_Roughness Tails_Hydration 1];


%% DEFINE LAYERS PRESENT IN EACH SAMPLE

%SYNTAX:
%output = [LAYER1; LAYER2; LayerN];

    switch contrast
        case (1)
            output =  [dTAILS; HG]; 
        case (2)
            output =  [ddTAILS_AFTER; HG_AFTER];
        case {3}
            output =  [dTAILS; HG];
        case {4}
            output =  [dTAILS_AFTER; HG_AFTER];
        case {5}
            output =  [hTAILS; HG];
        case {6}
            output =  [hTAILS_AFTER; HG_AFTER];

    end

  
%% ROUGHNESS MICROSLICING
  
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
