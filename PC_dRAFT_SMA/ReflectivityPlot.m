function ReflectivityPlot(BayesFile,scale,plottype, ShowXAxisLabels, ShowYAxisLabels)
%{
This script takes the output of the Bayesian analysis file produced by
RasCAL, and plots the data with R error bars, and the fits with 95%
confidence intervals.
Data can be plotted as loglog, semilog, either as R vs Q or RQ^4 vs Q

    Input arguments:
    BayesFile = 'filename' containing RasCAL MCMC Output
    scale = 'semilog' or 'loglog'
    plottype = 'R' or 'RQ4'

These options are for hiding axes labels, which can be useful for producing
stacked multipanel figures. 1 shows the axis labels, 0 hides the axis
labels
    ShowXAxisLabels = 1 or 0
    Show YAxisLabels = 1 or 0
    

Input arguments can be defined either in the script body, or from the
terminal (if commented out in the script below)
%}

%% LOAD THE DATA

    BayesFile = 'Bayes_PC_Polymer.mat';
    output = load(BayesFile);
    output = output.output;
    nCurves = size(output.bestRefs,1);
    close all

%% CHOOSE PLOTTING OPTIONS
% Comment these if intending to define from command line

    scale = 'loglog';
    plottype = 'RQ4';
    ShowXAxisLabels = 1;
    ShowYAxisLabels = 1;

%% SETUP FIGURE AXES

    figure(1);
    clf;
    hold on;
    
    set(gca,...
        'FontSize',14,...
        'FontName', 'Arial',...
        'FontWeight', 'bold',...
        'LineWidth', 1.5,...
        'YMinorTick', 'on',...
        'TickDir', 'both',...
        'YScale','log',...
        'Position',[0.13 0.19 0.75 0.75],'Box','on'...
        );
    
%% CHOOSE THE SCALE, EITHER LOG-LOG OR SEMI-LOG

    if strcmp(scale,'loglog') == 1  
        
        set(gca,'XScale','log','XMinorTick', 'on');   
        
    elseif strcmp(scale,'semilog') == 1   
        
        set(gca,'XScale','linear');        
        ax = gca;
        ax.XAxis.TickLabelFormat = '%.2f';     
        
    else        
        
        disp('ERROR - the scale parameter should be only semilog or loglog')
        return     
        
    end       
    
%% SET AXIS LABELS DEPENDING ON PLOTTING OPTIONS

    axis square;
    xlim([0.009 0.28]);
    Ang = char(197);
    
    if ShowXAxisLabels == 1 && ShowYAxisLabels ==1 
    
        if strcmp(plottype,'R') == 1
            ylim([1e-7 2]);
            h2 = ylabel(['Reflectivity']);
            set(h2, 'Units');
            h = xlabel(['{\itQ_z} / ' Ang '^-^1']);
            set(h, 'Units');
        elseif strcmp(plottype,'RQ4') == 1
            ylim([1e-13 1e-7]);
            h2 = ylabel(['R{\itQ_z}^4 / ' Ang '^{-4}']);
            set(h2, 'Units');
            h = xlabel(['{\itQ_z} / ' Ang '^-^1']);
            set(h, 'Units');
        else
            disp('ERROR - the plottype parameter should only be R or RQ4')
            return
        end
        
    elseif ShowXAxisLabels == 1
        
        if strcmp(plottype,'R') == 1
            ylim([1e-7 2]);
            set(gca, 'yticklabel', '');
            h = xlabel(['{\itQ_z} / ' Ang '^-^1']);
            set(h, 'Units');
        elseif strcmp(plottype,'RQ4') == 1
            ylim([1e-13 1e-7]);
            set(gca, 'yticklabel', '');
            h = xlabel(['{\itQ_z} / ' Ang '^-^1']);
            set(h, 'Units');
        else
            disp('ERROR - the plottype parameter should only be R or RQ4')
            return
        end
        
    elseif ShowYAxisLabels == 1
        
        if strcmp(plottype,'R') == 1
            ylim([1e-7 2]);
            h2 = ylabel(['Reflectivity']);
            set(h2, 'Units');
            set(gca, 'xticklabel', '');
        elseif strcmp(plottype,'RQ4') == 1
            ylim([1e-13 1e-7]);
            h2 = ylabel(['R{\itQ_z}^4 / ' Ang '^{-4}']);
            set(h2, 'Units');
            set(gca, 'xticklabel', '');
        else
            disp('ERROR - the plottype parameter should only be R or RQ4')
            return
        end
        
    else
        
        if strcmp(plottype,'R') == 1
            ylim([1e-7 2]);
            set(gca, 'yticklabel', '');
            set(gca, 'xticklabel', '');
        elseif strcmp(plottype,'RQ4') == 1
            ylim([1e-13 1e-7]);
            set(gca, 'yticklabel', '');
            set(gca, 'xticklabel', '');
        else
            disp('ERROR - the plottype parameter should only be R or RQ4')
            return
        end
        
    end

    
    
%% CHOOSE COLOURS AND ORDER TO PLOT CURVES

    blue = [0 0 255] ./ 255;
    lightblue = [0 204 255] ./ 255;
    red = [255 0 0] ./ 255;
    orange = [255 128 0] ./ 255;
    green = [0 175 0] ./ 255;
    lightgreen = [124 255 0] ./ 255;
    
    
    colours = [blue; lightblue; green; lightgreen; red; orange];
    
    order = [1 2 3 4 5 6];
    
    offsets = [1 1 1 1 1 1];
    
%% LOOP TO ADD PLOTS

    for n = 1:length(order)
        i = order(n);
        thisData = output.data{i};
        thisFit = output.bestRefs{i};
        minCurve = output.refBounds{1}{i};
        maxCurve = output.refBounds{2}{i};
        col = colours(n,:);
        offset = offsets(i);
        addPlot(thisData, thisFit, minCurve, maxCurve, col, offset, plottype);
        hold on
    end
end

function addPlot(thisData, thisFit, minCurve, maxCurve, col, offset, plottype)
%{
Function to add plots of reflectivity data, errors, best fit and 95% 
confidence intervals
%}
%% TRANSFORM REFLECTED INTENSITIES DEPENDING ON PLOTTING OPTIONS

    if strcmp(plottype,'R') == 1
        
        %Data
        xData = thisData(:,1);
        yData = thisData(:,2) .* offset;
        yErr = thisData(:,3);
        
        %Fit
        xFit = thisFit(:,1);
        yFit = thisFit(:,2) .* offset;
        minCurve = minCurve .* offset;
        maxCurve = maxCurve .* offset;
        xp = [xData',fliplr(xData')];
        yp = [minCurve',fliplr(maxCurve')];
    
    elseif strcmp(plottype,'RQ4') == 1
        
        %Data
        xData = thisData(:,1);
        yData = thisData(:,2) .*(xData.^4) .* offset;
        yErr = thisData(:,3).*(xData.^4) .* offset;
        
        %Fit
        xFit = thisFit(:,1);
        xFit4 = xFit .^ 4;
        yFit = thisFit(:,2) .*(xFit4) .* offset;
        minCurve = minCurve .* offset;
        maxCurve = maxCurve .* offset;
        xp = [xData',fliplr(xData')];
        yp = [minCurve',fliplr(maxCurve')].*(xp.^4);
        
    else
        
        disp('ERROR - the plottype parameter should only be R or RQ4')
        return
        
    end
%% ADD DATA TO FIGURE    

    figure(1);
    
    %Define a patch, representing 95% confidence intervals
    hp = patch(xp,yp,1,...
        'FaceAlpha',0.4,...
        'facecolor',col,...
        'linestyle','none'...
        );
    
    %Define data as an errorbar object
    hData = errorbar(xData,yData,yErr,'ok',...
        'MarkerSize', 5,...
        'MarkerEdgeColor', [0 0 0],...
        'MarkerFaceColor', col...
        );
    
    %Define a plot representing the best fit
    hFit = plot(xFit,yFit,'-',...
        'Linewidth', 2,...
        'color', col...
        );
    
    %Add data, fit and confidence intervals to figure
    set(hFit);
    set(hData);
    set(hp);

end
