function SLDplot(bayesfile, ShowXAxisLabels, ShowYAxisLabels)
%{
This script takes the output of the Bayesian analysis file produced by
RasCAL, and SLD profiles with 95% confidence intervals.
Data can be plotted as with or without axis labels

    Input arguments:
    bayesfile = 'filename.mat'
    xlabels = 1 or 0
    ylabels = 1 or 0

Input arguments can be defined either in the script body, or from the
terminal (if commented out in the script)
%}

%% LOAD THE DATA
% Comment bayesfile if intending to define from command line

    bayesfile = 'Bayes_PC_SMALP_Nanodiscs.mat';
    output = load(bayesfile);
    output = output.output;
    close all

%% SHOW AXES LABELS?
% Useful for producing fingures without labels so they can be stacked as
% multipanel figures with shared axes
% Comment these if intending to define from command line
%
% 1 = show axis labels
% 0 = do not show axis labels

    ShowXAxisLabels = 1;
    ShowYAxisLabels = 1;
    
%% SETUP THE AXES

    figure(1);
    clf;
    hold on;
    
    xlim([-15 40]);
    ylim([-1 6.5]);
    
    set(gca,...
        'FontSize',14,...
        'FontName', 'Arial',...
        'FontWeight', 'bold',...
        'LineWidth', 1.5,...
        'TickDir', 'both',...
        'Position',[0.13 0.19 0.75 0.75],...
        'Layer','top',...
        'Box', 'on'...
        );
    
    ax = gca;
    ax.XAxis.TickLabelFormat = '%.0f';
    ax.YAxis.TickLabelFormat = '%.0f';
    
    axis square;
    
    Ang = char(197);
 
%% SETUP AXES LABELS

    if ShowXAxisLabels == 1 && ShowYAxisLabels == 1
        
        hy = ylabel(['SLD / \times10^{-6} ' Ang '^-^2']);
        hx = xlabel(['Distance from interface / ' Ang '']);
        set(hy, 'Units');
        set(hx,'Units');
        
    elseif ShowXAxisLabels == 1
        
        hx = xlabel(['Distance from interface / ' Ang '']);
        set(hx, 'Units');
        set(gca, 'yticklabel', '');
        
    elseif ShowYAxisLabels == 1
        
        hy = ylabel(['SLD / \times10^{-6} ' Ang '^-^2']);
        set(hy, 'Units');
        set(gca, 'xticklabel', '');
        
    else
        
        set(gca, 'yticklabel', '',...
            'xticklabel', '');
        
    end

%% CHOOSE COLOURS AND LINE ORDER

    blue = [0 0 255] ./ 255;
    lightblue = [0 204 255] ./ 255;
    red = [255 0 0] ./ 255;
    orange = [255 128 0] ./ 255;
    green = [0 175 0] ./ 255;
    lightgreen = [124 255 0] ./ 255;

    colours = [red; orange; green; lightgreen; blue; lightblue];
    
    order = [1 2 3 4 5 6];
    
    linestyles = {'-', '-', '--', '-', '-', '-'};
    
    offset = 99;
    
%% Setup the loop to add plots   

    for n = 1:length(order)
        i = order(n);
        thisSLD = output.bestSlds{i};
        minCurve = output.sldBounds{1}{i} .* 1e6;
        maxCurve = output.sldBounds{2}{i} .* 1e6;
        col = colours(i,:);
        Linestyle = linestyles{i};
        addPlot(thisSLD, minCurve, maxCurve, col, offset, Linestyle);
        
    end
end

function addPlot(thisSLD, minCurve, maxCurve, col, offset, Linestyle)
%{
    Function to add plots of SLD profiles of best fit and associated 95% confidence
    intervals
%}
%% 
    xData = thisSLD(:,1) - offset;
    yData = thisSLD(:,2) * 1e6;
    
    if length(xData) < length(minCurve);
        minCurve = minCurve(1:length(xData));
        maxCurve = maxCurve(1:length(xData));
    else
        xData = xData(1:length(minCurve));
    end
    
    yp = [minCurve',fliplr(maxCurve')];
    xp = [xData',fliplr(xData')];
    
    hp = patch(xp, yp, 1,...
        'FaceAlpha', 0.3,...
        'facecolor', col,...
        'linestyle', 'none'...
        );
           
    hFit = plot(xData, yData, 'k',...
        'Linewidth', 2, ...
        'linestyle', Linestyle, ...
        'color', col ...
        );
    
    set(hFit);
    set(hp);

end

