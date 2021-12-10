function CorrelationPlot(BayesFile)
%{
This script generates a correlation plot of all fitted parameters using the
MCMC chains from RasCALs Bayesian analysis. Histograms are plotted to show
the posterior distribution of fitted parameters, and 2-dimensional
histograms are plotted to show correlations between fitted paramers.
%}
%% LOAD THE DATA

    BayesFile = 'Bayes_PC_RAFT_SMALP_Nanodiscs.mat';
    output = load (BayesFile);
    output = output.output;
    close all

%% SET NUMBER OF PARAMETERS

    nParams = 6;
    order = [1 2 3 4 5 6];
    figure(1); clf
    axis square

%% FORMAT PARAMETER NAMES

    Ang = char(197);
    param1 = ['d_{t} / ' Ang ''];
    param2 = ['d_{h}/ ' Ang ''];
    param3 = ['d_{t}^{eq} / ' Ang ''];
    param4 = ['d_{h}^{eq} / ' Ang ''];
    param5 = ['\chi_{lipid ex}'];
    param6 = ['\chi_{pol}^{eq}'];

    axes_label_names = {param1, param2, param3, param4, param5, param6};

%% SET UP LOOP TO PLOT ALL SUBPLOTS IN CORRECT POSITIONS

    for n = 1:nParams;
        i = order(n);
        thisParam = output.chain(:,i);     
        subplot(nParams,nParams,(i*(nParams+1))-nParams);
        addDistribution(thisParam, nParams, axes_label_names, i);
        
        for m = i:nParams;
            
           if m<nParams;
               j=order(m+1);
               nextParam = output.chain(:,j);
               subplot(nParams,nParams,(i+((j-1)*nParams)));
               addCorrelationPlot(thisParam, nParams, nextParam, axes_label_names, i, j);
           else 
               continue
           end
           
        end
    end
end


function addDistribution (thisParam, nParams, axes_label_names, i);
%{
Plot a histogram showing the posterior distribution of each parameter
%}
%% SETUP THE HISTOGRAM
    h = histogram(thisParam, 50,...
        'binlimits', [min(thisParam), max(thisParam)],...
        'DisplayStyle', 'bar',...
        'Normalization', 'count',...
        'EdgeColor', 'none'...
        );
%% SETUP THE HISTOGRAM AXES

    set(gca,...
        'Layer', 'top',...
        'FontSize', 11,...
        'FontName', 'Arial',...
        'FontWeight', 'bold',...
        'LineWidth', 1.5,...
        'XMinorTick', 'off',...
        'YMinorTick', 'off',...
        'TickDir', 'out',...
        'YTick', [],...
        'TickLength', [0.02, 0.01]...
        );  
        
    axis square;
    xlim([min(thisParam), max(thisParam)]);
    ylim([0,max(h.Values)]);

%% Add AXES LABELS DEPENDING ON THE SUBPLOT POSITION

    if i >= nParams
        set(gca, 'yticklabel', '');
        xlabel(axes_label_names{i});
    else
        set(gca,'yticklabel', '', 'xticklabel', '');
    end

end

%==========================================================================
%Plot bivariate histograms of each pair of parameters and display as a 2D
%heatmap
%==========================================================================

function addCorrelationPlot (thisParam, nParams, nextParam, axes_label_names, i, j)
%{
Plot bivariate histograms of each pair of parameters and display as a 2D
heatmap
%}
%% PLOT EACH PARAMETER CHAIN AGAINST THE NEXT PARAMETER IN THE LIST
    
    histogram2(thisParam,nextParam, 50,...
        'DisplayStyle','tile',...
        'EdgeColor','none',...
        'Normalization','count'...
        );
    
    %Axes Limits
    xlim([min(thisParam), max(thisParam)]);
    ylim([min(nextParam), max(nextParam)]); 
    
    axis square;
    grid off;
    
%% SETUP AXES
    set(gca,...
        'Layer', 'top',...
        'FontSize',11,...
        'FontName', 'Arial',...
        'FontWeight', 'bold',...
        'LineWidth', 1.5,...
        'XMinorTick', 'off',...
        'YMinorTick', 'off',...
        'TickDir', 'both',...
        'TickLength', [0.02, 0.01]...
        );  
    
%% Add AXES LABELS DEPENDING ON THE SUBPLOT POSITION

         if i<=1 && j>=nParams
              set(gca);
              xlabel(axes_label_names{i});
              ylabel(axes_label_names{j});
         elseif j>=nParams
             set(gca,'YTicklabel', []);
             xlabel(axes_label_names{i});
         elseif i<=1;
             set(gca,'XTicklabel', []);
             ylabel(axes_label_names{j});
         else
             set(gca,'XTicklabel', [], 'YTicklabel', []);
         end  
end





