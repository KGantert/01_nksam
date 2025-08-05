% Loop over different models
for tt = 1:4
    % Load default parameters
    load('parameters_updated.mat');

    % CALL DYNARE
    % ---------------------------------------------------------------------
    if tt == 1
        % Call Dynare NK model
        eval(['dynare nk_textbook.dyn -Dreplic_number=' num2str(replicnumber) ...
            ' -Dsimul_replic_number=' num2str(simreplicnumber) ' -Ddrop_number=' num2str(dropnumber) ...
	        ' -Dapprox_order=' num2str(approxorder) ' -Dirf_periods=' num2str(irfperiod) ...
	        ' -Dsim_periods=' num2str(simperiod) ' noclearall;']);
    elseif tt == 2
        % Call Dynare NK-SaM model
        eval(['dynare nk_sam.dyn -Dreplic_number=' num2str(replicnumber) ...
            ' -Dsimul_replic_number=' num2str(simreplicnumber) ' -Ddrop_number=' num2str(dropnumber) ...
	        ' -Dapprox_order=' num2str(approxorder) ' -Dirf_periods=' num2str(irfperiod) ...
	        ' -Dsim_periods=' num2str(simperiod) ' noclearall;']);
    elseif tt == 3
        % Call Dynare NK housework model
        eval(['dynare nk_housework.dyn -Dreplic_number=' num2str(replicnumber) ...
            ' -Dsimul_replic_number=' num2str(simreplicnumber) ' -Ddrop_number=' num2str(dropnumber) ...
	        ' -Dapprox_order=' num2str(approxorder) ' -Dirf_periods=' num2str(irfperiod) ...
	        ' -Dsim_periods=' num2str(simperiod) ' noclearall;']);
    elseif tt == 4
        % Call Dynare time allocation model
        eval(['dynare nk_time.dyn -Dreplic_number=' num2str(replicnumber) ...
            ' -Dsimul_replic_number=' num2str(simreplicnumber) ' -Ddrop_number=' num2str(dropnumber) ...
	        ' -Dapprox_order=' num2str(approxorder) ' -Dirf_periods=' num2str(irfperiod) ...
	        ' -Dsim_periods=' num2str(simperiod) ' noclearall;']);
    end

    % RETRIEVE SIMULATED DATA
    % ---------------------------------------------------------------------
    for ii = 1:length(var_sim_select)
        var_ind_SAM(ii)     = find(cellstr(M_.endo_names)==var_sim_select(ii));
    end
    % Simulated data
    y_sim_COMP = transpose(oo_.endo_simul(var_ind_SAM,:));

    % CALCULATE IRFs
    % ---------------------------------------------------------------------
    for ii = 1:length(M_.exo_names)
        for hh = 1:length(var_sim_select)
            eval(['irf_COMP_' char(M_.exo_names(ii)) '.' char(var_sim_select(hh)) '(:,tt)' ... 
                    '= transpose(oo_.irfs.' char(var_sim_select(hh)) '_' char(M_.exo_names(ii)) ').*100;']);
        end
    end

    % CALCULATE SIMULATED SECOND MOMENTS
    % ---------------------------------------------------------------------
    % Detrending
    [~,y_cyc_COMP]   = one_sided_hp_filter(y_sim_COMP(dropnumber+1:end,:),1600);
    % Relative Standard Deviations
    for hh = 1:size(y_cyc_COMP,2)
        for jj = 1:size(y_cyc_COMP,2)
            relstd_COMP(hh,jj,tt)   = std(y_cyc_COMP(:,hh)) ./ std(y_cyc_COMP(:,jj));
        end
    end
    % Contemporaneous Correlations
    corr_COMP(:,:,tt)    =  corrcoef(y_cyc_COMP);

    % Clear structures from the previous loop
    clearvars M_ oo_;

end

%% ------------------------------------------------------------------------
% IRF FIGURE FOR TFP SHOCK 
figure('Name','IRFs to a TFP Shock')
% Set up tiled layout depending on number of variables
tiledlayout(ceil(length(var_select)./4), 4);
% Loop IRF figures over all selected variables
for ii = 1:length(var_select)
    % Create next tile in the loop
    nexttile;
    % Plot all four different model IRFs in one graph
    hold on;
    eval(['area(irf_COMP_eA.' char(var_select(ii)) '(:,1),"LineStyle","none","FaceColor",[0.7 0.7 0.7])']);
    eval(['plot(irf_COMP_eA.' char(var_select(ii)) '(:,2),"r-","LineWidth",1)']);
    eval(['plot(irf_COMP_eA.' char(var_select(ii)) '(:,3),"b-","LineWidth",1)']);
    eval(['plot(irf_COMP_eA.' char(var_select(ii)) '(:,4),"k-","LineWidth",1)']);
    yline(0,'r');
    hold off;
    axis tight;
    title(char(title_select(ii)));
    ylabel('%-Dev.')
    xlabel('Quarters');  
end
% Print legend with the gamma values for the different models
leg = legend("NK w/o Home Production", "NK-SaM w/o Home Production", "NK Model", "NK-SaM Model", ...
                'Location', 'southoutside', 'Orientation', 'horizontal', 'NumColumns', 2);
leg.Layout.Tile = 'South';
%title(leg,'Model Type')
% Save figure
set(gcf, 'Units', 'centimeters');
set(gcf, 'Position', [0 0 width height]);
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperSize', [width height]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 width height]);
fontsize(gcf, 7,"points")
print(gcf, 'figures/irf_robust_comparison_tfp.png', '-dpng', '-vector');

%% ------------------------------------------------------------------------
% IRF FIGURE FOR MONETARY POLICY SHOCK 
figure('Name','IRFs to a Monetary Policy Shock')
% Set up tiled layout depending on number of variables
tiledlayout(ceil(length(var_select)./4), 4);
% Loop IRF figures over all selected variables
for ii = 1:length(var_select)
    % Create next tile in the loop
    nexttile;
    % Plot all four different model IRFs in one graph
    hold on;
    eval(['area(irf_COMP_eM.' char(var_select(ii)) '(:,1),"LineStyle","none","FaceColor",[0.7 0.7 0.7])']);
    eval(['plot(irf_COMP_eM.' char(var_select(ii)) '(:,2),"r","LineWidth",1)']);
    eval(['plot(irf_COMP_eM.' char(var_select(ii)) '(:,3),"b","LineWidth",1)']);
    eval(['plot(irf_COMP_eM.' char(var_select(ii)) '(:,4),"k","LineWidth",1)']);
    yline(0,'r');
    hold off;
    axis tight;
    title(char(title_select(ii)));
    ylabel('%-Dev.')
    xlabel('Quarters');  
end
% Print legend with the gamma values for the different models
leg = legend("NK w/o Home Production", "NK-SaM w/o Home Production", "NK Model", "NK-SaM Model", ...
                'Location', 'southoutside', 'Orientation', 'horizontal', 'NumColumns', 2);
leg.Layout.Tile = 'South';
%title(leg,'Model Type')
% Save figure
set(gcf, 'Units', 'centimeters');
set(gcf, 'Position', [0 0 width height]);
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperSize', [width height]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 width height]);
fontsize(gcf, 7,"points")
print(gcf, 'figures/irf_robust_comparison_policy.png', '-dpng', '-vector');

%% ------------------------------------------------------------------------
% IRF FIGURE FOR SEARCH EFFORT SHOCK 
figure('Name','IRFs to a Search Effort Shock')
% Set up tiled layout depending on number of variables
tiledlayout(ceil(length(var_select)./4), 4);
% Loop IRF figures over all selected variables
for ii = 1:length(var_select)
    % Create next tile in the loop
    nexttile;
    % Plot all four different model IRFs in one graph
    hold on;
    eval(['area(irf_COMP_eD.' char(var_select(ii)) '(:,1),"LineStyle","none","FaceColor",[0.7 0.7 0.7])']);
    eval(['plot(irf_COMP_eD.' char(var_select(ii)) '(:,2),"r","LineWidth",1)']);
    eval(['plot(irf_COMP_eD.' char(var_select(ii)) '(:,3),"b","LineWidth",1)']);
    eval(['plot(irf_COMP_eD.' char(var_select(ii)) '(:,4),"k","LineWidth",1)']);
    yline(0,'r');
    hold off;
    axis tight;
    title(char(title_select(ii)));
    ylabel('%-Dev.')
    xlabel('Quarters');  
end
% Print legend with the gamma values for the different models
leg = legend("NK w/o Home Production", "NK-SaM w/o Home Production", "NK Model", "NK-SaM Model", ...
                'Location', 'southoutside', 'Orientation', 'horizontal', 'NumColumns', 2);
leg.Layout.Tile = 'South';
%title(leg,'Model Type')
% Save figure
set(gcf, 'Units', 'centimeters');
set(gcf, 'Position', [0 0 width height]);
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperSize', [width height]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 width height]);
fontsize(gcf, 7,"points")
print(gcf, 'figures/irf_robust_comparison_search.png', '-dpng', '-vector');

%% ------------------------------------------------------------------------
% IRF FIGURE FOR MATCHING EFFICIENCY SHOCK 
figure('Name','IRFs to a Matching Efficiency Shock')
% Set up tiled layout depending on number of variables
tiledlayout(ceil(length(var_select)./4), 4);
% Loop IRF figures over all selected variables
for ii = 1:length(var_select)
    % Create next tile in the loop
    nexttile;
    % Plot all four different model IRFs in one graph
    hold on;
    eval(['area(irf_COMP_eT.' char(var_select(ii)) '(:,1),"LineStyle","none","FaceColor",[0.7 0.7 0.7])']);
    eval(['plot(irf_COMP_eT.' char(var_select(ii)) '(:,2),"r","LineWidth",1)']);
    eval(['plot(irf_COMP_eT.' char(var_select(ii)) '(:,3),"b","LineWidth",1)']);
    eval(['plot(irf_COMP_eT.' char(var_select(ii)) '(:,4),"k","LineWidth",1)']);
    yline(0,'r');
    hold off;
    axis tight;
    title(char(title_select(ii)));
    ylabel('%-Dev.')
    xlabel('Quarters');  
end
% Print legend with the gamma values for the different models
leg = legend("NK w/o Home Production", "NK-SaM w/o Home Production", "NK Model", "NK-SaM Model", ...
                'Location', 'southoutside', 'Orientation', 'horizontal', 'NumColumns', 2);
leg.Layout.Tile = 'South';
%title(leg,'Model Type')
% Save figure
set(gcf, 'Units', 'centimeters');
set(gcf, 'Position', [0 0 width height]);
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperSize', [width height]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 width height]);
fontsize(gcf, 7,"points")
print(gcf, 'figures/irf_robust_comparison_efficiency.png', '-dpng', '-vector');

%% ------------------------------------------------------------------------
% IRF FIGURE FOR EIS SHOCK 
figure('Name','IRFs to an EIS Shock')
% Set up tiled layout depending on number of variables
tiledlayout(ceil(length(var_select)./4), 4);
% Loop IRF figures over all selected variables
for ii = 1:length(var_select)
    % Create next tile in the loop
    nexttile;
    % Plot all four different model IRFs in one graph
    hold on;
    eval(['area(irf_COMP_eP.' char(var_select(ii)) '(:,1),"LineStyle","none","FaceColor",[0.7 0.7 0.7])']);
    eval(['plot(irf_COMP_eP.' char(var_select(ii)) '(:,2),"r","LineWidth",1)']);
    eval(['plot(irf_COMP_eP.' char(var_select(ii)) '(:,3),"b","LineWidth",1)']);
    eval(['plot(irf_COMP_eP.' char(var_select(ii)) '(:,4),"k","LineWidth",1)']);
    yline(0,'r');
    hold off;
    axis tight;
    title(char(title_select(ii)));
    ylabel('%-Dev.')
    xlabel('Quarters');  
end
% Print legend with the gamma values for the different models
leg = legend("NK w/o Home Production", "NK-SaM w/o Home Production", "NK Model", "NK-SaM Model", ...
                'Location', 'southoutside', 'Orientation', 'horizontal', 'NumColumns', 2);
leg.Layout.Tile = 'South';
%title(leg,'Model Type')
% Save figure
set(gcf, 'Units', 'centimeters');
set(gcf, 'Position', [0 0 width height]);
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperSize', [width height]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 width height]);
fontsize(gcf, 7,"points")
print(gcf, 'figures/irf_robust_comparison_eis.png', '-dpng', '-vector');