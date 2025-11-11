function plot_multipliers(full, reg, dataSpec, stats_reg,cfg)

T_B = stats_reg.T_B; 
horizons = cfg.modelSpec.horizons;
k = dataSpec.k;
if ~cfg.plot_irf, return; end

    fprintf(' Plot multipliers and CI...\n')
    fprintf('---------------------------------------------------------------\n');

    inityear = floor(min(dataSpec.years));
    initquar = (min(dataSpec.years)-inityear)/0.25;
    finayear = floor(max(dataSpec.years));
    finaquar = (max(dataSpec.years)-finayear)/0.25+1;
    lab_full_sample = sprintf('Proxy-SVAR, %d:Q%d-%d:Q%d',inityear,initquar,finayear,finaquar);
    
    breayear = floor(T_B);
    breaquar = (T_B-breayear)/0.25+1;
    init2ndreg = T_B + 0.25;
    inr2year = floor(init2ndreg) ;
    inr2quar = (init2ndreg-inr2year)/0.25+1;

    lab_r1 = sprintf('$1^{st}$ regime,  %d:Q%d-%d:Q%d',inityear,initquar,breayear,breaquar);
    lab_r2 = sprintf('$2^{nd}$ regime,  %d:Q%d-%d:Q%d',inr2year,inr2quar,finayear,finaquar);

    x = 0:horizons;
    irf_font_size = 20; irf_font_size_legend = 12; irf_line_width = 2;

    set(groot,'defaultAxesTickLabelInterpreter','latex');
    fig = figure(7); tiledlayout(k,1);

    % Tax multipliers
    nexttile;
    plot(x, full.Tax_Multiplier, 'k-', 'LineWidth', irf_line_width); hold on;
    plot(x, full.Tax_Multiplier_ub, 'k-.','HandleVisibility','off'); 
    plot(x, full.Tax_Multiplier_lb, 'k-.','HandleVisibility','off');
    p1 = plot(x, reg{1}.Tax_Multiplier, 'r-', 'LineWidth', irf_line_width); 
    fill([x fliplr(x)], [reg{1}.Tax_Multiplier_ub' fliplr(reg{1}.Tax_Multiplier_lb')], 'r', 'linestyle','none', 'facealpha',0.35,'HandleVisibility','off');
    p2 = plot(x, reg{2}.Tax_Multiplier, 'b-', 'LineWidth', irf_line_width);  
    fill([x fliplr(x)], [reg{2}.Tax_Multiplier_ub' fliplr(reg{2}.Tax_Multiplier_lb')], 'b', 'linestyle','none', 'facealpha',0.35,'HandleVisibility','off');
    yline(0);
    title('Tax multiplier','Interpreter','latex'); grid on; set(gca,'FontSize',irf_font_size); axis tight;
    legend({lab_full_sample, lab_r1, lab_r2}, 'Interpreter','latex','FontSize',irf_font_size_legend,'Location','best');

    % G multipliers if k>1
    if k>1
        nexttile;
        plot(x, full.G_Multiplier, 'k-', 'LineWidth', irf_line_width); hold on;
        plot(x, full.G_Multiplier_ub, 'k-','HandleVisibility','off'); 
        plot(x, full.G_Multiplier_lb, 'k-','HandleVisibility','off');
        plot(x, reg{1}.G_Multiplier, 'r-', 'LineWidth', irf_line_width);
        fill([x fliplr(x)], [reg{1}.G_Multiplier_ub' fliplr(reg{1}.G_Multiplier_lb')], 'r', 'linestyle','none', 'facealpha',0.35);
        plot(x, reg{2}.G_Multiplier, 'b-', 'LineWidth', irf_line_width,'HandleVisibility','off');
        fill([x fliplr(x)], [reg{2}.G_Multiplier_ub' fliplr(reg{2}.G_Multiplier_lb')], 'b', 'linestyle','none', 'facealpha',0.35,'HandleVisibility','off');
        yline(0);
        title('Spending multiplier','Interpreter','latex'); grid on; set(gca,'FontSize',irf_font_size); axis tight;
    end

    xlabel('Horizon','Interpreter','latex');
    if exist('setmyfig_out2','file') == 2
        setmyfig_out2(fig);
    end
    if ~isfolder('results/figures'), mkdir('results/figures'); end
    print(fig, fullfile('results','figures',[cfg.serial_id '_figure_1']), '-depsc');
    saveas(fig, fullfile('results','figures',[cfg.serial_id '_figure_1.pdf']));
    print(fig, fullfile('results','figures',[cfg.serial_id '_figure_1']), '-dpng', '-r300');

    fprintf('Figure with multipliers is saved in results/figures/%s\n',['id' cfg.serial_id '_figure_1.pdf']);



end
