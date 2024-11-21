function create_LCModel_fit_figures()

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Creates publication-ready LCModel fit figures from directory of COORD outputs 
% 
%
% Inputs: LCModel COORD files to use in figure creation 
% 
% Function run prompts user to select root folder containing directory of
% COORD files for which to create fit figures
% 
% Outputs: .png, .eps, and potentially .pdf (see plot inputs) files
% displaying fitted data, baseline, basis functions, final model, and
% residual with corresponding legend labels 
% 
% Author: Kelley Swanberg (Lunds universitet, 2024) 
% swanberg@post.harvard.edu 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %% Define plot inputs
    space_value = 35; 
    baseline_plot_color = '#00D1D1'; 
    fit_plot_color = 'r'; 
    fit_data_color = 'k'; 
    plot_dimensions = [0 0 500 1280]; 
    save_as_png = 1; 
    save_as_eps = 1; 
    save_as_pdf = 0; 
    
    %% Prompt user to select root directory containing folder of COORD files for which to write .sh script
    current_folder = uigetdir(); 
    
    % Find COORD files for which to write script 
    list_of_cases_struct = dir(current_folder);
    list_of_cases_struct_clean = list_of_cases_struct(~ismember({list_of_cases_struct.name},{'.','..'}));
    list_of_cases = {list_of_cases_struct_clean.name}'; 
    % https://se.mathworks.com/matlabcentral/answers/431023-list-all-and-only-files-with-no-extension#answer_348046
    list_of_cases_notcoord = contains(list_of_cases, '.');
    list_of_cases(list_of_cases_notcoord) = []; 

    % Now we have the list of cases in memory as a vertical vector 
    num_coords = length(list_of_cases);  
    
    % Prepare table for outputs 
    % combined_array = zeros(num_coords, 5); 
    % combined_table = array2table(combined_array); 
    % row_names = {num_coords};
    
    %% For each COORD file in directory create a publication-ready fit figure
    for ii = 1:num_coords
    
        % Find and load data coords
        case_name = list_of_cases{ii}; 
    
        fid=fopen(case_name,'r');
        f=fread(fid,'*char')';
        fclose(fid);
    
        f_lines = splitlines(f); 
    
        % Find first line of data
        start_index = find(contains(f_lines,'points on ppm-axis'));
        end_index = find(contains(f_lines,'phased data points follow'));
        num_data_lines = end_index - start_index - 1; 
    
        % Find num_metabs
        metabs_line_first = find(contains(f_lines,'in following concentration table'));
        metabs_line_last = find(contains(f_lines,'in following misc. output table'));
        num_metabs = metabs_line_last - metabs_line_first - 2; 
    
        % Find metab names 
        metab_search_first = metabs_line_first+2; 
        metab_search_last = metabs_line_last-1;
        metab_name_list = {num_metabs};
    
        for jj = metab_search_first : metab_search_last
            metabstring = f_lines(jj); 
            metabstring_clean = regexprep(metabstring,'[+\s]',' ');
            metabstring_split = split(metabstring_clean); 
            metabstring_split_clean = metabstring_split(~cellfun('isempty',metabstring_split)); 
            metab_name_list_index = jj - metab_search_first + 1; 
            metab_name =  metabstring_split_clean(end); 
            metab_name_list{metab_name_list_index} = metab_name; 
        end
    
        % Remove repeat values from metab name list 
        metab_name_list_trans = metab_name_list'; 
        metab_name_list = table2cell(unique(cell2table(metab_name_list_trans), 'rows')); 
        
        num_metabs_to_plot = length(metab_name_list); 
        
        % Find and prepare fit ppm axis 
        ppm_axis_line_first = start_index+1; 
        ppm_axis_line_last = end_index-1; 
        ppm_axis_lines = f_lines(ppm_axis_line_first:ppm_axis_line_last); 
        
        ppm_axis_vector = coord_lines_to_vector(ppm_axis_lines); 
        
        % Find and prepare data coords 
        data_coords_line_first = end_index+1; 
        data_coords_line_last = data_coords_line_first + num_data_lines-1; 
        data_coords_lines = f_lines(data_coords_line_first:data_coords_line_last);
    
        data_coords_vector = coord_lines_to_vector(data_coords_lines); 
        
        % Find and load fit coords 
        fit_coords_line_first = data_coords_line_last+2; 
        fit_coords_line_last = fit_coords_line_first + num_data_lines-1; 
        fit_coords_lines = f_lines(fit_coords_line_first:fit_coords_line_last);
        
        fit_coords_vector = coord_lines_to_vector(fit_coords_lines); 
        
        % Find and load baseline coords 
        baseline_coords_line_first = find(contains(f_lines,'background values follow'))+1;
        baseline_coords_line_last = baseline_coords_line_first + num_data_lines-1; 
        baseline_coords_lines = f_lines(baseline_coords_line_first:baseline_coords_line_last);
        
        baseline_coords_vector = coord_lines_to_vector(baseline_coords_lines); 

        % Calculate residual coords 
        residual_coords_vector = data_coords_vector - fit_coords_vector; 
        
        num_points = length(ppm_axis_vector); 
        metab_coords_vector_array = zeros(num_metabs_to_plot, num_points); 
    
        % Calculate metabolite coords
        for kk=1:num_metabs_to_plot
            metab_coords_vector = zeros(1, num_points); 
            metab_coords_line_first = []; 
        
            string_to_find = sprintf(' %s   Conc.', metab_name_list{kk});
            metab_coords_line_first = find(startsWith(f_lines,string_to_find))+1;
        fit_coords_lines
            % There is almost certainly a better way using regular expressions
            % Three spaces
            if isempty(metab_coords_line_first)
                string_to_find = sprintf(' %s   Conc.', metab_name_list{kk});
                metab_coords_line_first = find(startsWith(f_lines,string_to_find))+1;
            end

            % Four spaces
            if isempty(metab_coords_line_first)
                string_to_find = sprintf(' %s    Conc.', metab_name_list{kk});
                metab_coords_line_first = find(startsWith(f_lines,string_to_find))+1;
            end
        
            % Five spaces
            if isempty(metab_coords_line_first)
                string_to_find = sprintf(' %s     Conc.', metab_name_list{kk});
                metab_coords_line_first = find(startsWith(f_lines,string_to_find))+1;
            end
        
            % Six spaces
            if isempty(metab_coords_line_first)
                string_to_find = sprintf(' %s      Conc.', metab_name_list{kk});
                metab_coords_line_first = find(startsWith(f_lines,string_to_find))+1;
            end
        
            % Seven spaces
            if isempty(metab_coords_line_first)
                string_to_find = sprintf(' %s       Conc.', metab_name_list{kk});
                metab_coords_line_first = find(startsWith(f_lines,string_to_find))+1;
            end
        
            if ~isempty(metab_coords_line_first)
                metab_coords_line_last = metab_coords_line_first + num_data_lines-1; 
                metab_coords_lines = f_lines(metab_coords_line_first:metab_coords_line_last);
                metab_coords_vector = coord_lines_to_vector(metab_coords_lines); 
            end
        
            metab_coords_vector_array(kk, :) = metab_coords_vector; 
        
        end 
    
        % Plot fit figure
        figure(1); 
        
        % Plot baseline 
        baseline_coords_vector_plot = baseline_coords_vector+space_value*(num_metabs_to_plot + 2);
        plot(ppm_axis_vector, baseline_coords_vector_plot,'LineWidth',1,'Color',baseline_plot_color,'DisplayName','base');
        hold on; 

        % Plot data
        data_coords_vector_plot = data_coords_vector+space_value*(num_metabs_to_plot + 3); 
        plot(ppm_axis_vector, data_coords_vector_plot,'LineWidth',1,'Color',fit_data_color,'DisplayName','data'); 
        
        % Plot fit
        fit_coords_vector_plot = fit_coords_vector+space_value*(num_metabs_to_plot + 3);
        plot(ppm_axis_vector, fit_coords_vector_plot,'LineWidth',1,'Color',fit_plot_color,'DisplayName','fit');

        % Plot metabs 
        cmap = cool(num_metabs_to_plot);
        for jj=1:num_metabs_to_plot
            metab_coords_vector_plot = metab_coords_vector_array(jj, :) + space_value*(num_metabs_to_plot + 1) - space_value*jj; 
            plot(ppm_axis_vector, metab_coords_vector_plot,'LineWidth',1,'Color',cmap(jj, :),'DisplayName',metab_name_list{jj});
        end
        
        % plot residual 
        plot(ppm_axis_vector, residual_coords_vector,'LineWidth',1,'Color','#808080','DisplayName','res');
        
        set(gca, 'XDir','reverse'); 
        hold off; 
        legend('Location','eastoutside'); 
        legend boxoff; 
        box off; 
        ylimits = ylim;
        ylim_lower = -2*space_value; 
        ylim_upper = ylimits(2); 
        set(gca,'ylim',[ylim_lower ylim_upper]);
        xlim("tight"); 
        xlabel('Chemical shift (ppm)'); 
        set(gca,'ytick',[]); 
        set(gcf,'Position', plot_dimensions); 
        set(gcf,'Color',[1 1 1]); 
    
        % Save as PNG
        if save_as_png == 1
            filename_png = sprintf('%s_fitfig.png', case_name); 
            saveas(gcf, filename_png); 
        end

        % Save as EPS
        if save_as_eps == 1
            filename_eps = sprintf('%s_fitfig.eps', case_name); 
            exportgraphics(gcf,filename_eps,'BackgroundColor','none','ContentType','vector');
        end

        % Save as PDF 
        if save_as_pdf == 1
            filename_pdf = sprintf('%s_fitfig.pdf', case_name); 
            exportgraphics(gcf,filename_pdf,'BackgroundColor','none','ContentType','vector');
        end
    
        fprintf('Exported figure %d of %d!\n', ii, num_coords); 
    end 

    close all; 
end

function [coord_vector] = coord_lines_to_vector(coord_lines)

    coord_lines_split = split(coord_lines(1:end-1)); 
    coord_lines_split_1d = reshape(coord_lines_split.',1,[]);

    % Remove empty cells: https://se.mathworks.com/matlabcentral/answers/
    % 209-how-do-i-remove-the-empty-cells-from-a-vector-of-cells#answer_253
    coord_lines_split_1d_clean = coord_lines_split_1d(~cellfun('isempty',coord_lines_split_1d))'; 

    coord_lines_split_last = split(coord_lines(end)); 
    coord_lines_split_last_clean = coord_lines_split_last(~cellfun('isempty',coord_lines_split_last)); 

    coord_vector =str2double([coord_lines_split_1d_clean; coord_lines_split_last_clean]); 

end