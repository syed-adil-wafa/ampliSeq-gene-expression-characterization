%% Characterization of gene expression profiles of human iPSCs and iPSC-derived neurons using AmpliSeq

close all % close all figures
clear % clear workspace
clc % clear command window
rng(1); % seed random number generator for reproducibility

%% Ion AmpliSeq configurations

% Define the # of samples in each Ion AmpliSeq preparation
ampliseq_config.number_of_samples = 8;

% Define the gene expression start and end rows in data files
ampliseq_config.start_row = 2;
ampliseq_config.end_row = 20813;

% Define the gene name and NCBI name
ampliseq_config.gene_name_col = 1; % gene name column
ampliseq_config.ncbi_name_col = 12; % NCBI name column

% Define the sample start and end columns
ampliseq_config.sample_start_col = 3;
ampliseq_config.sample_end_col = ...
    ampliseq_config.sample_start_col + ampliseq_config.number_of_samples - 1;

%% Sample label information

% Define columns for each cell type
ampliseq_config.cell_cols = {...
    'IonXpress_051',... % control iPSCs
    'IonXpress_052',... % patient iPSCs
    'IonXpress_048',... % control neurons
    'IonXpress_049'}; % patient neurons

% Define labels for each cell type
ampliseq_config.labels_condition = {...
    'control iPSCs',... % control iPSCs
    'patient iPSCs',... % patient iPSCs
    'control neurons',... % control neurons
    'patient neurons'}; % patient neurons

%% Data extraction

cd data % change directory to extract data

% Retrieve all .xlsx files in the data folder
extraction.files = dir('*.xlsx');

% Create empty variables
extraction.raw_data = [];
extraction.col_names = [];

% Extract data from each file
for file = 1:size(extraction.files, 1)
    
    % Get the name of each individual file
    extraction.current_file = extraction.files(file).name;
    
    % Read data from each file
    [~, ~, extraction.file_data] = xlsread([pwd, ['/', extraction.current_file]]);
    extraction.raw_data = [extraction.raw_data, cell2mat(extraction.file_data(...
        ampliseq_config.start_row:size(extraction.file_data, 1),...
        ampliseq_config.sample_start_col:ampliseq_config.sample_end_col))];
    
    % Get gene and NCBI names of all 20,812 genes
    if file == size(extraction.files, 1)
        % Get gene names
        extraction.gene_names = extraction.file_data(...
            ampliseq_config.start_row:size(extraction.file_data, 1),...
            ampliseq_config.gene_name_col);
        % Get NCBI names
        extraction.ncbi_names = extraction.file_data(...
            ampliseq_config.start_row:size(extraction.file_data, 1),...
            ampliseq_config.ncbi_name_col);
    end
    
    % Get the column names for each file
    extraction.col_names = [extraction.col_names, extraction.file_data(...
        1, ampliseq_config.sample_start_col:ampliseq_config.sample_end_col)];
    
end

%% Data pre-processing

% Suppress warnings associated with data analysis
w = warning('on', 'all');
id = w.identifier;
warning('off', id);

cd ../helper-functions % change directory to enable helper functions

% Pre-process gene name strings
extraction.gene_names = gene_name_preprocess(...
    extraction.gene_names, extraction.ncbi_names);

% Extract required columns
columns = [];
for i = 1:size(ampliseq_config.cell_cols, 2)
    for j = 1:size(extraction.col_names, 2)
        if contains(extraction.col_names(:, j), ampliseq_config.cell_cols(:, i)) == 1
            columns = [columns, j];
        end
    end
end

% Create a DataMatrix object of the extracted data
processing.dm = bioma.data.DataMatrix(...
    extraction.raw_data(:, columns),... % limit to required data
    'RowNames', extraction.gene_names,...
    'ColNames', ampliseq_config.labels_condition);

% Make a copy of the data and filter values for clustergram
% Filter out genes with absolute expression levels below 10th percentile
[~, processing.dm_copy] = genelowvalfilter(processing.dm);

% Convert values <=1.0 to 1.0 for downstream log10 transformation
for rows = 1:size(processing.dm_copy, 1)
    for columns = 1:size(processing.dm_copy, 2)
        if processing.dm_copy.(rows)(columns) <= 1
            processing.dm_copy.(rows)(columns) = 1;
        end
    end
end

% Transform the data using log10 transformation
processing.dm_copy = log10(processing.dm_copy);

% Filter out genes with variance below 10th percentile
[~, processing.dm_copy] = genevarfilter(processing.dm_copy);

% Filter out genes with <10rpm across samples
filter = [];
processing.dm_low_expression = processing.dm_copy' < 1;
for gene = 1:size(processing.dm_low_expression, 2)
    if sum(processing.dm_low_expression(:, gene)) == size(processing.dm_copy, 2)
        filter = [filter; gene];
    end
end
processing.dm_copy(filter, :) = [];

% Convert values <=1.0 to 1.0 for downstream log10 transformation
for rows = 1:size(processing.dm, 1)
    for columns = 1:size(processing.dm, 2)
        if processing.dm.(rows)(columns) <= 1
            processing.dm.(rows)(columns) = 1;
        end
    end
end

% Transform the data using log10 transformation
processing.dm = log10(processing.dm);

%% Data visualization

% Construct clustered dendrograms of gene expression data
plots.heatmap = dm_clustergram(processing.dm_copy);
addTitle(plots.heatmap, 'Dendrogram of differential cellular gene expression');

% Define gene columns for markers (subtract 1 to account for header)
genes.markers = [13348;... % POU5F1 (OCT4)
    11229;... % NANOG
    16142;... % SOX2
    6268;... % FUT4
    18012;... % TUBB3
    10247;... % MAP2
    11286;... % NCAM1
    11423;... % NES (neural progenitor marker)
    6492;... % GFAP
    17841] - 1; % TSC2

plots.gene_data = processing.dm.(genes.markers)(':'); % numeric data

% Get labels and preserve the order of labels
plots.gene_labels = categorical(cellstr(processing.dm.RowNames(genes.markers)));
plots.gene_labels = reordercats(plots.gene_labels,...
    {'GFAP',... % glial marker
    'TSC2',... % TSC2 gene
    'NES',... % neural progenitor marker
    'POU5F1', 'NANOG', 'SOX2', 'FUT4',... % pluripotency markers
    'TUBB3', 'MAP2', 'NCAM1'}); % neuronal markers

% Plot bar graph of markers and their expression across cells
figure; hold on; box on;
set(gca, 'TickDir', 'in');
title('Differential cellular expression of pluripotency and neuronal markers');
ylabel('log_{10}_ (value)');
bar(plots.gene_labels, plots.gene_data);
legend(ampliseq_config.labels_condition, 'location', 'northwest');

% Compute correlation across specific pairs of lines and plot scatter plots
plots.pairs = [1, 2; 3, 4; 1, 3; 2, 4];

figure;
suptitle('Scatter plots of control and patient iPSCs and neurons');
for pair = 1:size(plots.pairs, 1)
    
    % Compute correlation
    X = processing.dm.(':')(plots.pairs(pair, 1));
    Y = processing.dm.(':')(plots.pairs(pair, 2));
    rho = corr(X, Y, 'Type', 'Pearson');
    
    % Plot scatter plots and add least-squares line
    subplot(2, 2, pair);
    hold on; box on;
    set(gca, 'TickDir', 'in');
    scatter(X, Y);
    h1 = lsline;
    h1.Color = 'k';
    h1.LineWidth = 1;
    
    % Add title and labels
    title(strcat('Pearson''s \itr\rm\bf=', num2str(rho)));
    xlabel(strcat('log_{10}_ (value) of', {' '},...
        ampliseq_config.labels_condition(plots.pairs(pair, 1))));
    ylabel(strcat('log_{10}_ (value) of', {' '},...
        ampliseq_config.labels_condition(plots.pairs(pair, 2))));
    
end

%% Remove unneeded variables from workspace
remove_variables = {'cols',...
    'columns',...
    'file',...
    'filter',...
    'gene_data',...
    'gene',...
    'h1',...
    'i',...
    'id',...
    'j',...
    'pair',...
    'rho',...
    'rows',...
    'w',...
    'X',...
    'Y'};
clear(remove_variables{1, :})
clear 'remove_variables'