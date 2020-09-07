% Plot clustergram for DataMatrix object
function [dm_clustergram] = dm_clustergram(datamatrix)

% Get the numeric data from the DataMatrix object
dm = (datamatrix.(':')(':'));
row_names = cellstr(datamatrix.RowNames);
column_names = datamatrix.ColNames;

% Plot clustergram
dm_clustergram = clustergram(dm,...
    'Standardize', 'Row',...
    'RowLabels', row_names,...
    'ColumnLabels', column_names,...
    'RowPdist', 'correlation',...
    'ColumnPdist', 'correlation',...
    'Colormap', 'redbluecmap',...
    'DisplayRatio', [1/6, 1/5]);

% Set the font size to 6 for axis labels
set(0, 'ShowHiddenHandles', 'on');
all_handles = get(0, 'Children');
h = findall(all_handles, 'Tag', 'HeatMapAxes');
set(h, 'FontSize', 6);

end