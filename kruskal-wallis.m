% -*- coding: utf-8 -*-
% Created on Fri Apr 12 02:27:40 2024
% @author: GBAGUIDI Mannondé D.

function correction_factor = rank_correction(rm)
    % Calculates the rank correction factor for a list of ranks.
    %
    % Args:
    %   rm: A vector of ranks.
    %
    % Returns:
    %   The rank correction factor, a scalar between 0 and 1. A value closer to 1
    %   indicates less correction is needed.
    %
    % Example usage:
    %   ranks = [1, 2, 2, 3, 1];
    %   correction_factor = rank_correction(ranks);
    %   disp(correction_factor);  % Output: 0.9

    sv = sort(rm);  % Sort the input vector of ranks
    
    % Create a histogram of rank counts
    count_values = unique(sv);  % Unique rank values
    count = histc(sv, count_values);  % Count occurrences of each rank value

    n = numel(sv);  % Get the number of elements

    corr = 0.0;  % Initialize correction factor

    % Calculate correction for tied ranks
    for i = 1:numel(count_values)
        if count(i) > 1
            corr = corr + (count(i)^3 - count(i));
        end
    end

    if n < 2
        correction_factor = 1.0;  % If less than 2 elements, correction factor is 1
    else
        correction_factor = 1.0 - corr / (n^3 - n);  % Calculate rank correction factor
    end
end


% Definition of the KW_test function, which takes a pandas DataFrame as a parameter.
function [p_value, H] = KW_test(data, interpretation, alpha)
    % Vérifier si interpretation n'est pas fourni
    if nargin < 2
        interpretation = true; % Valeur par défaut
    end
    
    % Vérifier si alpha n'est pas fourni
    if nargin < 3
        alpha = 0.05; % Valeur par défaut
    end
    % Check validity of interpretation
if ~(islogical(interpretation) || isscalar(interpretation))
    error('interpretation must be a logical scalar (true or false)');
end

% Check validity of alpha
if alpha > 0 && alpha < 1
    % Valid alpha value
else
    error('alpha must be a float number between 0 and 1, e.g., 0.05');
end

    % Extract variable names (group) from the DataFrame.
    listeV = data.Properties.VariableNames;
    % Determine the number of groups or variables in the DataFrame and store it in 'p'.
    p = length(listeV);
    if p < 2
    error('Need at least two groups in stats.kruskal()');
    end
    % Initialize an empty list 'n' to store the number of observations per group.
    n = [];
    % Loop through each variable in the DataFrame.
    for i = 1:p
        % Select the i-th column (variable) from the data.
        d = data(:, i);

        % Convert the selected column to a cell array.
        cel = table2cell(d);

        % Vertically concatenate the cell array into a single array.
        M = vertcat(cel{:});

        % Remove NaN values from the concatenated array.
        Mcl = M(~isnan(M));

        % Determine the number of non-NaN elements.
        [x, ~] = size(Mcl);

        % Append the number of non-NaN elements to the 'n' list.
        n = [n, x];
    end
         % Check for at least one data point in each group
    for i = 1:p
    if n(i) == 0
        error('Need at least one data in each group');
    end
    end


    % Calculate the total number of observations ('N') by summing the 'n' list.
    N = sum(n);

    % Convert the entire DataFrame to a cell array.
    cel = table2cell(data);

    % Vertically concatenate the cell array into a single array.
    M = vertcat(cel{:});

    % Remove NaN values from the concatenated array.
    Mcl = M(~isnan(M));

    % Compute the tied ranks for the non-NaN elements ('Rg').
    Rg = tiedrank(Mcl(:), 'average');

    % Initialize empty lists and variables for storing intermediate results.
    groupe = [];
    indi = 1;
    R = [];

    % Loop through each group and calculate the sum of ranks ('R') for each group.
    for i = 1:p
        groupe = Rg(indi:indi + n(i) - 1);
        R = [R, sum(groupe)];
        indi = indi + n(i);
    end

    % Calculate the sum of squares of ranks ('S') weighted by the group sizes.
    ScR = R.^2;
    S = sum(ScR./n);
    % Set the output format to 'long'.
    format long;

    % Compute the factor for the H statistic.
    facto = 12 / (N * (N + 1));

    % Compute the H statistic.
    H = facto * S - 3 * (N + 1);
    H = H/rank_correction(Rg);

    % Calculate the degrees of freedom for the chi-squared distribution.
    df = p - 1;

    % Compute the p-value using the chi-squared cumulative distribution function.
    p_value = 1 - chi2cdf(H, df);

    % Uncomment the following lines to display the results.
    if interpretation
     disp(['                      Kruskal-Wallis test']);
     disp(['                      ___________________']);
     disp(['                      ___________________']);
     disp([' H = ' + string(H) + ',        df = ' + string(df) + ',       p-value = ' + string(p_value) ]); 
     if (p_value <= alpha)
        disp(' Alternative hypothesis: True')
        disp(' So the medians are not all equal. ')
     else  
        disp(' Null hypothesis: True ')
        disp(' So the medians are all equal.')
     end
    end
end
