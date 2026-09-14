% Function to count carbon atoms in a chemical formula
function n_C = countCarbons(formula)
    % Initialize the carbon atom counter
    n_C = 0;
    
    % Find carbon atoms in the formula using regular expressions
    tokens = regexp(formula, 'C(\d*)', 'tokens');
    
    % Iterate over the regular expression results
    for k = 1:length(tokens)
        % If a number is found, use that number; otherwise, assume 1
        if isempty(tokens{k}{1})
            n_C = n_C + 1;  % If there is no number, there is 1 carbon atom
        else
            n_C = n_C + str2double(tokens{k}{1});  % If there is a number, add the value
        end
    end
end