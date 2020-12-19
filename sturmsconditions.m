function conds = sturmsconditions(tabsize, numroots, selections, vals)
% Find the necessary and sufficient conditions of a sign 
% table that permit exactly 'n' number of distinct roots.
%
% SYNTAX_____________________________________________________
%   conds = sturmsconditions(tsize, n)
%   conds = sturmsconditions(__, sel, vals)
%
% DESCRIPTION________________________________________________
%   conds = sturmsconditions(tsize, n) returns all possible
%       conditions 'conds' of a sign table with size 'tsize',
%       corresponding to 'n' number of roots. The conditions
%       are in a logically simplified, symbolic form in the
%       variable 't'. The subscript of 't' indicates the cell
%       position in the table.
%
%   conds = sturmsconditions(__, sel, vals) presets the cell
%       of the sign table at the positions selected by the
%       matrix 'sel', to the sign values 'vals'.
%
% INPUT ARGUMENTS____________________________________________
%   tsize:  vector
%       Size of the sign table of interest. The number of 
%       rows should be 2. E.g. dim = [2 x] where x is any
%       positive integer.
%
%   n:      number
%       The desired number of roots of the polynomial.
%
%   sel:    matrix (optional arg.)
%       Selection matrix where size(sel) == 'tsize'. Each 1 
%       selects the corresponding cell and 0 ignores it.
%       E.g. sel = [1 0 1; 0 1 0] selects the cells (1,1),
%       (1,3) and (2,1).
%
%   vals:   vector (required if 'sel' is given)
%       Sign values of the selected cells. The vector should
%       given in column-major format. E.g. given the
%       'sel' above, vals = [1 -1 0] assigns sign(t11) = 1,
%       sign(t21) = -1 and sign(13) = 0.
%
% OUTPUT_____________________________________________________
%   cond:   logical symbolic relation
%       Simplified logical expression indicating the possible
%       sign values of (some) cells.
%

% Since the last term of the Sturm's sequence is always a constant, its
% limits at the boundaries of the interval are equal and should be either
% positive or negative. Even if it was a zero term, it wouldn't have any
% effect on the number of sign change.
lastcol = [1 -1; 1 -1];
% The other columns can be any 2x1 combination of -1, 0 and 1 in the most
% general case of a Sturm's table.
cols = combvec(1:-1:-1, 1:-1:-1);

% Get the all possible repetitive combinations of these columns for the
% 2x3 sign table without including the last column. At this step, each 
% possible table is represented as a column-major vector.
possibletab = combvec(cols, cols);
for i = 3:tabsize(2)
    % If the table is larger than 2x3, append the remaining columns,
    if i < tabsize(2)
        possibletab = combvec(possibletab, cols);
    else % append the last column.
        possibletab = combvec(possibletab, lastcol);
    end
end

if exist('selections', 'var')
    % Flatten the selection matrix to match the format of possible tables.
    % Then get the indices of the selected cells.
    selections = selections(:) == 1;
    % Flatten the sign values. It should be in column-major format.
    vals = vals(:);
    
    % Find the indices of the tables having the values of all selected
    % cells equal to the specified 'vals'.
    possibleindices = all(possibletab(selections, :) == vals)...
                    & ones(2*tabsize(2), 1);
    % Get the eligible tables.
    possibletab = possibletab(possibleindices);
end

% Reshape the column-major vectors into ordirary 2xM sign tables.
possibletab = reshape(possibletab, tabsize(1), tabsize(2), []);

% Smooth the table before the numerical differentiation. This procedure is
% analoguos to edge detection algorithm in an image.
smoothtab = sign(smoothdata(possibletab, 2, 'gaussian', [0 1]));
% Differentiate the smooth table through the rows. If there is a sign
% change from 1 to -1 or vice versa, the differentiated value becomes 2,
% otherwise it is 0. By counting the number of 2's in the rows, we find 
% the number of sign changes in each row. Subtracting these numbers, we
% find the number of distinct real roots of the polynomial.
signchanges = sum(abs(diff(smoothtab, 1, 2)) == 2, 2).*[1; -1];
% Find the indices of the tables having exactly the desired number of 
% roots and reshape the eligible tables.
possibleindices = (sum(signchanges, 1) == numroots) & (ones(tabsize));
possibletab = reshape(possibletab(possibleindices),...
                      tabsize(1), tabsize(2),[]);

% Logically simplify the possible tables into symbolic relations.
conds = simplify_conditions(possibletab, selections, vals);
end

function conditions = simplify_conditions(tabs, selections, vals)
% Generate a symbolic table with integer elements.
tab = sym('t', size(tabs(:,:,1)), 'integer');
% The last column of the table can be either positive or negative.
assumeAlso(tab(1, end) == 1 | tab(1, end) == -1)
% The other elements can also be zero.
assumeAlso(-1 <= tab(:, 1:end-1) <= 1)
% Moreover, the cells of the last column should be equal.
assumeAlso(tab(1, end) == tab(2, end))

% Select the desired cells for value assignments.
selections = selections(:) == 1;
vals = vals(:);
coltab = tab(:);
% Override the assumptions on the selected cells. There should be no
% possibility other than the assigned values.
assume(coltab(selections, :) == vals)

conditions = false;
for i = 1:size(tabs, 3)
    % Get a possible table and assign its sign values to the symbolic
    % table. Then flatten the table for easier iteration in the next loop.
    condition = tab == tabs(:,:,i);
    condition = condition(:);
    andedliterals = true;
    for j = 1:length(condition)
        % All conditions in the table should be satisfied together.
        % Logically AND each condition.
        andedliterals = andedliterals & condition(j);
    end
    % Each table is valid under different sign configurations and only one
    % of them is possible for a specific model. Logically OR them.
    conditions = simplify(conditions | andedliterals, 500);
end
% Simplify all conditions to remove the redundancies and reduce the number
% of conditions.
conditions = simplify(conditions, 1000);
end