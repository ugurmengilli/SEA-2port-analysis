function [table, seq] = sturms(p, var, lims, simplification)
%STURMS Apply Sturm's Theorem given a univariate polynomial
%   to determine the number of real roots. Let the number of
%   sign changes in the first row be n and the number of sign
%   changes in the second row be m. Then abs(n - m) gives the 
%   number of real roots of the polynomial.
%
%   tab = STURMS(P, VAR, LIMS) returns sign table given the
%       polynomial P, variable VAR and the limits
%       lim = [low_lim, high_lim].
%
%   [tab, seq] = STURMS(P, VAR, LIMS) returns sign table and
%       the Sturm's chain/sequence for convenience.
%   
%INPUTS:
%   p: Symbolic univariate polynomial in variable var.
%
%   var: Free variable of the polynomial.
%
%   lims: The interval [low_lim, high_lim] to be evaluated
%       for the number of roots.
%
%   simplification: The number of simplification steps for 
%       each term in the Sturm's sequence. Default is 500.
%
%OUTPUT:
%   table: 2 x n sign table where n <= deg(p)+1
%
%   seq: n x 1 Sturm's chain whose limits are evaluated.
%
if ~exist('simplification', 'var')
    simplification = 500;
end
% First two polynomials in Sturm's chain/sequence:
seq(1) = p;
seq(2) = diff(p, var);

i = 3;
while true
    
    [~, remainder] = quorem(seq(i-2), seq(i-1), var);
    remainder = simplify(remainder, simplification);
    
    if remainder == 0
        % then the sequence ends with the previous term.
        last_term = seq(i-1);

        if polynomialDegree(last_term, var)
            % The last term is a polynomial of any order in variable
            % 'var', which indicates that the polynomial p has some
            % repeating roots. These result in a series of zeros in
            % the Sturm's table and prevent the sign check if one of
            % the limits in 'lim' is actually a repeating root.
            % Therefore, the sequence should be updated dividing
            % each polyinomial in the sequence by these roots.
            % Factor out the last term:
            %
            % e.g. last_term = c (x + 1)^2*(x^2 + x + 2)
            %      factors   = [c, x+1, x+1, x^2+x+2]
            %
            repeating_roots = factor(last_term, var);
            % Clear the constant term (if any!) since it's not a root.
            repeating_roots = repeating_roots(...
                polynomialDegree(repeating_roots, var) > 0);

            disp("Repeating root(s) of the sequence:")
            disp(repeating_roots)

            seq = quorem(seq, prod(repeating_roots), var);
        else
            break;
        end
        break;
    end

    seq(i) = -remainder;
    i = i + 1;
end

table = [sign(limit(seq, lims(1)));
         sign(limit(seq, lims(2)))];

seq = seq';
end
