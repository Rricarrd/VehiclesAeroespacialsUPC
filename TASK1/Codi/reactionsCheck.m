function [w_error,mass] = reactionsCheck(M,F,g,DirichlettDOF,dir)

    % Calculation of the mass in kg
    mass = sum(diag(M))/3;
    
    %Check if Dirichlett forces in the direction are equal to weight
    if dir == 'X'
        w_error = repmat([1 0 0 0 0 0],1,6)*F(DirichlettDOF)+mass*g;
    elseif dir == 'Y'
        w_error = repmat([0 1 0 0 0 0],1,6)*F(DirichlettDOF)+mass*g;
    elseif dir == 'Z'
        w_error = repmat([0 0 1 0 0 0],1,6)*F(DirichlettDOF)+mass*g;
    end

end

