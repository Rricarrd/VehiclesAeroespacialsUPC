function fN = fNCalc(g, nnodes,M,NeumannDOF, dir)
% Fucntion that calculates the gravity body force from the mass matrix (M)
%
% Args:
%
% Returns:
%

if dir == 'X'
    accDOF = sparse( g*[1,0,0,0,0,0]); % Gravity acceleration vector applied to the Y axis (X,Y,Z,R1,R2,R3)
elseif dir == 'Y'
    accDOF = sparse( g*[0,1,0,0,0,0]); % Gravity acceleration vector applied to the Y axis (X,Y,Z,R1,R2,R3)
elseif dir == 'Z'
    accDOF = sparse( g*[0,0,1,0,0,0]); % Gravity acceleration vector applied to the Y axis (X,Y,Z,R1,R2,R3)
end

acc = repmat(accDOF, 1, nnodes); % Assemply of the vector for each node

fN = acc * M; % Force calculation

fN = fN(NeumannDOF)';
end

