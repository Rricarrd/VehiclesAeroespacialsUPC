function [X, f_values] = freqResponse(values,f_max,nnodes,NeumannDOF,MODES,M,MNN,KNN,damping,FREQ)
%TEST
% Frequencies
f_values = linspace(0.1,f_max*2,values);

% Forces
accDOF = sparse( 1*[1,0,0,0,0,0]); % Gravity acceleration vector applied to the Y axis (X,Y,Z,R1,R2,R3)
acc = repmat(accDOF, 1, nnodes); % Assemply of the vector for each node
fN_freq = acc * M; % Force calculation
fN_freq = fN_freq(NeumannDOF)';

F_xi = MODES'*fN_freq;

m = diag(MODES'*MNN*MODES);
k = diag(MODES'*KNN*MODES);

xi = zeros(size(MODES,2),length(f_values));

for j = 1:length(FREQ)
    b = 2.*m.*FREQ.*damping;
    w = 2*pi*f_values;

    for i = 1:length(f_values)
        xi(j,i) = F_xi(j)./(-w(i)^2.*m(j)+k(j)+1i.*b(j).*w(i));
    
    end
end

% Modal method
X = MODES*xi;
end

