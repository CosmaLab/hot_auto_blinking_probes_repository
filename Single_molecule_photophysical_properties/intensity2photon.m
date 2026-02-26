function photons = intensity2photon(intensity, params)
% INTENSITY2PHOTON Converts intensity values to photon counts
% 
% Input:
%   intensity - Intensity values in counts
%   params - Camera parameters, can be struct or array
%       If struct: params.gain, params.QE
%       If array: params(1) = gain, params(2) = QE
%
% Output:
%   photons - Estimated photon counts

% Check if params is a struct or array
if isstruct(params)
    % Extract parameters from struct
    gain = params.gain;        % electrons/ADU
    QE = params.QE;            % Quantum efficiency
else
    % Extract parameters from array
    gain = params(1);          % electrons/ADU
    QE = params(2);            % Quantum efficiency
end

% Convert intensity to photon counts
electrons = intensity * gain;  % Convert counts to electrons
photons = electrons / QE;      % Convert electrons to photons

end