function x_SpecEnvY = modifySpectralEnvelope(x,y,parameter)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% check parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin < 3
    parameter = [];
end
if nargin<2
    error('Please specify input data x and y.');
end

if ~isfield(parameter,'anaHop')
    parameter.anaHop = 64;
end
if ~isfield(parameter,'win')
    parameter.win = win(1024,1); % sin window
end
if ~isfield(parameter,'filterLength')
    parameter.filterLength = 24;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% some pre calculations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
anaHop = parameter.anaHop;
filLen = parameter.filterLength;
numOfChan = size(x,2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% modify the spectral envelopes channel wise
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x_SpecEnvY = zeros(size(x));            % Initialize output

for c = 1 : numOfChan                   % loop over channels
xC = x(:,c);
yC = y(:,c);

% compute the STFT
parStft.anaHop = anaHop;
parStft.win = parameter.win;
X = stft(xC,parStft);
Y = stft(yC,parStft);

% compute spectral envelopes
envX = compEnv(X,filLen);
envY = compEnv(Y,filLen);
% filLen = 510;
% envX = compTE(X,filLen);
% envY = compTE(Y,filLen);

X_specEnvY = X./envX.*envY;

% istft
parIstft.synHop = anaHop;
parIstft.win = parameter.win;
parIstft.zeroPad = 0;
parIstft.numOfIter = 1;
parIstft.origSigLen = size(x,1);
x_SpecEnvY(:,c) = istft(X_specEnvY,parIstft);

end

end

function env = compEnv(X,filLen)

    kern = win(filLen,2); % Hann window
    env = conv2(abs(X),kern,'same');

    % Scale the envelope such that the largest value of the envelope
    % coincides with the largest value of the spectrum (this is a
    % heuristic). We therfore first normalize the envelope such that the
    % largest value is 1 and afterwards multiply with the largest value of
    % the respective spectral frame.

    env = env ./ (eps + repmat(max(env),size(env,1),1)); % normalization
    env = env .* repmat(max(abs(X)),size(abs(X),1),1); % scaling

    % avoid values close to zero

    env(env<10^-2) = 10^-2;

end

function env = compTE(X,filLen)
filLen = 512 - 20;

% filLen = round(1.66 * filLen);

% numOfBins = size(X,1);
% sizeOfZeroPadding = round(1.66*numOfBins) - numOfBins + mod(round(1.66*numOfBins - numOfBins),2);

numOfIter = 100;

A = log(abs(X));
V = -Inf(size(X));
cInd = floor(size(X,1)/2+1);
for i = 1 : numOfIter
    A = max(A,V);
    C = fft(A);
    C(cInd-filLen:cInd+filLen,:) = 0;
    V = real(ifft(C));
end

env = exp(A);

% %% log freq axis
% % log axis
% linAxis = 1:size(X,1);
% logAxis = exp((linAxis-1)/linAxis(end)*5)-1;
% logAxis = logAxis/logAxis(end)*(linAxis(end)-1)+1;
%
%
% A = log(abs(X));
% A = interp1(linAxis,A,logAxis);
% % A = abs(X);
% V = -Inf(size(X));
% cInd = floor(size(X,1)/2+1);
% for i = 1 : numOfIter
%     A = max(A,V);
%
% %     Xpadded = [A; zeros(sizeOfZeroPadding, size(X,2))];
% %     Xshifted = fftshift(Xpadded,1);
% %     Xwin = fftshift(repmat(hann(size(Xshifted,1)),1,size(X,2)) .* Xshifted,1);
% %     cInd = floor(size(Xwin,1)/2+1);
% %     C = fft(abs(Xwin));
%     C = fft(A);
%     C(cInd-filLen:cInd+filLen,:) = 0;
%     V = real(ifft(C));
% %     V = V(1:end-sizeOfZeroPadding,:);
% end
%
% A1 = interp1(logAxis,A,linAxis);
%
% filLen = 512 - 15;
%
% A = log(abs(X));
% % A = interp1(linAxis,A,logAxis);
% % A = abs(X);
% V = -Inf(size(X));
% cInd = floor(size(X,1)/2+1);
% for i = 1 : numOfIter
%     A = max(A,V);
%
% %     Xpadded = [A; zeros(sizeOfZeroPadding, size(X,2))];
% %     Xshifted = fftshift(Xpadded,1);
% %     Xwin = fftshift(repmat(hann(size(Xshifted,1)),1,size(X,2)) .* Xshifted,1);
% %     cInd = floor(size(Xwin,1)/2+1);
% %     C = fft(abs(Xwin));
%     C = fft(A);
%     C(cInd-filLen:cInd+filLen,:) = 0;
%     V = real(ifft(C));
% %     V = V(1:end-sizeOfZeroPadding,:);
% end
%
% A2 = A;
%
% A = max(A1,A2);
%
% env = exp(A);
%
% %%

end
