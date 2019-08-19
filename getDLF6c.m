function [frekvens,DLF] = getDLF6c(timeVector,forceVector,damping,cutoff)
% Funktion som givet en tidsvektor och en kraftvektor returnerar ett
% responsspektra.
%
% INPUT
% timeVector   = Mx1-vektor med tidsvektorn
% forceVector  = Mx1-vektor med kraftvektorn
% damping      = Dämpningskoeff. c/c_critical = damping. Ex 0.05 för 5%
%                dämpning
% cutoff       = Brytfrekvens. Om högre än vad som samplingteoremet medger
%                påverkar denna parameter inte
% OUTPUT
% frekvens     = Nx1-vektor med frekvenser (Hz) där N är den minsta värdet
%                av 'cutoff' eller f_max=1/(2*dt)
% DLF          = Nx1-vektor med smax_transient/smax_statisk
%
%

global dtVector

frekvens = 0;
DLF = 0;

if size(timeVector) ~= size(forceVector)
    disp('Tidsvektor och kraftvektor olika storlek');
    return
end

tmpDiffVec = diff(timeVector);
tmpMaxDt = max(tmpDiffVec);
tmpMinDt = min(tmpDiffVec);
tmpMeanDt = mean(tmpDiffVec);
dt = tmpMeanDt;    % FIXA SENARE - MER ROBUST
maxFreq = round(1/(2*dt));

if cutoff < maxFreq, maxFreq = cutoff; end


if tmpMinDt <= 0
    disp('Tidsvektor inte stigande');
    return
end


% Frekvensspann som ska kontrolleras
frekvens=1:1:maxFreq;
wn=frekvens*2*pi;
DLF = zeros(size(wn));

% Initialdata för DLF

k=10e3;                         % Fjäderkonstant
mn=k./wn.^2;                   % beräknar systemets massa (egenfrekv: w=rot(k/m))
smax=max(abs(forceVector))/k;  % Förskjutning



% Förfina kraftvektorn för att undvika numeriska problem
dtVector = [dt,dt/2,dt/4];
lastSplit = 4;
timeVectorRefined = cell(5,1);
forceVectorRefined = cell(5,1);
for i = 1:length(dtVector)
    tmp1 = timeVector(1):dtVector(i):timeVector(length(forceVector));
    timeVectorRefined{i} = tmp1;
    forceVectorRefined{i} = interp1(timeVector,forceVector,tmp1);
end



%%%%%     LOOPar varje frekvens     %%%%%
 for i = 1:length(wn)  
    
    DLFtry_0 = -9999;   % Startgissning för DLF sätts
    
    % Loopa igenom olika förfiningar av kraftvektorn. När konvergens nåtts
    % avslutas denna loop
    for ii = 1:length(dtVector)
        maxVal = 0; 
        x = 0;
        xp = 0;
        
        forceVec = forceVectorRefined{ii};
        
        %%%%%     BERÄKNAR DYNAMISK RESPONS     %%%%%
        for j = 1:length(forceVec)    % Tidssteg 
               
            Ft=forceVec(j);
        
            % Med dämpning
            xpp=(Ft-k*x-damping*2*mn(i)*wn(i)*xp)/mn(i);
        
            % Beräknar rörelseekvationen
            xp = xp + xpp*dtVector(ii);
            x = x + xp*dtVector(ii);

            % Om absolutbeloppet av x är större än tidigare värden - ersätt
            if abs(x) > maxVal, maxVal = abs(x); end 
        end     % j-loop
       
        % Kollar om lösning konvergerad
        DLFtry_1 = maxVal/smax;
        err = abs(DLFtry_0 - DLFtry_1)/abs(DLFtry_0);
        if err < 0.001
%             fprintf('Frekvens %d hittades efter %d iterationer (err=%1.4f)\n',i,ii,err);
            break; 
        else
            if ii == 22
                fprintf('\n(Lösning ej konvergerad efter 22 iterationer, f=%d av %d, err=%1.4f)--> ODE45 -->\n',i,length(wn),err); 
                break;
            % Om lösning ej konvergerad, utöka dtVector    
            elseif ii == length(dtVector)-1
%                 fprintf('\n(Lösning ej konvergerad, f=%d av %d, err=%1.4f)--> Utökar tidsvektor från längd %d till %d (dt/%d)',i,length(wn),err,length(dtVector),length(dtVector)+1,lastSplit + 2); 
                lastSplit = lastSplit + 2;
                dtVector = [dtVector,dt/lastSplit];
                tmp1 = timeVector(1):dtVector(ii+2):timeVector(length(forceVector));
                timeVectorRefined{i} = tmp1;
                forceVectorRefined{ii+2} = interp1(timeVector,forceVector,tmp1);
                
                
%                 Tspan = [timeVector(1) timeVector(length(timeVector))]; % Solve from t=1 to t=5
%                 options = odeset('RelTol',1e-3);
%                 [T Y] = ode45(@(t,y) springEq1(t,y,timeVectorRefined{ii},forceVec,mn(i),wn(i),damping),Tspan,[0 0],options); % Solve ODE
%                 fprintf('OK');
%                 maxVal = max(abs(Y(:,1)));
            end
            
            DLFtry_0 = DLFtry_1;
        end

    end   % ii-loop
    
    DLF(i)=maxVal/smax;
end                         % Frekvens 




function dy = springEq1(t,y,timeVec,forceVec,mn,wn,damping)

Ft=interp1(timeVec,forceVec,t);

dy = zeros(2,1);
dy(1) = y(2);
dy(2) = Ft/mn-2*damping*wn*y(2)-wn^2*y(1);
