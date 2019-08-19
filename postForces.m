function postForces(varargin)
%
%
% postForces.m
% Script som genererar DLF (Dynamisk LastFaktor) på tidshistoriekrafter givna i RELAP5s
% .str-format. En högsamplad .str-fil innehållande krafter läses in
% varefter denna samplas ned i olika steg. En sammanfattande plott
% generere
%
% Köra i GUI-läge
%
% 1. starta med att skriva 'postForces' i MATLAB
% 2. Lokalisera .str-fil att läsa in
% 3. Ange basnamn för .pipestress-fil att generera (vid avbryt skrivs inga filer)
% 4. Ange basnamn för .ps (postscript) att generera (vid avbryt genereras inga plottar)
% 5. Ange tmin, tmax, brytfrekvens, dämpkvot och kraftamplitud över vilken inga DLF:er kommer genereras
% 6. Ange nedsampling. Sampling varje tidssteg (1X) är default
% 7. Vänta...
%
% Köra i batchläge
%
% 1. Starta genom att skriva 'postForces('-str','WHPT1.str',...,'tmin',50,...,etc)
%    Om inte samtliga parametrar anges
%
% Parameter         Beskrivning
% -str              .str-fil innehållande krafter
% -ps               Basnamn för plottar (nedsamplingar döps  XXXX_2X.ps)
% -peps             Basnamn för kraftfiler (nedsamplingar döps XXX_2X.pipestress)
% tmin              Starttid i sekunder
% tmax              Sluttid i sekunder
% tsamp             Samplingsfrekvens som grundfilen ska läsas in med
% cutoff            Brytfrekvens
% damp              Dämpkvot (0.05 betyder 5% av kritisk dämpning)
% forceLimit        Gräns när DLF ska utföras (ingen DLF under denna kraftamplitud)
% sampfactor        Heltal som visar nedsamplingen (2 samplas vartannat steg).
%                   Kan anges flera gånger.
% writeXLS          1 om sammanställning ska skrivas. 0 annars
% makePipestress	0 om pipestressfil INTE ska genereras (slippa pop-up)
% plotFigures       0 om plottar INTE ska genereras (slippa pop-up)
% selectForces      1 om du vill välja vilka krafter som ska köras
% title             Titel på plott
% verbosity         0, 1 eller 2 beroende på hur mycket text som ska skrivas ut
%
%
%
% Kom ihåg: Snabba upp
%    - Kolla ett brutalt stort fall (0.1 ms, 200 krafter exempelvis. Även
%    på vanlig dator)
%    - Skriv ut spektrum för olika längd på lasten (på orginalfilen eller
%    en nedsamplad version. Kanske bättre på en nedsamplad version)
%
% Version 1.00      Daniel Edebro
% Nytt script
%
%
%
versionNumber = 1.00;


global Data


locateInputFile = 1;
locatePlotFile = 1;
locatePepsFile = 1;
selectForces = 0;
setPlotReduction = 1;
plotFigures = 1;
makePipestress = 1;
writeXLSsummary = 1;
lastPath = pwd;
verbosity = 0;
titleString = '';

tmin = [];
tmax = [];
tsamp = [];
cutoff = [];
damp = [];
discardForceAmp = [];

cutoffDef = 300;
dampDef = 0.05;
discardForceAmpDef = 100;
multVec = 1;
screening = 10;

% Parametrar som styr utseendet på plottarna
colorVector = [1,0,0;0,0,0;0,0,1;0,1,0;0.5,0.5,0];


fprintf('postForces version %1.2f\n\n',versionNumber);
fprintf('1. LÄSER ARGUMENT OCH KONTROLLERAR INDATA\n');


% Skriv ut inläst syntax (underlättar felsökning)
fprintf('   syntax: postForces(');
for i = 1:size(varargin,2);
    if ischar(varargin{i})
        tmpStr = ['''',num2str(varargin{i}),''''];
    else
        tmpStr = num2str(varargin{i});
    end
    
    if i < size(varargin,2);
        fprintf('%s,',tmpStr);
    else
        fprintf('%s',tmpStr);
    end
end
fprintf(')\n\n');


% Läs in argument i 'varargin'
fprintf('   argument:\n'); 
for i = 1:size(varargin,2);
    
    % om 'inputfilename' läses antas att nästa argument är filnamn
    if strcmpi(varargin{i},'-str')
        stripFilePath = '';
        stripFile = varargin{i+1};
        locateInputFile = 0;
        i = i + 1;%#ok
        fprintf('      .str-fil: %s\n',stripFile);
    % om filnamn för jämförelseplott        
    elseif strcmpi(varargin{i},'-ps')
        plotFilename = varargin{i+1};
        locatePlotFile = 0;
        i = i + 1;%#ok
        fprintf('      .ps-fil: %s\n',plotFilename);
    % plotta figurer        
    elseif strcmpi(varargin{i},'plotFigures')
        plotFigures = varargin{i+1};
        i = i + 1;%#ok
        fprintf('      plotFigures: %d\n',plotFigures);
    elseif strcmpi(varargin{i},'-peps')
        tmpPepsPath = varargin{i+1};
        
        % Parse:a sökvägen
        periods = strfind(tmpPepsPath,'\');
        if isempty(periods)
            pipestressfilePath = ''; 
            pipestressfile = tmpPepsPath;
        else
            pipestressfilePath = tmpPepsPath(1:periods(length(periods)));
            pipestressfile = tmpPepsPath(1+periods(length(periods)):length(tmpPepsPath));
        end
        locatePepsFile = 0;
        i = i + 1;%#ok
        fprintf('      .pipestress-fil: %s\n',pipestressfile);    
    % gör pipestressfiler
    elseif strcmpi(varargin{i},'makePipestress')
        makePipestress = varargin{i+1};
        i = i + 1;%#ok
        fprintf('      makePipestress: %d\n',makePipestress);
    % skriv sammanfattande XLS
    elseif strcmpi(varargin{i},'writeXLS')
        writeXLSsummary = varargin{i+1};
        i = i + 1;%#ok
        fprintf('      writeXLS: %d\n',writeXLSsummary);
    % tmin
    elseif strcmpi(varargin{i},'tmin')
        tmin = varargin{i+1};
        i = i + 1;%#ok
        fprintf('      tmin: %1.4f\n',tmin);
    % tmax
    elseif strcmpi(varargin{i},'tmax')
        tmax = varargin{i+1};
        i = i + 1;%#ok
        fprintf('      tmax: %1.4f\n',tmax);  
    % tsamp
    elseif strcmpi(varargin{i},'tsamp')
        tsamp = varargin{i+1};
        i = i + 1;%#ok
        fprintf('      tsamp: %1.4f\n',tsamp);      
    % cutoff
    elseif strcmpi(varargin{i},'cutoff')
        cutoff = varargin{i+1};
        i = i + 1;%#ok
        fprintf('      cutoff: %d\n',cutoff);
    elseif strcmpi(varargin{i},'title')
        titleString = varargin{i+1};
        i = i + 1;%#ok
        fprintf('      title: %s\n',titleString);
    % damp
    elseif strcmpi(varargin{i},'damp')
        damp = varargin{i+1};
        i = i + 1;%#ok
        fprintf('      cutoff: %1.4f\n',damp);        
    % forceLimit
    elseif strcmpi(varargin{i},'forceLimit')
        discardForceAmp = varargin{i+1};
        i = i + 1;%#ok
        fprintf('      forceLimit: %1.2f N\n',discardForceAmp);        
    % sampFactor
    elseif strcmpi(varargin{i},'sampFactor')
        multVec = [multVec,varargin{i+1}];%#ok
        setPlotReduction = 0;
        fprintf('      sampFactor: %d\n',varargin{i+1});
        i = i + 1;%#ok
    elseif strcmpi(varargin{i},'selectForces')
        selectForces = 1;
        fprintf('      selectForces: 1\n');
        
    % verbosity. 0 = minimalt med utskrifter, 1 = fler utskrifter, 
    % 2 = ännu fler
    elseif strcmpi(varargin{i},'verbosity')
        verbosity = varargin{i+1};
        fprintf('      verbosity: %d\n',varargin{i+1});
        i = i + 1;%#ok
    elseif strcmpi(varargin{i},'help')
        fprintf('Syntax: Display help\n');
        return
    end
end
fprintf('\n');


ElapsedTime = 0;


% Definera en temporär sökväg för att snabba upp skrivning av
% pipestressfiler när man jobbar mot en nätverksdisk
tempPath = 'c:\temp\';


% Om flaggan 'locateInputFile' är sann så ska fil att läsa specificeras
if locateInputFile == 1
    fprintf('   Manuellt angiven indatafil:\n'); 
    [stripFile, stripFilePath] = uigetfile({'*.str','RELAP5 Stripfil (*.str)';'*.*','Alla filer (*.*)'},'Öppna stripfil...');
    if stripFile == 0, return; end
    if ~strcmp(stripFilePath(length(stripFilePath):length(stripFilePath)),'\'), stripFilePath = [stripFilePath,'\']; end 
    lastPath = stripFilePath;
    fprintf('      .str-fil:\t\t%s\n\n',[stripFilePath,stripFile]);
end

% Om flaggan 'makePipestress' är sann så ska pipestressfilers sökväg och filnamn definieras
% Sker då parametern '-peps' inte anges som argument eller argumentet
% 'makePipestress',0 anges
if makePipestress == 1 && locatePepsFile == 1
        [pipestressfile, pipestressfilePath] = uiputfile({'*.pipestress','Pipestress load data (*.pipestress)';'*.*','Alla filer (*.*)'},'Pipestress: Spara som...',lastPath);
        if pipestressfile == 0 
            makePipestress=0; 
        else
            if ~strcmp(pipestressfilePath(length(pipestressfilePath):length(pipestressfilePath)),'\'), pipestressfilePath = [pipestressfilePath,'\']; end
            lastPath = pipestressfilePath;
        end
end

% Om flaggan 'makePipestress' är sann så ska postscriptfilens sökväg och filnamn definieras
% Sker då parametern '-ps' inte anges som argument eller argumentet
% 'plotFigures' = 0
if locatePlotFile == 1
    [plotfile, plotfilePath] = uiputfile({'*.ps','Post script file (*.ps)';'*.*','Alla filer (*.*)'},'Plot: Spara som...',lastPath);
    if plotfile == 0
        plotFigures=0; 
    else
        if ~strcmp(plotfilePath(length(plotfilePath):length(plotfilePath)),'\'), plotfilePath = [plotfilePath,'\']; end
        %lastPath = plotfilePath;
    end
    
    plotFilename = [plotfilePath,plotfile];
end



% Kontrollera input
fprintf('   Kontrollerar .str-fil...\n');
if exist([stripFilePath,stripFile],'file')
    fprintf('      finns...OK\n');
else
    fprintf('      finns...Error, saknas\n');
    return;
end

f=dir([stripFilePath,stripFile]);
fprintf('      storlek: %1.4f MB\n',f.bytes*1e-6);
[t1,t2] = getTimeVector([stripFilePath,stripFile]); 
tmin_file = t1;
tmax_file = t2;
fprintf('      t1=%1.4f s och t2=%1.4f s\n\n',t1,t2);

allParamsGiven = 0;
if ~isempty(tmin), t1 = tmin; end
if ~isempty(tmax), t2 = tmax; end

% Kontrollera om samtliga parametrar givna. Om så är fallet öppnas inte
% inputdialogen där tmin,tmax,brytfrekvens, etc anges.
if ~isempty(tmin) && ~isempty(tmax) && ~isempty(damp) && ~isempty(cutoff) && ~isempty(discardForceAmp) && ~isempty(tsamp), allParamsGiven = 1; end

if allParamsGiven == 0
    fprintf('   Läser manuellt angivna parametrar:\n\n');
    xInt=inputdlg({'tmin (s)','tmax (s)','tsamp','Cutoff (Hz)','Damping coeff','Min force (N)'},'Intervall:',1,{num2str(t1),num2str(t2),'1',num2str(cutoffDef),num2str(dampDef),num2str(discardForceAmpDef)});

    if isempty(xInt), return; end
    tmin = str2double(xInt{1});
    tmax = str2double(xInt{2});
    tsamp = str2double(xInt{3});
    cutoff = str2double(xInt{4});
    damp = str2double(xInt{5});
    discardForceAmp = str2double(xInt{6});
end

% Om ingen nedsampling angivits - visa dialogruta
if setPlotReduction == 1
    downSamp = {'2','3','4','5','6','7','8','9','10'};
    [s,v] = listdlg('PromptString','Välj nedsampling:','SelectionMode','multiple','ListString',downSamp);
    if v ~= 0
        multVec = [multVec,str2double(downSamp(s))];
    end
end




fprintf('   Kontrollerar angivna parametrar:\n');
if tmin < tmin_file, fprintf('      tmin < tmin_file\n'); else fprintf('      tmin: OK\n');end
if tmax > tmax_file, fprintf('      tmax > tmin_file\n'); else fprintf('      tmax: OK\n');end
fprintf('      Tidssampling: ');
for i = 1:length(multVec), fprintf('%dX dt0,',multVec(i)); end
fprintf('\n');
    
fprintf('\n');


% Om 'selectForces' lika med 1 visas listruta där samtliga krafter visas
% varefter användaren själv kan välja vilka som ska plottas.
if selectForces == 1
    rawData = readStrfile([stripFilePath,stripFile],tmin,tmax,1,15);
    plotvarList = cell(1,size(rawData,2)-2);
    for i = 3:size(rawData,2), plotvarList{1,i-2}=sprintf('%s-%s',rawData{1,i},rawData{2,i}); end
    [forcesToProcess,v] = listdlg('PromptString','Välj krafter att plotta:','SelectionMode','multiple','ListString',plotvarList);
    if v == 0, return; end
end






tic
fprintf('2. LÄSER .STR-FIL (%s)\n',[stripFilePath,stripFile]);
rawData = readStrfile([stripFilePath,stripFile],tmin,tmax,tsamp,10000000);
ElapsedTime=ElapsedTime+toc;

timeVec = cell2mat(rawData(3:size(rawData,1),2));
fprintf('\n');

tmin = timeVec(1);
tmax = timeVec(length(timeVec));
dt = mean(diff(timeVec));


% Om figurer ska plottas, parsa sökvägen för att konstruera filnamn för
% nedsamplingar utifrån basnamnet på plotfilen
%
% Filnamn.ps    - Sammanfattande plott
% Filnamn_1X.ps - Orginaldatans plott
% Filnamn_4X.ps - Nedsamplingens plott
if plotFigures == 1
    periods = strfind(plotFilename,'.');
    if isempty(periods)
        plotfilenamePrefix = plotFilename; 
        plotfilenameSuffix = '';
    else
        plotfilenamePrefix = plotFilename(1:periods(length(periods))-1);
        plotfilenameSuffix = plotFilename(periods(length(periods)):length(plotFilename));
    end
end






Data = struct([]);
Nmax = length(multVec);   % antal nedsamplingar


% Om pipestressfiler ska genereras - öppna dessa filer för output
% Dessa namnges utifrån basnamnet på samma sätt som för plottfilen
% (ex Filnamn.pipestress, Filnamn_1X.pipestress, etc)
% 
if makePipestress == 1
    periods = strfind(pipestressfile,'.');
    if isempty(periods)
        filenamePrefix = pipestressfile; 
        filenameSuffix = '';
    else
        filenamePrefix = pipestressfile(1:periods(length(periods))-1);
        filenameSuffix = pipestressfile(periods(length(periods)):length(pipestressfile));
    end
    
    % Loopa igenom 'multVec' och öppna pipestressfiler för input
    fidPeps = zeros(1,Nmax);
    for i = 1:Nmax
        fidPeps(i) = fopen(sprintf('%s%s_%dX%s',tempPath,filenamePrefix,multVec(i),filenameSuffix),'w');
        if fidPeps(i) == -1, fprintf('   Error, kan inte öppna ''%s%s_%dX%s'' för output\n',tempPath,filenamePrefix,multVec(i),filenameSuffix); makePipestress = 0; end 
    end
end



% Loopa igenom alla krafter
%
% Yttre loop - alla krafter
%      Inre loop - alla nedsamplingar
%           Generera DLF
%           Spara i 'Data'
%           Plotta
%           Skriv pipestressfiler
%           Skriv till den sammanfogade .ps-filen
%           Skriv till de individuella .ps-filerna (Filnamn_2X.ps)
%
fprintf('3. LOOPAR KRAFTER\n');

freqPlot = zeros(size(multVec));
timePlot = zeros(size(multVec));

if plotFigures == 1, newFig = figure(1); end


% Om inga krafter väljs och samtliga plottas - sätt vilka
if selectForces == 0
    totalNumberOfForces = size(rawData,2)-2;
    force1 = 1;   
    force2 = totalNumberOfForces;

    forcesToProcess = force1:force2;
end

currForcePos = 0;     % Räknare var kraften ska placeras i 'Data'
for i = 1:length(forcesToProcess)
    tic
    currForceInd = forcesToProcess(i);
    currForcePos = currForcePos + 1;
    currPlotvar = rawData{1,currForceInd+2};
    currPlotnumStr = rawData{2,currForceInd+2};
    currPlotnum = str2double(currPlotnumStr);
    
    fprintf('   Plot %d av %d: %s-%s',i,length(forcesToProcess),currPlotvar,currPlotnumStr);
    forceVec = cell2mat(rawData(3:size(rawData,1),currForceInd+2));
    
    maxF = max(abs(forceVec));
    
    % Om maximal kraftamplitud mindre än brytamplitud - beräkna ingen DLF
    if maxF < discardForceAmp
        fprintf(' (mindre än %1.2f N, skippar DLF)\n',discardForceAmp);
        calcDLF = 0; 
    else
        fprintf(' (kör DLF)\n');
        calcDLF = 1; 
    end
    
    % Loopa igenom samtliga nedsamplingar i 'multVec'
    for j = 1:Nmax
        tic
        timeVecN = tmin:dt*multVec(j):tmax;
        
        if verbosity > 0, fprintf('      nersampling %d: ',j-1);end
        forceVecTmp = forceVec;
        for k = 1:length(forceVec)-multVec(j)
            forceVecTmp(k) = mean(forceVec(k:k+multVec(j)-1));
        end
        
        forceVecN = interp1(timeVec,forceVecTmp,timeVecN);
        
        % Hämta frekvensvektor första gången så denna säkert ska finnas
        % tillgänglig även om DLF inte ska beräknas
        if i == 1
            tmpForce = sin(2*pi*15*timeVecN)+sin(2*pi*25*timeVecN);
            [frekvens,~] = getDLF6c(timeVecN,tmpForce,damp,cutoff);
            Data(j).frequency = frekvens;
        end
        
        % Om DLF ska beräknas
        if calcDLF == 1
            [frekvens,DLF] = getDLF6c(timeVecN,forceVecN,damp,cutoff);
            
            % Skala om lastfaktorn utifrån maximal kraftamplitud för
            % urspungskraften för att undvika fel om toppar samplas bort
            maxFN = max(abs(forceVecN));
            DLF = DLF*maxFN/maxF;
        else   % Om DLF inte ska beräknas, sätt värdet till 1
            frekvens = Data(j).frequency;
            DLF = ones(size(frekvens));
        end

        % Spara värden i 'Data'
        Data(j).maxVal(currForcePos) = max(forceVecN);
        Data(j).minVal(currForcePos) = min(forceVecN);
        Data(j).plotvar{currForcePos} = currPlotvar;
        Data(j).plotnum{currForcePos} = currPlotnumStr;
        Data(j).plotnum2(currForcePos) = currPlotnum;
        Data(j).DLF(currForcePos,:) = DLF;
        
        
        % 
        if plotFigures == 1
            plot1 = subplot(2,1,1);
            if j == 1, hold off; else hold on; end    % Rensar plott om första
            timePlot(j) = plot(timeVecN,forceVecN);
            title([currPlotvar,'-',currPlotnumStr]);
               
            set(gca,'XLim',[tmin,tmax]);
            % Om maxamplitud under gränsvärde - sätt y-axelns värde till
            % +/-gränsvärdet för att ge snyggare figurer (slippa brus kring
            % nollan)
            if calcDLF == 0, set(gca,'YLim',[-discardForceAmp,discardForceAmp]); end
            
           
            plot2 = subplot(2,1,2);
            if j == 1, hold off; else hold on; end    % Rensar plott om första
            freqPlot(j) = plot(frekvens,DLF);
            
            % Plotta felgränser för DLF (ex +- 20%)
            hold on
            if j == Nmax
                tmpFreq = Data(1).frequency;
                tmpDLF = Data(1).DLF(currForcePos,:);
                freqPlot(j+1) = plot(tmpFreq,tmpDLF*(100+screening)/100);
                freqPlot(j+2) = plot(tmpFreq,tmpDLF*(100-screening)/100);
            end
            
        end
        

        
        % Om pipestress-fil ska genereras
        if makePipestress == 1
            fprintf(fidPeps(j),'TFUN CV_%s\n',currPlotnumStr);
            for k = 1:length(timeVecN)
                fprintf(fidPeps(j),'%1.5f %1.3f\n',timeVecN(k),forceVecN(k));
            end
            fprintf(fidPeps(j),'-1.E11\n');
        end
        
        timeLag = toc;
        ElapsedTime=ElapsedTime+timeLag;
        if verbosity > 0, fprintf('  Klar efter %1.3f sekunder\n',timeLag);end
        
    end
    
    
    if plotFigures == 1
        
        
        
        legendText = cell(1,length(multVec));
        for j = 1:Nmax
            set(timePlot(j),'LineWidth',1);
            set(timePlot(j),'LineStyle','-');
            set(timePlot(j),'MarkerSize',3);
            set(timePlot(j),'Color',colorVector(j,:));
            if j>1, set(timePlot(j),'Visible','off'); end
            
            set(freqPlot(j),'Color',colorVector(j,:));
            set(freqPlot(j),'LineWidth',1);
            set(freqPlot(j),'LineStyle','-');
            legendText{j} = sprintf('t_s_a_m_p = %dX dt_0',multVec(j));
            
        end
        
        legendText{Nmax+1} = sprintf('t_s_a_m_p = %dX dt_0 +/- %1.0f%%',multVec(1),screening);
        %legendText{Nmax+2} = sprintf('dt/%d - 10%%',multVec(1));
        
        % Sätt plottegenskaper för felgränserna
        set(freqPlot(length(multVec)+1),'LineStyle','.');
        set(freqPlot(length(multVec)+2),'LineStyle','.');
        set(freqPlot(length(multVec)+1),'Color',colorVector(1,:));
        set(freqPlot(length(multVec)+2),'Color',colorVector(1,:));
        set(freqPlot(length(multVec)+1),'MarkerSize',3);
        set(freqPlot(length(multVec)+2),'MarkerSize',3);
        
        title(plot1,sprintf('%s-%s',Data(1).plotvar{currForcePos},Data(1).plotnum{currForcePos}));
        title(plot2,'DLF');
        plot1legend = legend(plot1,legendText{1},'Location','Best');%#ok
        plot2legend = legend(plot2,legendText,'Location','Best');
        xlabel(plot1,'Tid (s)');
        xlabel(plot2,'Frekvens (Hz)');
        ylabel(plot1,'Kraft (N)');
        ylabel(plot2,'DLF x_m_a_x/x_m_a_x_,_s_t_a_t_i_c (-)');
        
        % Släck DLF-plott om ingen beräknats
        if calcDLF == 0
            for k = 1:Nmax
                set(freqPlot(k),'Visible','off');
            end
            set(freqPlot(Nmax+1),'Visible','off');
            set(freqPlot(Nmax+2),'Visible','off');
            set(plot2legend,'Visible','off');
        end
        
        
        % Printa som en PostScript-fil. Ny bilder bifogas i slutet
        set(gcf,'PaperUnits','centimeters')
        set(gcf,'PaperType','a4')   
        set(gcf, 'PaperOrientation', 'landscape')
        set(gcf,'PaperPosition',[0.634517 0.634517 28.4084 19.715]);
        % set(gcf,'PaperPositionMode','auto');
        set(gcf, 'Visible', 'off')
        
        if i == 1
            annotation(newFig,'textbox',[0.1000 0.9168 0.7200 0.06077],'String',titleString,'EdgeColor',[0 0 0],'LineStyle','none','HorizontalAlignment','center','FontSize',14,'FontName','Arial');
        end
    
        % Om första grafen kör INTE 'append'
        if i == 1
            print('-dpsc','-loose','-r300',plotFilename)
        else
            print('-dpsc','-loose','-r300','-append',plotFilename)
        end
        
        % Skriv ut individuella plottar
        for j = 1:Nmax
            % Släck alla time-history- och DLF-plottar
            for k = 1:Nmax
                set(timePlot(k),'Visible','off');
                set(freqPlot(k),'Visible','off');
            end
            set(freqPlot(Nmax+1),'Visible','off');
            set(freqPlot(Nmax+2),'Visible','off');
            set(timePlot(j),'Visible','on');
            if calcDLF == 1, set(freqPlot(j),'Visible','on');set(freqPlot(Nmax+1),'Visible','on');set(freqPlot(Nmax+2),'Visible','on'); end
            if j > 1
                set(timePlot(j),'Color',colorVector(2,:));
                set(freqPlot(j),'Color',colorVector(2,:));
            end
            
            legend(timePlot(j),legendText{j},'Location','Best')
            if calcDLF == 1
                legend([freqPlot(j),freqPlot(Nmax+1)],{legendText{j},legendText{Nmax+1}},'Location','Best');
            else
                set(plot2legend, 'Visible', 'off');
            end
            
            
            plotFilenameInd = sprintf('%s_%dX%s',plotfilenamePrefix,multVec(j),plotfilenameSuffix);
            if i == 1
                print('-dpsc','-loose','-r300',plotFilenameInd)
            else
                print('-dpsc','-loose','-r300','-append',plotFilenameInd)
            end
            
        end
    
        % hold off;bar([1,2],[13,2],'FaceColor',[0,0,1]);hold on;bar([1,2],[-4,-6],'FaceColor',[1,0,0]);
    end
    
end



fprintf('\n4. AVSLUTAR\n');



% Stäng fid:ar och kopiera pipestressfiler från temporär mapp
if makePipestress == 1
    % Stäng fid:ar
    for i = 1:Nmax
        fclose(fidPeps(i));
        src = sprintf('%s%s_%dX%s',tempPath,filenamePrefix,multVec(i),filenameSuffix);
        trg = sprintf('%s%s_%dX%s',pipestressfilePath,filenamePrefix,multVec(i),filenameSuffix);
        fprintf('   Flyttar "%s" från temporär mapp\n      %s --> %s\n',sprintf('%s_%dX%s',filenamePrefix,multVec(i),filenameSuffix),src,trg);
        success = movefile(src,trg);
        if success == 0, fprintf('      Error, misslyckades att flytta pipestressfil från "%s" till "%s"\n',src,trg); end
    end
end

% Skriv sammanfattning
if writeXLSsummary == 1 && plotFigures == 1
    fid = fopen([plotfilenamePrefix,'.xls'],'w');
    if fid == -1
        fprintf('\n   Error, misslyckades att skriva till %s.xls\n      Skriver till _ver2.xls\n',plotfilenamePrefix); 
        fid = fopen([plotfilenamePrefix,'_ver2.xls'],'w');
        if fid == -1, fprintf('\n   Error, misslyckades att skriva till %s_ver2.xls\n',plotfilenamePrefix);  return; end
    end
    
    
    forceCnt = length(Data(1).plotvar);
    
    % Antal datapunkter och tidssteg
    fprintf(fid,'\n');
    fprintf(fid,'POSITIVA MAXAMPLITUDER\t');
    for i = 1:Nmax, fprintf(fid,'\t'); end
    fprintf(fid,'\tNEGATIVA MAXAMPLITUDER\n\t');
    for i = 1:length(multVec)
        fprintf(fid,'%dX dt0\t',multVec(i));
    end
    fprintf(fid,'\t');
    for i = 1:length(multVec)
        fprintf(fid,'%dX dt0\t',multVec(i));
    end
    fprintf(fid,'\n');
    
    for i = 1:forceCnt
        fprintf(fid,'%s-%d\t',Data(1).plotvar{i},Data(1).plotnum2(i));
        for j = 1:Nmax, fprintf(fid,'%f\t',Data(j).maxVal(i)); end 
        fprintf(fid,'\t');
        for j = 1:Nmax, fprintf(fid,'%f\t',Data(j).minVal(i)); end 
        fprintf(fid,'\n');
    end
    
    % Skriv DLF
    fprintf(fid,'\n');
    fprintf(fid,'DLF\n');
    fprintf(fid,'\t');
    for i = 1:Nmax, fprintf(fid,'%dX dt0\t',multVec(i)); end
    fprintf(fid,'\n');
    
    for i = 1:forceCnt
        [DLFmax,DLFmaxInd] = max(Data(1).DLF(i,:));
        if DLFmax == 1, tmpFlag = 1; else tmpFlag = 0; end
            
        if tmpFlag == 0
            fprintf(fid,'%s-%d\t%1.2f\t',Data(1).plotvar{i},Data(1).plotnum2(i),DLFmax);
        else
            fprintf(fid,'%s-%d\t\t',Data(1).plotvar{i},Data(1).plotnum2(i));
        end
 
        for j = 2:Nmax
            fMax = length(Data(j).DLF(i,:));
            if sum(Data(j).DLF(i,1:fMax)>Data(1).DLF(i,1:fMax)*(100-screening)/100) == fMax && sum(Data(j).DLF(i,1:fMax)<Data(1).DLF(i,1:fMax)*(100+screening)/100) == fMax && tmpFlag == 0
                fprintf(fid,'\t');
            elseif sum(Data(j).DLF(i,1:fMax)>Data(1).DLF(i,1:fMax)*(100-2*screening)/100) == fMax && sum(Data(j).DLF(i,1:fMax)<Data(1).DLF(i,1:fMax)*(100+2*screening)/100) == fMax && tmpFlag == 0
                fprintf(fid,'%1.0f %% < X < %1.0f %%\t',screening,2*screening);
            elseif tmpFlag == 1
                fprintf(fid,'\t');
            else
                fprintf(fid,'Check \t');
            end
        end
        fprintf(fid,'\t @ %d Hz\n',Data(1).frequency(1,DLFmaxInd));
    end
    
    fprintf(fid,'\n');
    fprintf(fid,'\n');
    
    for i = 1:Nmax
        fprintf(fid,'DLF för tsamp = %dX dt0\n',multVec(i));
        freqCnt = length(Data(i).frequency);
        
        fprintf(fid,'Fmax =\t');
        for k = 1:forceCnt, fprintf(fid,'%1.2f\t',max(abs(Data(1).maxVal(1,k)),abs(Data(1).minVal(1,k)))); end
        fprintf(fid,'\n');
        
        fprintf(fid,'Frekvens (Hz)\t');
        for k = 1:forceCnt, fprintf(fid,'%s-%d\t',Data(i).plotvar{1,k},Data(i).plotnum2(1,k)); end        
        fprintf(fid,'\n');
        
        for j = 1:freqCnt
            fprintf(fid,'%1.0f\t',Data(i).frequency(1,j));
            for k = 1:forceCnt
                fprintf(fid,'%1.2f\t',Data(i).DLF(k,j));
            end
            fprintf(fid,'\n');
        end
        fprintf(fid,'\n');
    end
    
    % Skriv ut DLF:er uttryckt i kraft (Newton)
    for i = 1:Nmax
        fprintf(fid,'DLF för tsamp = %dX dt0\n',multVec(i));
        freqCnt = length(Data(i).frequency);
        
        fprintf(fid,'Frekvens (Hz)\t');
        for k = 1:forceCnt, fprintf(fid,'%s-%d\t',Data(i).plotvar{1,k},Data(i).plotnum2(1,k)); end        
        fprintf(fid,'\n');
        
        for j = 1:freqCnt
            fprintf(fid,'%1.0f\t',Data(i).frequency(1,j));
            for k = 1:forceCnt
                Fmax = max(abs(Data(1).maxVal(1,k)),abs(Data(1).minVal(1,k)));
                fprintf(fid,'%1.2f\t',Data(i).DLF(k,j)*Fmax);
            end
            fprintf(fid,'\n');
        end
        fprintf(fid,'\n');
    end
    fclose(fid);
end



% Skriv ut syntax för att lättare kunna upprepa samma körning
fprintf('Syntax\n');
fprintf('>postForces(');
fprintf('''-str'',''%s''',[stripFilePath,stripFile])

if plotFigures == 1
    fprintf(',''-ps'',''%s''',plotFilename)
else
    fprintf(',''plotFigures'',0');
end

if makePipestress == 1
    fprintf(',''-peps'',''%s''',[stripFilePath,stripFile])
else
    fprintf(',''makePipestress'',0');
end

if writeXLSsummary == 1, fprintf(',''writeXLS'',1'); else fprintf(',''writeXLS'',0'); end


fprintf(',''tmin'',%1.4f',tmin);
fprintf(',''tmax'',%1.4f',tmax);
fprintf(',''damp'',%1.4f',damp);
fprintf(',''cutoff'',%1.0f',cutoff);
fprintf(',''forceLimit'',%1.2f',discardForceAmp);
for i=2:Nmax, fprintf(',''sampFactor'',%d',multVec(i)); end
fprintf(')\n');


fprintf('  \n\nKlar efter totalt %1.3f sekunder\n',ElapsedTime);






function [t1,t2] = getTimeVector(filename)
% Funktion som hämtar första och sista tidsrecord från stripfil
%
% Input
% filename = filnamn att öppna
%
% Output
% t1 = första tidsrecord
% t2 = sista tidrecord
%
fid=fopen(filename);
 
t1=[];
t2=9999;


% Loop igenom stripfilen och hitta första respektiva sista instansen av
% 'plotrec' och kontrollera vilka tider dessa korresponderar mot (t1 och t2)
while 1
    tline = fgetl(fid);
    if ~ischar(tline), break, end
    
    % Om första 'plotrec' - sätt lässtart 50 000 bytes från slutet
    if ~isempty(strfind(tline,'plotrec')) && isempty(t1)
        res=textscan(tline, '%s');
        t1 = str2double(res{1}{2});
        fseek(fid,-50000,1);
    % Om andra 'plotrec' X antal bytes från slutet
    elseif ~isempty(strfind(tline,'plotrec'))
        res=textscan(tline, '%s');
        t2 = str2double(res{1}{2});
    end
 end
fclose(fid);
 


