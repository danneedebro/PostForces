function output = readStrfile(filename,minTime,maxTime,tsamp,maxPlotrec)
%
% readStrfile   ver 1.0
%
% Script som läser in RELAP5 stripfil (plotvariabler i textformat)
%
%
% syntax: readStripFileGUI(filname,minTime,maxTime,tsamp,plotrecCnt)
% 
% Input
%       filename        = sökväg till .str-fil (string)
%       minTime         = starttid (double)
%       maxTime         = sluttid (double)
%       tsamp           = samplingsfrekvens (integer)
%       maxPlotrec      = Antal plotrec att läsa (integer)
%
%
%
% Daniel Edebro
%
% ver 1.0
% Nytt script för att användas av 'postForces.m'
%

% Indatakontroll


returnOutput = 1;
plotRecCnt = 0;


fid = fopen(filename);

% Define cells to store data in
data = cell(0);
    
% loop through file in search for number of plot-variables
% indicated by number after the string 'plotinf'
tmpString = '';


firstLine = 1;
run = 1;
rad = 1;
sampInd = tsamp - 1;
tic
dataArray = [];
fileCount = 0;
blockSize = 2000;
randomString = sprintf('%04d%s',round(rand(1)*1000),datestr(now,'MMSS'));
% time = 0;


% ------- start while-loop ---------
while run
    tline = fgetl(fid);  % Läs in rad

    % Om sista raden
    if ~ischar(tline) || plotRecCnt > maxPlotrec
        tline = 'plotrec'; 
        run=0;
        lastRow = rad;
    end
    
    
    % ------- start if-block 1 ---------
    % Om 'plotinf' läses (plotinformation, integer)
    if ~isempty(strfind(tline,'plotinf'))
        readPlotinf = strread(tline,'%s');
        tmpString = '';
       
    elseif ~isempty(strfind(tline,'plotalf'))
        tmpString = tline;
        
    % Om 'plotnum' (variabelnummer, integer) läses ska tidigare inläst data lagras
    % som 'plotalf' (plotnamn, sträng)
    elseif ~isempty(strfind(tline,'plotnum')) 
        readPlotalf = strread(tmpString,'%s');
        N = size(readPlotalf,1);
        
        for i = 1:N, data{1,i} = readPlotalf{i,1}; end   % läs in 'plotalf' till rad 1 i 'data'
            tmpString = tline;   % reset 'tmpString' 
            dataArray = zeros(blockSize,size(data,2)-1);
            
    % Om 'plotrec' läses (datapunkter).
    elseif ~isempty(strfind(tline,'plotrec'))
        % Om första 'plotrec' ska tidigare inläst data lagras som 'plotnum'
        plotRecCnt = plotRecCnt + 1;
            
        % ------- start if-block 2 ---------
        if firstLine == 1
            readPlotnum = strread(tmpString,'%s');
            N = size(readPlotnum,1);
            for i = 1:N, data{2,i} = readPlotnum{i,1}; end   % läs in 'plotnum' till rad 2 i 'data'
            tmpString = tline;   % reset 'tmpString' 
            fprintf('Läser (%d st plotvariabler)         \n',N-2);
            firstLine = 0;

            % Om 'plotrec' läses. Processa tidigare inläst datarad
        else
            sampInd = sampInd + 1;

            readPlotrec = strread(tmpString,'%s');
            N = size(readPlotrec,1);
            

            time = str2double(readPlotrec{2});
            %'disp(sprintf('sampInd=%f, time=%f, minTime=%f, maxTime=%f',sampInd,time,minTime,maxTime))
                
            if isempty(time), disp('error'); end

            
            
            % ------- start if-block 3 ---------
            % Spara bara data om tiden är inom aktuellt område.
            if time >= minTime && time <= maxTime  && sampInd == tsamp
                sampInd = 0;
                for i = 2:N
                    dataArray(rad,i-1) = str2double(readPlotrec{i,1});
                end
                rad = rad + 1;
                if rad == blockSize + 1
                    fileCount = fileCount + 1;
                    save(sprintf('tmpFileName_%s_%04d.mat',randomString,fileCount),'dataArray');
                    dataArray = zeros(blockSize,size(data,2)-1);
                    rad = 1;
                    % Plotta progress-bar. 
                    fprintf('\b\b\b\b\b\b\b\b\b\b.t=%7.2f\n',time);
                end
                 
                tmpString = tline;  % reset 'tmpString' 
            elseif time > maxTime
                run = 0;
            else
                   
                tmpString = tline;  % reset 'tmpString'
                %'sampInd = sampInd - 1;
                if sampInd > tsamp, sampInd = 0; end
            end
            % ------- slut if-block 3 ---------
            
        end
        % ------- slut if-block 2 ---------
        
    % Om blank rad    
    elseif isempty(tline)
        tmpString = tmpString;  
    else
        tmpString = [tmpString,tline];  % utöka 'tmpString'
    end % while
    % ------- slut if-block 1 ---------
end
% ------- slut while-loop ---------

fclose(fid);
    
fileCount = fileCount + 1;
save(sprintf('tmpFileName_%s_%04d.mat',randomString,fileCount),'dataArray');
    
lastRow = rad;
        
    
% Redovisa tidsåtgång
tid = toc;
timeHour = floor(tid/3600);
timeMinute = floor( (tid - timeHour*3600)/60 );
timeSecond = floor( (tid - timeHour*3600 - timeMinute*60) );
fprintf('\nTidsåtgång (hh:mm:ss): %02d:%02d:%02d\n',timeHour,timeMinute,timeSecond);

    
antalRader = fileCount*blockSize-(blockSize-lastRow) - 1;
antalCols = size(data,2);
data{2+antalRader,antalCols} = [];
    
% Loopa igenom och läs in alla temporära filer och lagra dessa i 'data'
for i = 1:fileCount
    load(sprintf('tmpFileName_%s_%04d.mat',randomString,i),'dataArray');
    delete(sprintf('tmpFileName_%s_%04d.mat',randomString,i));
    M = 2+(i-1)*blockSize;
    if i < fileCount
        maxSize = blockSize;
    else
        maxSize = antalRader - blockSize*(fileCount-1);
    end
    for j = 1:maxSize
        for k = 1:size(dataArray,2)
            data{M+j,1+k} = dataArray(j,k);
        end
    end
end


% Loopa igenom 'data' och leta efter dubbla tidssteg
if size(data,1)>10
    time = data{3,2};
    doubleCnt = 0;
    failed = 0;
    for i = 4:size(data,1)-1
        time_prev = time;
        time = data{i,2};
        time_next = data{i+1,2};
    
        if time == time_prev
            doubleCnt = doubleCnt + 1;
            data{i,2} = (time_prev + time_next)/2;
            if time_prev == time_next, failed = failed + 1; end
        end
    end
    if doubleCnt > 0
        fprintf('\nWarning: %d dubbla tidssteg funna och ersatta (misslyckades = %d)\n',doubleCnt,failed);
    end
end


    
% returnera output om strängen 'output' given som argument. 
if returnOutput == 1
    output = data;
    clear data
end
        



