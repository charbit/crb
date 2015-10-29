%=============== settings.m ======================================
clear
% ================  IMPORTANT ====================================
% the different records form the sensors can not be aligned
% Here we cut tehm to the longest commun part and omit the rest.
%
%
%=================================================================

Fs_Hz = 20;

switch computer
    case 'GLNXA64'
        addpath /dvlscratch/SHI/users/charbit/ProjectIMS2015b/myjob/1TaskOnSensors/textes/6distConjointHMSC/fullprocess/ZZtoolbox/00pierrick/
    otherwise
        addpath /Users/maurice/etudes/ctbto/allJOBs2015/myjob/1TaskOnSensors/textes/6distConjointHMSC/fullprocess/ZZtoolbox/00pierrick/
end
%========== matlab format files are saved into the following directory
stationnumber         = 37;
directorydatafromIDC  = sprintf('../../../../AAdataI%i/',stationnumber);
%=== sensor locations
cdeloadX = sprintf('load(''../../../../sensorlocation/I%i.mat'');',stationnumber);
eval(cdeloadX)
Msensors               = size(xsensors_m.coordinates,1);
%====================================================================

%=== temporary files
temporary_gparse_dir = 'tempfiles/';
if not(exist(temporary_gparse_dir,'dir'))
    rmdir(temporary_gparse_dir,'s')
    mkdir(temporary_gparse_dir)
else
    commandrmparse = sprintf('!rm %s.*',temporary_gparse_dir);
    eval(commandrmparse);
end
if exist('gparse.wfdisc','file')
    !rm gparse*.*;
end
%========= source of data
data_source = 'testbed_archive';
user        = 'charbit';
password    = 'sqlmomo';
% channel     = '(''BDF'',''BDF'',''LWS'',''LWD'',''LKO'')';
channel     = '(''BDF'')';

yearstart   =  '2014';
monthstart  =  '09';
HMSstart    = '00:00:10';
yearend     =  '2014';
monthend    =  '09';
HMSend      = '23:50:10';
records     = cell(10,1);

ihclist =  xsensors_m.name;
for daystart_num    = 16:30% %1:2:25
    if daystart_num<10
        daystart    = ['0' num2str(daystart_num)];
        if daystart_num==9
            dayend  = '10';
        else
            dayend  = ['0' num2str(daystart_num+1)];
        end
    else
        daystart    = num2str(daystart_num);
        dayend      = num2str(daystart_num+1);
    end
    dayend = daystart;
    for indhc=1:Msensors
        elementtype = xsensors_m.name{indhc}(4);
        ihc = str2double(xsensors_m.name{indhc}(5));
        stations    = sprintf(' (''I%i%s%i'') ',stationnumber,elementtype, ihc);
        
        %=== clean temporary files
        commandclean = sprintf('!rm %s/*.*',temporary_gparse_dir);
        eval(commandclean)
        %=== extract data from the database
        h_starttime   = sprintf('%s/%s/%s %s',yearstart,monthstart,daystart, HMSstart);
        h_endtime     = sprintf('%s/%s/%s %s',yearend,monthend,dayend, HMSend);
        [~,starttime] = unix(['h2e ',h_starttime,' ofmt="%#"']);
        [~,endtime]   = unix(['h2e ',h_endtime,' ofmt="%#"']);
        starttime     = str2double(starttime);
        endtime       = str2double(endtime);
        wlength       = endtime-starttime;
        %============= pipeline.m =========================================
        %==== Write the query
        fid           = fopen(sprintf('%sgparse_temp.par',temporary_gparse_dir),'w');
        fprintf(fid, 'open data_source=%s user=%s password=%s\n', data_source, user, password);
        fprintf(fid, '%s\n',['query wfdisc select * from sel3.wfdisc where sta in ', ...
            stations, 'and chan in ', channel,'  and time between ',num2str(starttime),...
            '  and ',num2str(starttime+wlength),' order by sta,chan,time']);
        fprintf(fid, '%s\n','read waveforms');
        fprintf(fid, '%s\n','write waveforms');
        fclose(fid);
        %====
        disp('***************************************************************')
        disp('****************** query to data base *************************')
        %==== Gparse run
        unix('setenv ORACLE_HOME /cots/oracle/oracle-10.2;');
        unix('setenv D_LIBRARY_PATH $ORACLE_HOME/lib:$ORACLE_HOME/lib32;');
        commandunix = ...
            sprintf('unix(''/ctbto/ims/sm/local/linux/Geotool++/2.3.10/bin/gparse < %sgparse_temp.par;'');',...
            temporary_gparse_dir);
        eval(commandunix);
        if exist('gparse.wfdisc','file')
            commandmove = sprintf('!mv gparse*.* %s.',temporary_gparse_dir);
        else
            filenamesavemat=NaN;
            return
        end
        
        eval(commandmove);
        %====
        disp('***************************************************************')
        disp('***************** convert to Matlab format ********************')
        %====================== Convert to Matlab
        filewfdisc      = sprintf('%sgparse.wfdisc',temporary_gparse_dir);
        
        [filenamesavemat, records{ihc+1}, samprate] = ...
            convertCSStomatlab(filewfdisc,directorydatafromIDC);
        %====================== end ============================
        if isnan(filenamesavemat)
            display('***** .wfdisc does not exist');
        end
    end
    %======================================================
    % to clean etime, stime
    % ================  IMPORTANT =========================
    %======================================================
    
    signals = cell(Msensors,1);
    stime   = zeros(Msensors,1);
    etime   = zeros(Msensors,1);
    
    for im=1:Msensors
        %     filename1_ii = filenames(im).name;
        %     cdload = sprintf('sig = load(''%s%s'');',directorydatafromIDC,filename1_ii);
        %     eval(cdload)
        Lrecords   = length(records{im});
        LL_max=0;
        for ir =1:Lrecords
            if length(records{im}{ir}.data)>LL_max
                irsave = ir;
                LL_max=length(records{im}{ir}.data);
            end
        end
        signals{im} = [records{im}{irsave}.data];
        stime(im) = records{im}{irsave}.stime;
        etime(im) = records{im}{irsave}.etime;
    end
    stimeMAX = max(stime);
    etimeMIN = min(etime);
    signalsproc = cell(Msensors,1);
    
    for im=1:Msensors
        Lim = length(signals{im});
        if stime(im)<=stimeMAX, ds = fix((stimeMAX-stime(im))*Fs_Hz)+1; end
        if etime(im)>=etimeMIN, de = fix((etime(im)-etimeMIN)*Fs_Hz); end
        signalsproc{im} = signals{im}(ds:Lim-de);
    end
    
    N = size([signalsproc{:}],1);
    observations.data = zeros(N,Msensors);
    for im=1:Msensors
        observations.station = stationnumber;
        observations.data(:,im) = [signalsproc{im}];
        observations.xsensors_m = xsensors_m;
    end
    cdesave = sprintf('save %s%s%s%s.mat observations',directorydatafromIDC,yearstart,monthstart,daystart);
    cdesave
    eval(cdesave)
end

%=================================================================

