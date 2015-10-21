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
filenamesavemat = convertCSStomatlab(filewfdisc,savedirnamefull);
%====================== end ============================

