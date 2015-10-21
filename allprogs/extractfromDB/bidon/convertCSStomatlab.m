function filenamesavemat = convertCSStomatlab(filewfdisc,dirname)
%===============================================================
% convert data from wfdisc into Matlab format .mat
% results are saved in filenamesavemat
%
%===============================================================
fid   = fopen(filewfdisc,'r');
tline = fgetl(fid);
iline = 0;
while ischar(tline)
    iline           = iline+1;
    station{iline}  = strtrim(tline(1:5));
    chan{iline}     = strtrim(tline(6:13));
    stime(iline)    = str2double(tline(14:34));
    wfid(iline)     = str2double(tline(36:43));
    chanid(iline)   = str2double(tline(44:52));
    jdate(iline)    = str2double(tline(53:61));
    etime(iline)    = str2double(tline(62:79));
    nsamp(iline)    = str2double(tline(80:88));
    samprate(iline) = str2double(tline(89:100));
    calib(iline)    = str2double(tline(101:117));
    calper(iline)   = str2double(tline(118:135));
    instype{iline}  = tline(135:136);
    segtype{iline}  = tline(137:144);
    datatype{iline} = tline(144:146);
    clip{iline}     = tline(147:148);
    dir1{iline}     = strtrim(tline(149:209));
    dfile{iline}    = strtrim(tline(210:247));
    foff(iline)     = str2double(tline(248:257));
    commid(iline)   = str2double(tline(258:266));
    lddate{iline}   = tline(267:276);
    tline           = fgetl(fid);
end
fclose(fid);
DY              = num2str(jdate(1));
filenamesavemat = sprintf('%ssta%s_Y%s_D%s.mat',...
    dirname,station{1}(5),DY(1:4),DY(5:end));
% ratiorates    = samprate / min(samprate);
length_record   = fix(etime-stime) .* samprate;
% Read waveforms
fid             = fopen(filewfdisc(1:end-5),'r','b');
signal_temp     = fread(fid,'float32');
fclose(fid);
%========= save data
records = cell(length(wfid),1);
for is = 1 : length(wfid)
    offset              = foff(is)/4;
    signal_is           = signal_temp(offset + 1:offset+length_record(is));
    records{is}.data    = signal_is*calib(is);
    records{is}.Fs_Hz   = samprate(is);
    records{is}.stime   = stime(is);
    records{is}.etime   = etime(is);
    records{is}.station = station{is};
    records{is}.channel = chan{is};
end
cmdsave = sprintf(' save %s records samprate',filenamesavemat);
cmdsave
eval(cmdsave)
%===============================================================
