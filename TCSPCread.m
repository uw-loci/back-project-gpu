function [data, resolution_s, time] = TCSPCread(filename)
%TCSPCread  Reads HidraHarp ASCII measurement data files - 2014 SPADLab (Mauro B.)
fid = fopen(filename);
header_1 = fgetl(fid);
header_2 = fgetl(fid);

channel_number = fscanf(fid, '%d\r\n');

header_3 = fgetl(fid);
header_4 = fgetl(fid);
header_5 = fgetl(fid);
header_6 = fgetl(fid);
header_7 = fgetl(fid);

resolution_ns = fscanf(fid, '%f\t%*f\t%*f\t%*f\t\r\n');
resolution_s = resolution_ns / 1e9;

header_8 = fgetl(fid);

data = fscanf(fid, '%d\t%*d\t%*d\t%*d\t\r\n', channel_number);

time = (0:resolution_s:resolution_s*channel_number-resolution_s);

fclose(fid);