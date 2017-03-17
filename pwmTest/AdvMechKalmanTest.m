clear all
close all
clc
delete(instrfindall);	% closes all serial connections

arduino=serial('COM14','BaudRate',115200);

disp('connecting to arduino...') 
fopen(arduino);

x=linspace(1,1000,2e3);

disp('starting loop...')

for i=1:length(x)

	% read two float values from serial monitor
	y(i)=fscanf(arduino,'%f'); %#ok<SAGROW>
	z(i)=fscanf(arduino,'%f'); %#ok<SAGROW>
end


disp('making plot...')
plot(x,y);
axis([0,inf,0,1.3])
hold all
plot(x,z);

fclose(arduino);
delete(instrfindall);
