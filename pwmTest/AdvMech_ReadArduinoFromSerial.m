clear all
close all
clc

% closes all open serial connections
% make sure that the arduino serial monitor is not open!
delete(instrfindall);	

% change the serial parameters based on your setup
arduino=serial('COM14','BaudRate',115200); 

disp('connecting to arduino...') 
fopen(arduino);

% change the third parameter to the number of reads you want to take
% to control frequency, you can add a timestamp check before printing to the serial monitor
x=linspace(1,1000,2e3);

disp('starting loop...')

for i=1:length(x)

	% read two float values from serial monitor
	% change this to fit your set up
	% in the arduino code, just have 
	y(i)=fscanf(arduino,'%f'); %#ok<SAGROW>
	z(i)=fscanf(arduino,'%f'); %#ok<SAGROW>
end

% arduino data is available after this point

disp('making plot...')

% plot raw and filtered data, zoomed in
plot(x,y);
axis([0,inf,0,1.3])
hold all
plot(x,z);

% double check that the serial connection is closed
fclose(arduino);
delete(instrfindall);
