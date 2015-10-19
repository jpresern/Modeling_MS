%%  Fig2A stimuli generation. Increasing conditioning stimuli amplitude
%   while keeping the test stimuli the same (8.5 um) beginning at the 
%   fixed point


stim = [0, 0.8, 1.5, 2.2, 2.9, 3.6, 4.3, 5, 5.7, 6.4, 7.1, 7.8, 8.5];
% maxSlopeTime = 10.4;
% % maxSlopeTime = 5.3125;
% slopeTime = linspace (0, maxSlopeTime, 14);
slopeTime = [0,1,1.875,2.75,3.625,4.5,5.375,6.25,7.125,8,8.875,9.75,10.625];

t = NaN (13, 9);
amp = NaN (13, 9); 
for k = 1:13
    timeGen = [10, slopeTime(k), 10, maxSlopeTime+(maxSlopeTime - slopeTime(k)),0,maxSlopeTime, 10, 20, 10];
%     timeGen = [10, slopeTime(k), 10, 10-slopeTime(k),0,maxSlopeTime, 10, 20, 10];
    t (k,:) = timeGen;
    
    ampGen = [0, stim(k), 0, -stim(k) , 0, 8.5, 0, -8.5, 0];
    amp (k,:) = ampGen;

end;

Fig2A_Stimulus = struct('amp',amp,'t',t);

%  Drawing stimuli pair amp,t.

figure;
    for val2 = 1:13
        x0 = 0;
        y0 = 0;
        for val1 = 1:9
            x1 = x0 + Fig2A_Stimulus.t(val2,val1);
            y1 = y0 + Fig2A_Stimulus.amp(val2,val1);
            line([x0 x1],[y0 y1]);
            x0 = x1;
            y0 = y1;
        end;
%         pause;
    end;

%%  Fig3C. stimuli generation, each pair in the structure has different conditioning stimuli
%   14 different conditioning stimuli used, starting with 0.7 um
%   first test stimuli in all fourteen pairs is 0.0 um
%   each of the 14 batches contain 6 layers, each with its own delay
%   between conditioning stimuli and test stimuli.
%   structure with 14 pairs (amp,t) is created.

stim = [0, 0.8, 1.5, 2.2, 2.9, 3.6, 4.3, 5, 5.7, 6.4, 7.1, 7.8, 8.5];
stim2 = [0, 0.7, 1.4, 2.1, 2.8, 3.5, 4.2, 4.9, 5.6, 6.3, 7.0, 7.7,8.4];
% % slopeTime = linspace (0, 12.250, 13);
% slopeTime = linspace (0, 10.4, 14);
slopeTime = [0,1,1.875,2.75,3.625,4.5,5.375,6.25,7.125,8,8.875,9.75,10.625];
delay = [4, 20, 40, 100, 400, 1000];
Fig3C_Stimulus = struct([]);
for i = 1:9                        %   counts the different amplitudes used desenzitization Poke
%    k = 0;      
    t = NaN (13, 7, 6);
    amp = NaN (13, 7, 6); 
    for j = 1:6                   %   counts the different delays used
            l = i;             %    l is used to start writing at the correct line in t and amp
            for k = 13:-1:i    %   counts the different amplitudes used for morePoke(TM) - has to be done in reverse
                %   condition correctly matches the amplitudes used in
                %   Hao2010 experiments
                if i == 1
                    timeGen = [5, slopeTime(i), delay(j), slopeTime(13-k+1), 40, slopeTime(13-k+1),20];
                    ampGen = [0, stim(i), 0, stim(13-k+1) , 0, -stim(13-k+1), 0];
                else
                    timeGen = [5, slopeTime(i), delay(j), slopeTime(13-k+1), 40, slopeTime(13-k+1),20];
                    ampGen = [0, stim(i), 0, stim2(13-k+1) , 0, -stim2(13-k+1), 0];                    
                end;
                t (l,:,j) = timeGen;
                amp (l,:,j) = ampGen;
                l = l + 1;
           end;     
    
    end;
    Fig3C_Stimulus = [Fig3C_Stimulus, struct('amp',amp, 't',t)]; 
end;

%  Drawing stimuli pair amp,t.

figure;
for figg = 1: length(Fig3C_Stimulus)
subplot(3,5,figg);
for val3 = 1:6
    for val2 = 1:13
        x0 = 0;
        y0 = 0;
        for val1 = 1:6
            x1 = x0 + Fig3C_Stimulus(1,figg).t(val2,val1,val3);
            y1 = y0 + Fig3C_Stimulus(1,figg).amp(val2,val1,val3);
            line([x0 x1],[y0 y1]);
            x0 = x1;
            y0 = y1;
        end;
    end;
end;
end;


%%  Calculates stimuli for Hao 2010 Fig 8 - different slope velocities (=different slope times)

% stim = [0, 0.8, 1.5, 2.2, 2.9, 3.6, 4.3, 5, 5.7, 6.4, 7.1, 7.8, 8.5];
stim = 8.5;

slopeTime = [10.625   82.6250   64.6250 50.3125   46.6250   28.6250 21.250 17.6375 14.1313 10.625];



t = NaN (10, 5);
amp = NaN (10, 5); 
for k = [1:10]
    timeGen = [10, slopeTime(k), 40, slopeTime(3) , 20];
    t (k,:) = timeGen;
    ampGen = [0, stim, 0, -stim , 0];
    amp (k,:) = ampGen;

end;

Fig8A_Stimulus = struct('amp',amp,'t',t);

%  Drawing stimuli pair amp,t.

figure;
    for val2 = 1:10
        x0 = 0;
        y0 = 0;
        for val1 = 1:5
            x1 = x0 + Fig8A_Stimulus.t(val2,val1);
            y1 = y0 + Fig8A_Stimulus.amp(val2,val1);
            line([x0 x1],[y0 y1]);
            x0 = x1;
            y0 = y1;
        end;
%         pause;
    end;

%%  stimuli for Hao2010 Fig 6A-C
%   Increasing duration of conditioning stimuli in the poke-repoke protocol


dur = [5, 10, 20, 50, 100, 200, 500, 600, 1000];

stim = [0, 0.8, 1.5, 2.2, 2.9, 3.6, 4.3, 5, 5.7, 6.4, 7.1, 7.8, 8.5];
% stim2 = [0, 0.7, 1.4, 2.1, 2.8, 3.5, 4.2, 4.9, 5.6, 6.3, 7.0, 7.7,8.4];

%single velocity and single amplitude used
slopeTime = [0,1,1.875,2.75,3.625,4.5,5.375,6.25,7.125,8,8.875,9.75,10.625];

% Fig6A_Stimulus = struct([]);
t = NaN (13, 9, 9);
amp = NaN (13, 9, 9); 
for i = 1:9
    for k = 1:13
            timeGen = [10, slopeTime(k), dur(i), slopeTime(k), 0,slopeTime(k), 10,slopeTime(k), 10];
            t (k,:,i) = timeGen;
            ampGen = [0, stim(k), 0, -stim(k), 0, stim(k), 0, -stim(k), 0];
            amp (k,:,i) = ampGen;
    end;
end; 

Fig6A_Stimulus = struct('amp',amp,'t',t);

%  Drawing stimuli pair amp,t.

figure;
for val3 = 1:9
    subplot (3,3,val3)
    for val2 = 1:13
        x0 = 0;
        y0 = 0;
        for val1 = 1:9
            x1 = x0 + Fig6A_Stimulus.t(val2,val1,val3);
            y1 = y0 + Fig6A_Stimulus.amp(val2,val1,val3);
            line([x0 x1],[y0 y1]);
            x0 = x1;
            y0 = y1;
        end;
    end;
    pause;
end;

%%
%%  Fig3C. stimuli generation, each pair in the structure has different conditioning stimuli
%   14 different conditioning stimuli used, starting with 0.7 um
%   first test stimuli in all fourteen pairs is 0.0 um
%   each of the 14 batches contain 6 layers, each with its own delay
%   between conditioning stimuli and test stimuli.
%   structure with 14 pairs (amp,t) is created.

stim = [0, 0.8, 1.5, 2.2, 2.9, 3.6, 4.3, 5, 5.7, 6.4, 7.1, 7.8, 8.5, 9.2, 9.9, 10.6, 11.3, 12.0, 12.7, 13.4,14.1];
stim2 = [0, 0.7, 1.4, 2.1, 2.8, 3.5, 4.2, 4.9, 5.6, 6.3, 7.0, 7.7,8.4, 9.1, 9.8, 10.5,11.2,11.9,12.6,13.3,14.0];
% slopeTime = linspace (0, 10.625, 13);
slopeTime = linspace (0, 16.8, 22);
slopeTime = [0,1,1.875,2.75,3.625,4.5,5.375,6.25,7.125,8,8.875,9.75,10.625,11.50,12.375,13.25,14.125,15,15.875,16.75,17.625];
delay = [4, 20, 40, 100, 400, 1000];
Fig3C_Stimulus = struct([]);
for i = 1:9                        %   conditioning: counts the different amplitudes used desenzitization Poke
%    k = 0;      
    t = NaN (21, 7, 6);
    amp = NaN (21, 7, 6); 
    for j = 1:6                   %   counts the different delays used
            l = i;             %    l is used to start writing at the correct line in t and amp
            for k = 21:-1:i    %   counts the different amplitudes used for morePoke(TM) - has to be done in reverse
                %   condition correctly matches the amplitudes used in
                %   Hao2010 experiments
                if i == 1
                    timeGen = [5, slopeTime(i), delay(j), slopeTime(21-k+1), 40, slopeTime(21-k+1),20];
                    ampGen = [0, stim(i), 0, stim(21-k+1) , 0, -stim(21-k+1), 0];
                else
                    timeGen = [5, slopeTime(i), delay(j), slopeTime(21-k+1), 40, slopeTime(21-k+1),20];
                    ampGen = [0, stim(i), 0, stim2(21-k+1) , 0, -stim2(21-k+1), 0];                    
                end;
                t (l,:,j) = timeGen;
                amp (l,:,j) = ampGen;
                l = l + 1;
           end;     
    
    end;
    Fig3C_Stimulus = [Fig3C_Stimulus, struct('amp',amp, 't',t)]; 
end;

%  Drawing stimuli pair amp,t.

figure;
for figg = 1: length(Fig3C_Stimulus)
subplot(3,5,figg);
for val3 = 1:6
    for val2 = 1:21
        x0 = 0;
        y0 = 0;
        for val1 = 1:6
            x1 = x0 + Fig3C_Stimulus(1,figg).t(val2,val1,val3);
            y1 = y0 + Fig3C_Stimulus(1,figg).amp(val2,val1,val3);
            line([x0 x1],[y0 y1]);
            x0 = x1;
            y0 = y1;
        end;
    end;
end;
end;

%%
%%  Fig2A stimuli generation. Increasing conditioning stimuli amplitude
%   while keeping the test stimuli the same (8.5 um) beginning at the 
%   fixed point


stim = [0, 0.8, 1.5, 2.2, 2.9, 3.6, 4.3, 5, 5.7, 6.4, 7.1, 7.8, 8.5];
% maxSlopeTime = 10.4;
% % maxSlopeTime = 5.3125;
% slopeTime = linspace (0, maxSlopeTime, 14);
slopeTime = [0,1,1.875,2.75,3.625,4.5,5.375,6.25,7.125,8,8.875,9.75,10.625];
maxSlopeTime = max(slopeTime);

t = NaN (13, 9);
amp = NaN (13, 9); 
for k = 1:13
%     timeGen = [10, slopeTime(k), 10, maxSlopeTime ,0,maxSlopeTime, 10, 20, 10];
    timeGen = [10, slopeTime(k), 10, maxSlopeTime + slopeTime(14-k),0,maxSlopeTime, 10, 20, 10];
    t (k,:) = timeGen;
    
    ampGen = [0, stim(k), 0, -stim(k) , 0, 8.5, 0, -8.5, 0];
    amp (k,:) = ampGen;

end;

Fig2A_Stimulus = struct('amp',amp,'t',t);

%  Drawing stimuli pair amp,t.

figure;
    for val2 = 1:13
        x0 = 0;
        y0 = 0;
        for val1 = 1:9
            x1 = x0 + Fig2A_Stimulus.t(val2,val1);
            y1 = y0 + Fig2A_Stimulus.amp(val2,val1);
            line([x0 x1],[y0 y1]);
            x0 = x1;
            y0 = y1;
        end;
%         pause;
    end;

