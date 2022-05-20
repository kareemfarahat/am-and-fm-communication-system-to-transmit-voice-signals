sampling_frequensy = 48000;
recorded = audiorecorder(sampling_frequensy,8,1);
disp('Start speaking.')
recordblocking(recorded, 30);
disp('End of Recording.');
message = getaudiodata(recorded);
audiowrite("test message.mp3",message,sampling_frequensy);