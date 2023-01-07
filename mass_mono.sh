#!/bin/bash

for FILE in *.wav;
do
ffmpeg -i $FILE -map_channel 0.0.0 temp.wav;
mv temp.wav $FILE
done;
