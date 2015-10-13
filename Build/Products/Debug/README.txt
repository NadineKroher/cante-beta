CANTE 
=====

CANTE is a software for automatic flamenco singing transcription. This is a prototype binary for demonstration purposes to accompany the IEEE TASP submission "Automatic Transcription of Flamenco Singing from Polyphonic Music Recordings". 
Algorithm and code developed by Nadine Kroher, Universitat Pompeu Fabra, 2015. 


Usage
------

./cante -f <filename> for transcribing a single .wav file
./cante -r <foldername> for transcribing all wav files in a folder
./cante -f <filename> -m for transcribing monophonic recordings

The output .notes.csv folder is placed in the same location as the audio file. The format is specified as follows:

<Onset time[s]>, <Duration [s]>, <MIDI pitch>;

Each row corresponds to a note.

Furthermore, a midi file is created for each track. To facilitat synchronization, an A0 note is inserted at time 0.0sec. 

Version 
--------

Cante_v0_beta
Mac OSX
Tested under OSX 10.9.5