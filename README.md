# CANTE - A software tool for automatic transcription of flamenco singing.

VERSON cante_v0_beta
This is a beta release of the CANTE software. The algorithm estimates a note-level transcription of the singing voice from monophonic or accompanied flamenco recordings. 

CANTE was developed and implemented by Nadine Kroher in the context of the COFLA (cofla-project.com) research project. 

This repository contains the c++ source code as well as a static build for OSX (located in /Build/Products/Debug/). This is a beta version, tested under OSX 10.9.5.

 *  Copyright (C) 2015  Nadine Kroher
 *  www.cofla-project.com
 *  ginsonicsound@gmail.com
 

# Usage

./cante -f [filename] for transcribing a single .wav file

./cante -r [foldername] for transcribing all wav files in a folder

./cante -f [filename] -m for transcribing monophonic recordings

The output .notes.csv folder is placed in the same location as the audio file. The format is specified as follows:

[Onset time[s]], [Duration [s], [MIDI pitch];

Each row corresponds to a note.

Furthermore, a midi file is created for each track. To facilitate synchronization, an A0 note is inserted at time 0.0sec. 


# Dependencies
CANTE uses the ESSENTIA library (predominant melody extraction) and jdkmidi (source files included). The XCode is set up as a static build, i.e. the resulting binary should only depend on system libraries. 

# License
 *  CANTE is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  CANTE is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 
 *  You should have received a copy of the GNU General Public License
 *  along with Foobar.  If not, see <http://www.gnu.org/licenses/>. **

All rights reserved.

# Does CANTE work for other genres?
CANTE has been developed for and tested on the specific case of (accompanied) flamenco singing and is tailored to the genre's specific musical characteristics. Feel free to test it on other styles of accompanied singing and share your experiences / results. 
