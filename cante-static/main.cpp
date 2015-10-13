//
//  main.cpp
//  Cante_v1
//
//  Created by Nadine Kroher on 08/08/15.
//  Copyright (c) 2015 Nadine Kroher - Music Technology Group. All rights reserved.

#include <iostream>
#include <fstream>
#include "essentia.h"
#include "algorithmfactory.h"
#include "essentiamath.h"
#include <dirent.h>

#include "jdksmidi/world.h"
#include "jdksmidi/track.h"
#include "jdksmidi/multitrack.h"
#include "jdksmidi/filereadmultitrack.h"
#include "jdksmidi/fileread.h"
#include "jdksmidi/fileshow.h"
#include "jdksmidi/filewritemultitrack.h"


using namespace std;
using namespace essentia;
using namespace essentia::standard;
using namespace jdksmidi;

float pi=3.1416;

// EVALUATE DISTRIBUTION
float evaldist(float _in, float _mu, float _sigma){
    float out=1.0/(_sigma*sqrt(2*pi))*exp(-((_in-_mu)*(_in-_mu))/(2*_sigma*_sigma));
    return out;
}

void findpeaks(const vector<float>& x, vector<float>& _pks, vector<int>& _locs){
    
    if (x.size()<50){
        return;
    }
    
    vector<float>dx,ddx;
    dx.clear();
    ddx.clear();
    for (int i=1; i<x.size(); i++){
        dx.push_back(x[i]-x[i-1]);
    }
    int i=1;
    while (i<dx.size()){
        if (dx[i]==0 && dx[i-1]>0){
            int j=i;
            int count=0;
            while (dx[i]==0){
                i++;
                count++;
            }
            if (dx[i]<0){
                for (int t=j; t<(j+round(float(count)/2)); t++){
                    dx[t]=1;
                }
                for (int t=(j+round(float(count)/2)); t<i; t++){
                    dx[t]=-1;
                }
                i=j+round(float(count)/2);
            }
        }else{
            i++;
        }
    }
    for (int ii=1; ii<dx.size(); ii++){
        ddx.push_back(dx[ii]-dx[ii-1]);
    }
    
    // find locations where first derivative changes sign and second derivative is negative
    for (int ii=1; ii<ddx.size(); ii++){
        if ((dx[ii]*dx[ii-1])<0 && ddx[ii-1]<0){
            _pks.push_back(x[ii]);
            _locs.push_back(ii);
        }
    }
}

int main(int argc, const char * argv[]) {
    
    // GLOBAL PARAMETERS //
    bool polyphonic=true;
    
    // channel selection
    vector<float> audioSample;
    int frameSizeFeat=2048;
    int hopSizeFeat=2048;
    float f1Low=80.0;
    float f1High=400.0;
    float f2Low=500.0;
    float f2High=6000.0;
    vector<float> sbrL, sbrR;
    
    // melody extraction
    int frameSizePitch=1024;
    int hopSizePitch=128;
    
    //float minFreq=120.0;
    float minFreq=70.0; // for QBH
    float maxFreq=720.0;
    float vTh=0.2;
    vector <float> pitch, pitchConfidence;
    
    // bark bands
    int numBands=12;
    vector<float> frame, spectrum, bands;
    vector<vector<float> > voicedBands, unvoicedBands, allBands;
    
    // voicing detection
    vector<float> muV, muU, sigmaV, sigmaU, features;
    vector<float> voicing;
    vector<float> voicingFilt;
    
    // segmentation
    vector<int> startC, endC;
    vector<float> contour, contourFilt;
    vector<float> pks, pksI;
    vector<int> locs, locsI;
    float sigmaG=15;
    float order=1;
    int kernelSize=ceil(sigmaG*(3 + 0.25 * order - 2.5/((order-6)*(order-6)+(order-9)*(order-9))));
    int minIndex=-kernelSize/2;
    vector<float> gaussianKernel;
    for (int i=0; i<kernelSize; i++){
        int n=minIndex+i;
        float val=-n/(sigmaG*sigmaG)*exp((-n*n)/(2*sigmaG*sigmaG));
        gaussianKernel.push_back(val*1/(sigmaG*sqrt(2*3.146)));
    }
    vector<float> h;
    float meanRMS;
    vector<float> zscore;
    
    // pitch labelling
    float fref=440.0;
    vector<float> tuningDev;
    vector<float> tuningRef;
    for (int i=0; i<36; i++){
        tuningRef.push_back(-2400+i*100);
    }
    float Re=0;
    float Im=0;
    vector<float>  windowedFrame, frequencies, magnitudes,locChroma, spectrumC;
    vector<vector<float> > hpcp;
    vector<float> avChroma;
    vector<float> onsetR, durationR, MIDInoteR, velR;
    vector<float> localHist;
    vector<float> sumHist;
    vector<int> keepInd;
    vector<float> onset, duration, MIDInote;
    
    
    // ESSENTIA SETUP //
    essentia::init();
    AlgorithmFactory& factory = standard::AlgorithmFactory::instance();
    // smoothing filter
    Algorithm* mov15   = factory.create("MovingAverage", "size", 15);
    Algorithm* mov30   = factory.create("MovingAverage", "size", 30);
    Algorithm* mov100   = factory.create("MovingAverage", "size", 100);
    // channel selection
    Algorithm* frameCutterL = factory.create("FrameCutter", "frameSize", frameSizeFeat, "hopSize", hopSizeFeat);
    Algorithm* frameCutterR = factory.create("FrameCutter", "frameSize", frameSizeFeat, "hopSize", hopSizeFeat);
    Algorithm* wL     = factory.create("Windowing","type", "blackmanharris62");
    Algorithm* wR     = factory.create("Windowing","type", "blackmanharris62");
    Algorithm* specL  = factory.create("Spectrum");
    Algorithm* specR  = factory.create("Spectrum");
    
    // bark band extraction
    Algorithm* bark = factory.create("BarkBands", "numberBands",numBands);
    Algorithm* frameCutter = factory.create("FrameCutter", "frameSize", frameSizePitch, "hopSize", hopSizePitch);
    Algorithm* spec  = factory.create("Spectrum");
    
    
    // I/O
    vector<string> filenames;
    string audioFilename;
    
    if (argc <= 2) {
        cout << "ERROR: No audio file specified!" << endl;
        cout << "USAGE: " << endl;
        cout << "./cante -f <input file path>" << endl;
        cout << "./cante -r <input folder path>" << endl;
        return 0;
    }else if (argc >=2){
        if (strcmp(argv[1], "-f")==0){
            audioFilename=argv[2];
            filenames.push_back(audioFilename);
        }
        if (strcmp(argv[1],"-r")==0){
            cout << "recursive mode" << endl;
            string foldername = argv[2];
            string last=foldername.string::substr(foldername.size()-1, foldername.size());
            if (strcmp(last.c_str(),"/")!=0){
                foldername=foldername+"/";
            }
            DIR *dir;
            struct dirent *ent;
            if ((dir = opendir (foldername.c_str())) != NULL) {
                // print all the files and directories within directory
                while ((ent = readdir (dir)) != NULL) {
                    string name=ent->d_name;
                    int len=name.size();
                    if (len>4){
                        string ext=name.string::substr(len-4,4);
                        if (strcmp(ext.c_str(),".wav")==0){
                            cout << "file: " << foldername+name << endl;
                            filenames.push_back(foldername+name);
                        }
                    }
                }
                closedir (dir);
            } else {
                cout << "Directory not found" << endl;
                return 0;
            }
        }
        if ((strcmp(argv[1],"-f")!=0) && (strcmp(argv[1],"-r")!=0)){
            cout << "ERROR: Invalid input arguments!" << endl;
            cout << "USAGE: " << endl;
            cout << "./cante -f <input file path>" << endl;
            cout << "./cante -r <input folder path>" << endl;
            return 0;
        }
    }
    
    if (argc >= 4){
        if (strcmp(argv[3], "-m")==0){
            vTh=1.4;
            polyphonic=false;
            cout << "monophonic mode" << endl;
        }
    }
    
    
    // ITERATE OVER FILES;
    for (int f=0; f<filenames.size(); f++){
        
        audioFilename=filenames[f];
        cout << audioFilename << endl;
        
        // check if file exists...
        ifstream testAudio(audioFilename.c_str());
        if (!testAudio.good()){
            cout << "Audio file not found!" << endl;
            continue;
        }
        
        // LOAD AUDIO FILE
        cout << "loading audio file..." << endl;
        vector<StereoSample> audio;
        int numChannels=0;
        Real sampleRate=0;
        Algorithm* audioLoad = factory.create("AudioLoader",
                                              "filename", audioFilename);
        audioLoad->output("audio").set(audio);
        audioLoad->output("sampleRate").set(sampleRate);
        audioLoad->output("numberChannels").set(numChannels);
        audioLoad->compute();
        delete audioLoad;
        
        // SELECT CHANNEL
        cout << "extracting low-level descriptors..." << endl;
        
        int b1Low=round((0.5*frameSizeFeat*f1Low)/(0.5*sampleRate));
        int b1High=round((0.5*frameSizeFeat*f1High)/(0.5*sampleRate));
        int b2Low=round((0.5*frameSizeFeat*f2Low)/(0.5*sampleRate));
        int b2High=round((0.5*frameSizeFeat*f2High)/(0.5*sampleRate));
        
        audioSample.clear();
        
        if (numChannels==2){
            vector<float> left, right, frameL, frameR, windowedFrameL, windowedFrameR, spectrumL, spectrumR;
            // get channel data
            for (int i=1; i<audio.size(); i++){
                left.push_back(audio[i].left());
                right.push_back(audio[i].right());
            }
            // connect algorithms
            frameCutterL->reset();
            frameCutterR->reset();
            frameCutterL->input("signal").set(left);
            frameCutterR->input("signal").set(right);
            frameCutterL->output("frame").set(frameL);
            frameCutterR->output("frame").set(frameR);
            wL->input("frame").set(frameL);
            wR->input("frame").set(frameR);
            wL->output("frame").set(windowedFrameL);
            wR->output("frame").set(windowedFrameR);
            specL->input("frame").set(windowedFrameL);
            specR->input("frame").set(windowedFrameR);
            specL->output("spectrum").set(spectrumL);
            specR->output("spectrum").set(spectrumR);
            
            sbrL.clear();
            sbrR.clear();
            
            
            while(true){
                frameCutterL->compute();
                if (!frameL.size()){
                    break;
                }
                frameCutterR->compute();
                // sepctral band ratio
                wL->compute();
                wR->compute();
                specL->compute();
                specR->compute();
                float highMagL, lowMagL, lowMagR, highMagR;
                lowMagL=0;
                highMagL=0;
                lowMagR=0;
                highMagR=0;
                for (int i=b1Low; i<=b1High; i++){
                    lowMagL+=spectrumL[i];
                    lowMagR+=spectrumR[i];
                }
                for (int i=b2Low; i<=b2High; i++){
                    highMagL+=spectrumL[i];
                    highMagR+=spectrumR[i];
                }
                if (lowMagL!=0 && lowMagR!=0){
                    sbrL.push_back(20*log10(highMagL/lowMagL));
                    sbrR.push_back(20*log10(highMagR/lowMagR));
                }else{
                    sbrL.push_back(0.0);
                    sbrR.push_back(0.0);
                }
            }
            
            // select channel
            if (mean(sbrL, 0, sbrL.size()-1)>mean(sbrR, 0, sbrL.size()-1)){
                for (int i=0; i<left.size(); i++){
                    audioSample.push_back(left[i]);
                }
                cout << "selected left channel" << endl;
            }else{
                for (int i=0; i<left.size(); i++){
                    audioSample.push_back(right[i]);
                }
                cout << "selected right channel" << endl;
            }
        }else{
            for (int i=0; i<audio.size(); i++){
                audioSample.push_back(audio[i].left());
            }
        }
        
        // MELODY EXTRACTION
        cout << "extracting pitch..." << endl;
        Algorithm* melodia     = factory.create("PredominantMelody","sampleRate", sampleRate, "frameSize", frameSizePitch, "hopSize",hopSizePitch, "voiceVibrato", true, "minFrequency", minFreq, "maxFrequency", maxFreq,"voicingTolerance",vTh);
        pitch.clear();
        pitchConfidence.clear();
        melodia->input("signal").set(audioSample);
        melodia->output("pitch").set(pitch);
        melodia->output("pitchConfidence").set(pitchConfidence);
        melodia->compute();
        
        delete melodia;
        
        if(polyphonic){
        
            // CONTOUR FILTERING //
            cout << "filtering contours..." << endl;
            // extract bark bands
            frame.clear();
            spectrum.clear();
            bands.clear();
            frameCutter->reset();
            frameCutter->input("signal").set(audioSample);
            frameCutter->output("frame").set(frame);
            spec->input("frame").set(frame);
            spec->output("spectrum").set(spectrum);
            bark->input("spectrum").set(spectrum);
            bark->output("bands").set(bands);
            voicedBands.clear();
            unvoicedBands.clear();
            allBands.clear();
            int frameCount=0;
            while(true){
                frameCutter->compute();
                if (!frame.size()){
                    break;
                }
                spec->compute();
                bark->compute();
                if (pitch[frameCount]==0){
                    unvoicedBands.push_back(bands);
                }else{
                    voicedBands.push_back(bands);
                }
                if ((frameCount+1)<pitch.size()){
                    frameCount++;
                }
                allBands.push_back(bands);
            }
            
            while (allBands.size()<pitch.size()){
                allBands.push_back(bands);
            }
            while (pitch.size()<allBands.size()){
                pitch.push_back(0);
            }
            if (allBands.size()!=pitch.size()){
                cout << "ERROR: vectors of different length!" << endl;
            }
            
            // fit distributions
            muV.clear();
            muU.clear();
            sigmaV.clear();
            sigmaU.clear();
            features.clear();
            for (int j=0; j<numBands; j++){
                features.clear();
                for (int i=0; i<voicedBands.size(); i++){
                    features.push_back(voicedBands[i][j]);
                }
                muV.push_back(mean(features));
                sigmaV.push_back(stddev(features, mean(features)));
            }
            for (int j=0; j<numBands; j++){
                features.clear();
                for (int i=0; i<unvoicedBands.size(); i++){
                    features.push_back(unvoicedBands[i][j]);
                }
                muU.push_back(mean(features));
                sigmaU.push_back(stddev(features, mean(features)));
            }
            
            // voicing detection
            voicing.clear();
            for (int i=0; i<allBands.size(); i++){
                voicing.push_back(0.0);
            }
            for (int i=0; i<allBands.size(); i++){
                vector<float> xData;
                xData=allBands[i];
                float v=1.0;
                float u=1.0;
                for (int ii=0; ii<numBands; ii++){
                    v*=evaldist(xData[ii], muV[ii], sigmaV[ii]);
                    u*=evaldist(xData[ii], muU[ii], sigmaU[ii]);
                }
                if (u>v){
                    voicing[i]=0;
                }else{
                    voicing[i]=1;
                }
            }
            voicingFilt.clear();
            mov100->input("signal").set(voicing);
            mov100->output("signal").set(voicingFilt);
            mov100->compute();
            for (int i=0; i<voicingFilt.size(); i++){
                if (voicingFilt[i]<0.3){
                    voicingFilt[i]=0;
                }else{
                    voicingFilt[i]=1.0;
                }
            }
            
            
            // find contour start and end indices
            startC.clear();
            endC.clear();
            if (pitch[0]>0){
                startC.push_back(0);
            }
            for (int i=0; i<pitch.size()-1; i++){
                if (pitch[i+1]>0 && pitch[i]==0){
                    startC.push_back(i+1);
                }
                if (pitch[i+1]==0 && pitch[i]>0){
                    endC.push_back(i);
                }
            }
            if (endC.size()<startC.size()){
                endC.push_back(pitch.size()-1);
            }
            
            // eliminate contours below voicing threshold
            for (int i=0; i<startC.size(); i++){
                if (sum(voicingFilt, startC[i], endC[i])==0.0){
                    for (int ii=startC[i]; ii<=endC[i]; ii++){
                        pitch[ii]=0;
                    }
                }
            }
        
        }
        // find contour start and end indices
        startC.clear();
        endC.clear();
        if (pitch[0]>0){
            startC.push_back(0);
        }
        for (int i=0; i<pitch.size()-1; i++){
            if (pitch[i+1]>0 && pitch[i]==0){
                startC.push_back(i+1);
            }
            if (pitch[i+1]==0 && pitch[i]>0){
                endC.push_back(i);
            }
        }
        if (endC.size()<startC.size()){
            endC.push_back(pitch.size()-1);
        }
        
        
        // SEGMENT CONTOURS
        cout << "segmenting contours..." << endl;
        
        // re-segment
        // find contour start and end indices
        startC.clear();
        endC.clear();
        if (pitch[0]>0){
            startC.push_back(0);
        }
        for (int i=0; i<pitch.size()-1; i++){
            if (pitch[i+1]>0 && pitch[i]==0){
                startC.push_back(i+1);
            }
            if (pitch[i+1]==0 && pitch[i]>0){
                endC.push_back(i);
            }
        }
        if (endC.size()<startC.size()){
            endC.push_back(pitch.size()-1);
        }
        
        if (endC[endC.size()-1]>pitch.size()){
            cout << "ERROR: index exceeds pitch vector!" << endl;
        }
        if (!startC.size()){
            cout << "no vocal contours found!" << endl;
            return 0;
        }
        
        // segment based on local maxima
        for (int j=0; j<startC.size(); j++){
            if ((endC[j]-startC[j])>50){
                contour.clear();
                contourFilt.clear();
                for (int i=startC[j]; i<=endC[j]; i++){
                    contour.push_back(1200*3.322038403*log10 (pitch[i]/55.0));
                }
                mov30->reset();
                mov30->input("signal").set(contour);
                mov30->output("signal").set(contourFilt);
                mov30->compute();
                pksI.clear();
                locsI.clear();
                findpeaks(contourFilt, pksI, locsI);
                pks.clear();
                locs.clear();
                int minInd=25;
                for (int i=0; i<pksI.size(); i++){
                    float locMean = 0.0;
                    int minIndM=max(0,locsI[i]-10);
                    int maxIndM=min(int(contourFilt.size()-1), locsI[i]+10);
                    int count_j=0;
                    for (int j=minIndM; j<=maxIndM; j++){
                        locMean+=contourFilt[j];
                        count_j++;
                    }
                    locMean/=count_j;
                    if (pksI[i]>locMean){
                        pks.push_back(pksI[i]);
                        locs.push_back(locsI[i]);
                    }
                }
                
                for (int i=1; i<pks.size(); i++){
                    if (abs(pks[i]-pks[i-1])>80){
                        int ind=locs[i-1]+round(0.5*(float(float(locs[i])-float(locs[i-1]))));
                        if (ind>minInd && ind<contourFilt.size()-25){
                            pitch[startC[j]+ind-15]=0;
                            minInd=ind+25;
                        }
                    }
                }
            }
        }
        
        // find contour start and end indices
        startC.clear();
        endC.clear();
        if (pitch[0]>0){
            startC.push_back(0);
        }
        for (int i=0; i<pitch.size()-1; i++){
            if (pitch[i+1]>0 && pitch[i]==0){
                startC.push_back(i+1);
            }
            if (pitch[i+1]==0 && pitch[i]>0){
                endC.push_back(i);
            }
        }
        if (endC.size()<startC.size()){
            endC.push_back(pitch.size()-1);
        }
        
        
        // segment based on gaussian filter output
        for (int j=0; j<startC.size(); j++){
            
            if ((endC[j]-startC[j])>50){
                contour.clear();
                contourFilt.clear();
                for (int i=startC[j]; i<=endC[j]; i++){
                    contour.push_back(1200*log2(pitch[i]/440.0));
                }
                mov15->reset();
                mov15->input("signal").set(contour);
                mov15->output("signal").set(contourFilt);
                mov15->compute();
                float m=mean(contourFilt, 0, contourFilt.size()-1);
                for (int i=0; i<contourFilt.size(); i++){
                    contourFilt[i]-=m;
                }
                int numSamples=contourFilt.size();
                float y[numSamples+kernelSize-1];
                for ( int ii = kernelSize; ii < numSamples; ii++ )
                {
                    y[ii] = 0;                       // set to zero before sum
                    for ( int jj = 0; jj < kernelSize; jj++ )
                    {
                        y[ii] += contourFilt[ii - jj] * gaussianKernel[jj];
                    }
                }
                h.clear();
                for (int ii=kernelSize; ii<numSamples; ii++){
                    h.push_back(abs(y[ii]));
                }
                pks.clear();
                locs.clear();
                findpeaks(h, pks, locs);
                int xLim=0;
                for (int ii=0; ii<locs.size(); ii++){
                    if (pks[ii]>5.0 && locs[ii]>xLim && locs[ii]>25 && locs[ii]<numSamples-25){
                        pitch[startC[j]+locs[ii]+kernelSize/2]=0;
                        xLim=locs[ii]+25+kernelSize/2;
                    }
                }
            }
            
        }

        
        // find contour start and end indices
        startC.clear();
        endC.clear();
        if (pitch[0]>0){
            startC.push_back(0);
        }
        for (int i=0; i<pitch.size()-1; i++){
            if (pitch[i+1]>0 && pitch[i]==0){
                startC.push_back(i+1);
            }
            if (pitch[i+1]==0 && pitch[i]>0){
                endC.push_back(i);
            }
        }
        if (endC.size()<startC.size()){
            endC.push_back(pitch.size()-1);
        }
        
        
        
        // segment based on RMS
        for (int j=0; j<startC.size(); j++){
            if ((endC[j]-startC[j])>100){
                int rs=startC[j]*128;
                int re=endC[j]*128;
                rs=max(0, rs);
                re=min(int(audioSample.size()),re);
                meanRMS=0.0;
                int countRMS=0;
                for (int i=rs; i<=re; i++){
                    meanRMS+=audioSample[i]*audioSample[i];
                    countRMS++;
                }
                meanRMS/=countRMS;
                meanRMS=sqrt(meanRMS);
                contour.clear();
                for (int i=startC[j]; i<=endC[j]; i++){
                    float r=0.0;
                    float r1=(i-5)*128;
                    float r2=(i+5)*128;
                    int cRMS=0;
                    for (int ii=r1; ii<=r2; ii++){
                        r+=audioSample[ii]*audioSample[ii];
                        cRMS++;
                    }
                    r/=cRMS;
                    r=sqrt(r);
                    contour.push_back(-20*log10(r/meanRMS));
                }
                pks.clear();
                locs.clear();
                findpeaks(contour, pks, locs);
                int numSamples=contour.size();
                int xLim=0;
                for (int ii=0; ii<locs.size(); ii++){
                    if (pks[ii]>10.0){
                        pitch[startC[j]+locs[ii]]=0;
                        xLim=locs[ii]+25;
                    }
                }
            }
        }
        
        // find contour start and end indices
        startC.clear();
        endC.clear();
        if (pitch[0]>0){
            startC.push_back(0);
        }
        for (int i=0; i<pitch.size()-1; i++){
            if (pitch[i+1]>0 && pitch[i]==0){
                startC.push_back(i+1);
            }
            if (pitch[i+1]==0 && pitch[i]>0){
                endC.push_back(i);
            }
        }
        if (endC.size()<startC.size()){
            endC.push_back(pitch.size()-1);
        }
        
        
        // segment at contour dips
        for (int j=0; j<startC.size(); j++){
            
            if ((endC[j]-startC[j])>100){
                contour.clear();
                contourFilt.clear();
                for (int i=startC[j]; i<=endC[j]; i++){
                    contour.push_back(1200*log2(pitch[i]/440.0));
                }
                mov15->reset();
                mov15->input("signal").set(contour);
                mov15->output("signal").set(contourFilt);
                mov15->compute();
                zscore.clear();
                float m=mean(contourFilt, 0, contourFilt.size()-1);
                float s=stddev(contourFilt, m);
                for (int i=0; i<contourFilt.size(); i++){
                    zscore.push_back(-((contourFilt[i]-m)/s));
                }
                pks.clear();
                locs.clear();
                findpeaks(zscore, pks, locs);
                int minIndex=25;
                for (int ii=0; ii<locs.size(); ii++){
                    if (pks[ii]>2 && locs[ii]<zscore.size()-25 && locs[ii]>minIndex){
                        pitch[startC[j]+locs[ii]]=0;
                        minIndex=locs[ii]+25;
                    }
                }
            }
        }
        
        // find contour start and end indices
        startC.clear();
        endC.clear();
        if (pitch[0]>0){
            startC.push_back(0);
        }
        for (int i=0; i<pitch.size()-1; i++){
            if (pitch[i+1]>0 && pitch[i]==0){
                startC.push_back(i+1);
            }
            if (pitch[i+1]==0 && pitch[i]>0){
                endC.push_back(i);
            }
        }
        if (endC.size()<startC.size()){
            endC.push_back(pitch.size()-1);
        }
        
        // TUNING ESTIMATION
        cout << "estimate tuning..." << endl;
        float ftune;
        tuningDev.clear();
        for (int i=0; i<startC.size(); i++){
            contour.clear();
            contourFilt.clear();
            for (int ii=startC[i]; ii<=endC[i]; ii++){
                contour.push_back(1200*log2(pitch[ii]/440.0));
            }
            mov30->reset();
            mov30->input("signal").set(contour);
            mov30->output("signal").set(contourFilt);
            mov30->compute();
            for (int ii=0; ii<contourFilt.size(); ii++){
                tuningDev.push_back(contourFilt[ii]-round(contourFilt[ii]/100)*100);
            }
        }
        
        Re=0.0;
        Im=0.0;
        for (int i=0; i<tuningDev.size(); i++){
            tuningDev[i]*=(2*3.146/100);
        }
        for (int i=0; i<tuningDev.size(); i++){
            Re+=cos(tuningDev[i]);
            Im+=sin(tuningDev[i]);
            
        }

        float deltaTune=atan2(Im, Re);
        deltaTune*=100.0/(2.0*3.146);
        ftune=pow(10, deltaTune/(1200*3.32203))*fref;
        cout << "delta tune: " << deltaTune << endl;
        
        // GLOBAL PITCH PROBABILITY
        cout << "estimate global pitch probability..." << endl;
        frameCutter->reset();
        Algorithm* w     = factory.create("Windowing","type", "blackmanharris62");
        Algorithm* specC = factory.create("Spectrum");
        Algorithm* specPeaks = factory.create("SpectralPeaks");
        Algorithm* chroma = factory.create("HPCP", "referenceFrequency", ftune, "sampleRate", sampleRate,"weightType", "cosine", "harmonics", 50);
        
        frame.clear();
        windowedFrame.clear();
        spectrumC.clear();
        frequencies.clear();
        magnitudes.clear();
        locChroma.clear();
        hpcp.clear();
        avChroma.clear();
        
        frameCutter->reset();
        frameCutter->input("signal").set(audioSample);
        frameCutter->output("frame").set(frame);
        w->input("frame").set(frame);
        w->output("frame").set(windowedFrame);
        spec->input("frame").set(windowedFrame);
        spec->output("spectrum").set(spectrumC);
        specPeaks->input("spectrum").set(spectrumC);
        specPeaks->output("frequencies").set(frequencies);
        specPeaks->output("magnitudes").set(magnitudes);
        chroma->input("frequencies").set(frequencies);
        chroma->input("magnitudes").set(magnitudes);
        chroma->output("hpcp").set(locChroma);
        
        while (true) {
            frameCutter->compute();
            if (!frame.size()){
                break;
            }
            w->compute();
            spec->compute();
            specPeaks->compute();
            chroma->compute();
            hpcp.push_back(locChroma);
        }
        
        delete w;
        delete specC;
        delete specPeaks;
        delete chroma;
        
        for (int i=0; i<12; i++){
            avChroma.push_back(0.0);
        }
        for (int i=0; i<hpcp.size(); i++){
            for (int ii=0; ii<12; ii++){
                avChroma[ii]+=hpcp[i][ii];
            }
        }
        float s=sum(avChroma);
        for (int i=0; i<12; i++){
            avChroma[i]/=s;
        }
        
        // PITCH LABELLING
        onsetR.clear();
        durationR.clear();
        MIDInoteR.clear();
        cout << "assign pitch labels..." << endl;
        for (int i=0; i<startC.size(); i++){
            onsetR.push_back(float(startC[i])*float(hopSizePitch)/float(sampleRate));
            durationR.push_back(float(endC[i]-startC[i])*float(hopSizePitch)/float(sampleRate));
            // local histogram
            localHist.clear();
            sumHist.clear();
            for (int ii=0; ii<36; ii++){
                localHist.push_back(0.0);
            }
            for (int ii=startC[i]; ii<=endC[i];ii++){
                float c=1200*log2(pitch[ii]/ftune);
                int pBin=round(c/100)+24;
                pBin=max(0,pBin);
                pBin=min(35,pBin);
                localHist[pBin]+=+1;
            }
            sumHist.clear();
            for (int ii=0; ii<localHist.size(); ii++){
                sumHist.push_back(0.0);
            }
            for (int ii=1; ii<localHist.size()-1; ii++){
                sumHist[ii-1]+=0.1080*localHist[ii];
                sumHist[ii+1]+=0.1080*localHist[ii];
            }
            for (int ii=0; ii<localHist.size()-1; ii++){
                localHist[ii]+=sumHist[ii];
            }
            float totalHist=sum(localHist);
            for (int ii=0; ii<localHist.size()-1; ii++){
                localHist[ii]/=totalHist;
            }
            for (int ii=0; ii<localHist.size(); ii++){
                int k=ii % 12;
                localHist[ii]*=avChroma[k];
            }
            int maxBin=argmax(localHist);
            int cQuant=tuningRef[maxBin];
            float fQuant=(pow(10,cQuant/(1200*3.3220))*fref);
            MIDInoteR.push_back(round(12*log2(fQuant/fref)+69));
        }
        
        // POST-PROCESS NOTES
        cout << "note post-processing..." << endl;
        keepInd.clear();
        onset.clear();
        duration.clear();
        MIDInote.clear();
        
        float meanNote=mean(MIDInoteR, 0, MIDInoteR.size()-1);
        for (int i=0; i<durationR.size(); i++){
            if (MIDInoteR[i]>(meanNote+8)){
                MIDInoteR[i]-=12;
            }
            if (durationR[i]>0.05 && MIDInoteR[i]>(meanNote-12)){
                keepInd.push_back(i);
            }
        }
        for (int i=0; i<keepInd.size(); i++){
            MIDInote.push_back(MIDInoteR[keepInd[i]]);
            duration.push_back(durationR[keepInd[i]]);
            onset.push_back(onsetR[keepInd[i]]);
        }
        
        // eliminate temprary isolated notes
        for (int i=2; i<MIDInote.size()-1; i++){
            if (duration[i]<1.0){
                if ((onset[i]-onset[i-1])>1.0 && (onset[i+1]-onset[i])>1.0){
                    onset.erase(onset.begin()+i);
                    MIDInote.erase(MIDInote.begin()+i);
                    duration.erase(duration.begin()+i);
                }
            }
        }
        
        // WRITE CSV FILE //
        cout << "write text file..." << endl;
        int len=audioFilename.size();
        string csvFilename = audioFilename.string::substr(0,len-4)+".notes.csv";
        cout << csvFilename << endl;
        ofstream outFile;
        outFile.open(csvFilename);
        for (int ii=0; ii<onset.size(); ii++){
            outFile << onset[ii] << ", " << duration[ii] << ", " << MIDInote[ii]<< endl;
        }
        outFile.close();
        
        // WRITE MIDI FILE //
        cout << "write MIDI file" << endl;
        int return_code = -1;
        MIDITimedBigMessage m; // the object for individual midi events
        unsigned char chan, // internal midi channel number 0...15 (named 1...16)
        note, velocity, ctrl, val;
        MIDIClockTime t; // time in midi ticks
        MIDIClockTime dt; // time interval (1 second)
        int clks_per_beat = 100; // number of ticks in crotchet (1...32767)
        int num_tracks = 2; // tracks 0 and 1
        MIDIMultiTrack tracks( num_tracks );  // the object which will hold all the tracks
        tracks.SetClksPerBeat( clks_per_beat );
        int trk; // track number, 0 or 1
        
        t = 0;
        m.SetTime( t );
        // track 0 is used for tempo and time signature info, and some other stuff
        trk = 0;
        int tempo = 500000; // set tempo to 1 000 000 usec = 1 sec in crotchet
        // with value of clks_per_beat (100) result 10 msec in 1 midi tick
        // If no tempo is define, 120 beats per minute is assumed.
        m.SetTempo( tempo );
        tracks.GetTrack( trk )->PutEvent( m );
        
        // create midi events
        trk = 1;
        // create synch note
        t = 0;
        dt = 0.1*200; // (0.1 seconds...)
        m.SetTime( t );
        m.SetNoteOn( chan = 0, note = 33, velocity = 100 ); // A1
        tracks.GetTrack( trk )->PutEvent( m );
        m.SetTime( t += dt );
        m.SetNoteOff( chan, note, velocity );
        tracks.GetTrack( trk )->PutEvent( m );
        // add melody notes
        for (int i=0; i<MIDInote.size(); i++){
            t=onset[i]*200;
            dt=duration[i]*200;
            m.SetTime(t);
            m.SetNoteOn( chan = 0, note = MIDInote[i], velocity = 100 ); // A1
            tracks.GetTrack( trk )->PutEvent( m );
            m.SetTime( t += dt );
            m.SetNoteOff( chan, note, velocity );
            tracks.GetTrack( trk )->PutEvent( m );
        }
        
        
        string outfile_name = audioFilename.string::substr(0,len-4)+".mid";
        MIDIFileWriteStreamFileName out_stream( outfile_name.c_str() );
        
        // then output the stream like my example does, except setting num_tracks to match your data
        
        if( out_stream.IsValid() )
        {
            // the object which takes the midi tracks and writes the midifile to the output stream
            MIDIFileWriteMultiTrack writer( &tracks, &out_stream );
            
            // write the output file
            if ( writer.Write( num_tracks ) )
            {
                cout << "\nOK writing file " << outfile_name << endl;
                return_code = 0;
            }
            else
            {
                cerr << "\nError writing file " << outfile_name << endl;
            }
        }
        else
        {
            cerr << "\nError opening file " << outfile_name << endl;
        }
        
    } // end file loop
    
    // ESSENTIA CLEAN UP
    // smoothing filter
    delete mov15;
    delete mov30;
    delete mov100;
    // channel selection
    delete frameCutterL;
    delete frameCutterR;
    delete wL;
    delete wR;
    delete specL;
    delete specR;
    // bark bands
    delete frameCutter;
    delete spec;
    delete bark;
    
    
    essentia::shutdown();
    
    return 0;
}
