//Used .cpp extension to get nice syntax highlighting 
//Code to perform T1 like measurements
//_samplesparam_ should be used to specify the number of points in the waveforms
//Wait time can be set either here or by formatting the string

//Configure pi pulses
const piAmp      = _piAmp_;
const piTime     = _piTime_;
const piWidth    = _piWidth_;

//Timing control 
const qlifetime  = _qlifetime_;
const waitInc    = _waitInc_;
const markerPos  = _markerPos_;

//Control number of repetitions and variations
const averages   = _averages_;
const numConfig  = _numConfig_;

//Convert times to number of samples
const sampleRate      = AWG_RATE_DEFAULT;
const sequencerRate   = sampleRate/8;
var waitInc_cycles    = round(waitInc*sequencerRate);
var qlifetime_cycles  = round(qlifetime*sequencerRate);
var piTime_samples    = round(piTime*sampleRate/16)*16;
var piWidth_samples   = round(piWidth*sampleRate);
var markerPos_samples = round(markerPos*sampleRate);

wave w1        = gauss(piTime_samples,piTime_samples/2,piWidth_samples);
wave w1_final  = zeros(piTime_samples); //Gaussian starting at 0
wave mark_low  = marker(markerPos_samples,0);
wave mark_high = marker(piTime_samples-markerPos_samples,1);
wave w1_mark   = join(mark_low,mark_high);

cvar i;
for(i=0;i<piTime_samples;i++){
  w1_final[i] = w1[i]-w1[0]; //Remove offset from discrete gaussian
}

wave w1_markers = w1_final + w1_mark;

for (i=0 ; i<numConfig ; i++){
  cvar j;
  for (j=0 ; j<averages; j++){
    //Wait for trigger on channel 1
    waitDigTrigger(1)
    playWave(w1_markers);
    waitWave();
    wait(waitInc_cycles*i);
    playWave(w1_final);
    waitWave();
    wait(5*qlifetime_cycles);
  }
}