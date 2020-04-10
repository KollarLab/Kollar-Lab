//Used .cpp extension to get nice syntax highlighting 
//Code to perform T1 like measurements
//_samplesparam_ should be used to specify the number of points in the waveforms
//Wait time can be set either here or by formatting the string

const NumSamples = _NumSamples_;
const waitInc    = _waitInc_;
const markerPos  = _markerPos_;
const averages   = _averages_;
const numMeas    = _numMeas_;
const waitInit   = _waitInit_;

wave w1       = gauss(NumSamples,NumSamples/2,NumSamples/8);
wave w1_final = zeros(NumSamples); //Gaussian starting at 0
wave w1_mark  = marker(NumSamples,1); //Make rectangular window around gaussian pulse to be used as trigger for digitizer or something

cvar i;
for(i=0;i<NumSamples;i++){
  w1_final[i] = w1[i]-w1[0]; //Remove offset from discrete gaussian
}

wave w1_markers = w1_final + w1_mark;

for (i=0 ; i<numMeas ; i++){
  cvar j;
  for (j=0 ; j<averages; j++){
    //Wait for trigger on channel 1
    waitDigTrigger(1)
    playWave(w1_markers);
    waitWave();
    wait(waitInit+waitInc*i);
    playWave(w1_final);
    waitWave();
  }
}