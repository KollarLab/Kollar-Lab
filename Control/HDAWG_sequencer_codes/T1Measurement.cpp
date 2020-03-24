//Used .cpp extension to get nice syntax highlighting 
//Code to perform T1 like measurements
//_samplesparam_ should be used to specify the number of points in the waveforms
//Wait time can be set either here or by formatting the string

const NumSamples = 800;
const waitInc    = 1000;
const clearTime  = 1000;
const markerPos  = 300;

wave w1   = gauss(NumSamples,NumSamples/2,NumSamples/8);
wave w1f  = zeros(NumSamples); //Gaussian starting at 0
wave mark = marker(NumSamples,1); //Make rectangular window around gaussian pulse to be used as trigger for digitizer or something
cvar i;
for(i=0;i<NumSamples;i++){
  w1f[i]=w1[i]-w1[0]; //Remove offset from discrete gaussian
}
w1f=w1f+mark;
for (i=0;i<5;i++){
  playWave(w1f);
  waitWave();
  wait(waitInc*i);
  playWave(w1f);
  waitWave();
  wait(clearTime);
}