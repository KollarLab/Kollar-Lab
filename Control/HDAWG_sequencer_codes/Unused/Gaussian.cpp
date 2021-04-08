//Used .cpp extension to get nice syntax highlighting 
//Code to perform T1 like measurements
//_samplesparam_ should be used to specify the number of points in the waveforms
//Wait time can be set either here or by formatting the string

const NumSamples = 800;

wave w1   = gauss(NumSamples,NumSamples/2,NumSamples/8);
//wave w1f  = zeros(NumSamples); //Gaussian starting at 0
//wave mark = marker(NumSamples,1); //Make rectangular window around gaussian pulse to be used as trigger for digitizer or something
//cvar i;

//w1f=w1f+mark;
while(true){
    playWave(w1);
}
