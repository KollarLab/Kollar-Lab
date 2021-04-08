//Used .cpp extension to get nice syntax highlighting 
//Standard measurement code for two pulses and programmable timing between the two pulses
//Keeps measurement window at the same temporal position as pulse spacing varies
//Implements I/Q angle to provide full control over the pulse

//Configure tones
const Time = _Time_;

//Convert times to number of samples
const sampleRate      = 2.4e+9;
const sequencerRate   = sampleRate/8;
const Time_samples    = round(Time*sampleRate/16)*16;

wave sig  = zeros(Time_samples);

while(true){
  //Wait for trigger on channel 1
  waitDigTrigger(1);
  playWave(sig, sig);
  waitWave();
}