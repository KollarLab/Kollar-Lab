//Used .cpp extension to get nice syntax highlighting 

//Configure tones
const Time     = 100e-6;

//Convert times to number of samples
const sampleRate      = 2.4e+9;
const sequencerRate   = sampleRate/8;
const Time_samples    = Time*sequencerRate;
const ramp_samples    = 160;

//The -0.5 is critical to have a symmetric gaussian (so that the signal goes to 0 at both ends)
wave ramp = gauss(ramp_samples, ramp_samples/2-0.5, ramp_samples/8);
wave rise = cut(ramp,0,ramp_samples/2-1);
wave fall = cut(ramp,ramp_samples/2,ramp_samples-1);

while(true){
  //Wait for trigger on channel 1
  waitDigTrigger(1);
  //Always want to measure in the same spot, this timing shifts the first pulse around but keeps everything else fixed
  wait(5);
  playWave(rise);
  waitWave();
  //Wait the programmed amount of time for next pulse
  wait(Time_samples);
  playWave(fall);
  waitWave();
}