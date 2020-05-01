//Used .cpp extension to get nice syntax highlighting 
//Standard measurement code for two pulses and programmable timing between the two pulses
//Keeps measurement window at the same temporal position as pulse spacing varies

//Configure tones
const Amp      = _piAmp_;
const Time     = _piTime_;

//Ramp config
const frac = 10;

//Timing control 
const tau        = _tau_;
const max_time   = _max_time_;

//Convert times to number of samples
const sampleRate      = 2.4e+9;
const sequencerRate   = sampleRate/8;
const Time_samples    = round(piTime*sampleRate/16)*16;
const ramp_samples    = round(Time_samples/(16*frac))*16;

wave ramp = gauss(2*ramp_samples, ramp_samples, ramp_samples/2);
wave rise = cut(ramp,0,ramp_samples-1);
wave fall = cut(ramp,ramp_samples,2*ramp_samples-1);
wave sig  = ones(Time);
wave tone = Amp*join(rise, sig, fall);

const init_wait_cycles = round((max_time-tau-(1+2./frac)Time/2)*sequencerRate);
const pulse_sep_cycles = round((tau-(1+2./frac)Time)*sequencerRate);

while(true){
  //Wait for trigger on channel 1
  waitDigTrigger(1);
  //Always want to measure in the same spot, this timing shifts the first pulse around but keeps everything else fixed
  wait(init_wait_cycles);
  playWave(tone);
  waitWave();
  //Wait the programmed amount of time to measure a T2 for example
  wait(pulse_sep_cycles);
  playWave(tone);
  waitWave();
}
