//Used .cpp extension to get nice syntax highlighting 
//Standard measurement code for two pulses and programmable timing between the two pulses
//Keeps measurement window at the same temporal position as pulse spacing varies
//Implements I/Q angle to provide full control over the pulse

//Configure tones
const Amp  = _Amp_;
const Time = _Time_;
const Iamp = _Iamp_;
const Qamp = _Qamp_;

//Ramp config
const frac = _frac_;

//Convert times to number of samples
const sampleRate      = 2.4e+9;
const sequencerRate   = sampleRate/8;
const Time_samples    = round(Time*sampleRate/16)*16;
const ramp_samples    = max(round(Time_samples*frac/16)*16,32);

//The -0.5 is critical to have a symmetric gaussian (so that the signal goes to 0 at both ends)
wave ramp = gauss(ramp_samples, ramp_samples/2-0.5, ramp_samples/8);
wave rise = cut(ramp,0,ramp_samples/2-1);
wave fall = cut(ramp,ramp_samples/2,ramp_samples-1);
wave sig  = ones(Time_samples);
//Round the edges of the square pulse by adding a gaussian ramp up and down around it
wave tone = join(rise, sig, fall);

//Subtract offset from data (gaussian doesn't go to 0 at the ends initially)
const offset = ramp[0];
cvar i;
for (i=0; i<ramp_samples+Time_samples;i++){
  tone[i] = tone[i]-offset;
}

//Perform scaling to user specified height 
tone = Amp*tone/(1-offset);

while(true){
  //Wait for trigger on channel 1
  waitDigTrigger(1);
  playWave(Iamp*tone, Qamp*tone);
  waitWave();
}