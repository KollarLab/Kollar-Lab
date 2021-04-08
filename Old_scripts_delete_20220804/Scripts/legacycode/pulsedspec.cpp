//Used .cpp extension to get nice syntax highlighting 

//Configure 
const qbitTime = 100e-9;
const measTime = 100e-6;

//Convert times to number of samples
const sampleRate      = 2.4e+9;
const sequencerRate   = sampleRate/8;


//Measurement tone params
const meas_samples      = measTime*sequencerRate;
const meas_ramp_samples = 160;

//The -0.5 is critical to have a symmetric gaussian (so that the signal goes to 0 at both ends)
wave meas_ramp = gauss(meas_ramp_samples, meas_ramp_samples/2-0.5, meas_ramp_samples/8);
wave meas_rise = cut(meas_ramp,0,meas_ramp_samples/2-1);
wave meas_fall = cut(meas_ramp,meas_ramp_samples/2,meas_ramp_samples-1);


//Qubit pulse params
const frac = 0.1
const qbit_samples = qbitTime*sampleRate;
const ramp_samples = max(round(qbit_samples*frac/16)*16,32);

wave ramp  = gauss(ramp_samples, ramp_samples/2-0.5, ramp_samples/8);
wave rise  = cut(ramp,0,ramp_samples/2-1);
wave fall  = cut(ramp,ramp_samples/2,ramp_samples-1);
wave sig   = ones(Time_samples);
wave blank = zeros(ramp_samples+qubit_samples);

//Round the edges of the square pulse by adding a gaussian ramp up and down around it
wave qubit = join(rise, sig, fall);

//Subtract offset from data (gaussian doesn't go to 0 at the ends initially)
const offset = ramp[0];
cvar i;
for (i=0; i<ramp_samples+qubit_samples;i++){
  qubit[i] = qubit[i]-offset;
}

while(true){
  //Wait for trigger on channel 1
  waitDigTrigger(1);
  wait(5);
  //Qubit tone
  playWave(blank, sig);
  waitWave();
  wait(50);
  //Measurement tone
  playWave(rise);
  waitWave();
  wait(meas_samples);
  playWave(fall);
  waitWave();
  
  playWave(blank, blank);
  waitWave();
}