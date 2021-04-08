//Used .cpp extension to get nice syntax highlighting 

//Configure 
const measTime = _meas_window_;
const overlap = _overlap_;
//Timing control 
const max_time   = _max_time_;

//Convert times to number of samples
const sampleRate      = 2.4e+9;
const sequencerRate   = sampleRate/8;

//Measurement tone params
const meas_samples      = measTime*sequencerRate;
const overlap_samples   = overlap*sequencerRate;
const meas_ramp_samples = 160;

//The -0.5 is critical to have a symmetric gaussian (so that the signal goes to 0 at both ends)
wave meas_ramp = gauss(meas_ramp_samples, meas_ramp_samples/2-0.5, meas_ramp_samples/8);
wave meas_rise = cut(meas_ramp,0,meas_ramp_samples/2-1);
wave meas_fall = cut(meas_ramp,meas_ramp_samples/2,meas_ramp_samples-1);
wave rise_clean = zeros(meas_ramp_samples/2);
wave fall_clean = zeros(meas_ramp_samples/2);
wave blank = zeros(meas_ramp_samples/2);
wave top = ones(meas_ramp_samples/2);

const ramp_off = meas_ramp[0];
cvar i; 
for(i=0; i<meas_ramp_samples/2; i++){
  rise_clean[i] = meas_rise[i]-ramp_off;
  fall_clean[i] = meas_fall[i]-ramp_off;
}
const ramp_top = rise_clean[meas_ramp_samples/2];

const init_wait_cycles = round(max_time*sequencerRate);

while(true){
  //Wait for trigger on channel 1
  waitDigTrigger(1);
  wait(init_wait_cycles);
  //Measurement tone
  playWave(rise_clean, blank);
  waitWave();
  wait(overlap_samples);
  playWave(2, fall_clean);
  wait(meas_samples - overlap_samples);
  playWave(fall_clean, rise_clean);
  waitWave();
}