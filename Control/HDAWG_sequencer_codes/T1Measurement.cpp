//Used .cpp extension to get nice syntax highlighting 
//Standard measurement code for two pulses and programmable timing between the two pulses
//Keeps measurement window at the same temporal position as pulse spacing varies

//Configure pi pulses
const piAmp      = _piAmp_;
const piTime     = _piTime_;
const piWidth    = _piWidth_;

//Timing control 
const tau        = _tau_;
const meas_wait  = _meas_wait_;
const meas_time  = _meas_time_;
const max_time   = _max_time_;

//Convert times to number of samples
const sampleRate      = 2.4e+9;
const sequencerRate   = sampleRate/8;
const piTime_samples    = round(piTime*sampleRate/16)*16;
const piWidth_samples   = round(piWidth*sampleRate);

wave w1       = gauss(piTime_samples,piTime_samples/2,piWidth_samples);
wave w1_final = zeros(piTime_samples); //Gaussian starting at 0
wave mark     = marker(piTime_samples,1);
wave w2       = ones(32);
wave w2_2     = zeros(32);

cvar i;
for(i=0;i<piTime_samples;i++){
  w1_final[i] = w1[i]-w1[0]; //Remove offset from discrete gaussian
}

wave w1_markers = w1_final + mark;

const init_wait_cycles = round((max_time-tau-piTime/2)*sequencerRate);
const pulse_sep_cycles = round((tau-piTime)*sequencerRate);
const meas_wait_cycles = round(meas_wait*sequencerRate);
const meas_time_cycles = round(meas_time*sequencerRate); 

while(true){
  //Wait for trigger on channel 1
  waitDigTrigger(1);
  //Always want to measure in the same spot, this timing shifts the first pulse around but keeps everything else fixed
  wait(init_wait_cycles);
  playWave(w1_markers);
  waitWave();
  //Wait the programmed amount of time to measure a T2 for example
  wait(pulse_sep_cycles);
  playWave(w1_markers);
  waitWave();
  //Wait some amount of time until measurement
  wait(meas_wait_cycles);
  //Measurement window (simple TTL signal)
  playWave(2,w2);
  wait(meas_time_cycles);
  playWave(2,w2_2);
}