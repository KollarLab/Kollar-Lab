//Used .cpp extension to get nice syntax highlighting 

//Configure 
const qbit_sigma = _qsigma_;
const num_sigma  = _num_sigma_;
const min_q_time = qbit_sigma*num_sigma;
const measTime   = _meas_window_;

//Timing control 
const hold_time  = _hold_time_;
const max_time   = _max_time_;
const meas_delay = _meas_delay_;

//Convert times to number of samples
const sampleRate      = 2.4e+9;
const sequencerRate   = sampleRate/8;

//Measurement tone params
const meas_samples      = measTime*sampleRate;
const hold_samples      = hold_time*sampleRate;
const meas_ramp_samples = 160;

//The -0.5 is critical to have a symmetric gaussian (so that the signal goes to 0 at both ends)
wave meas_ramp = gauss(meas_ramp_samples, meas_ramp_samples/2-0.5, meas_ramp_samples/8);
wave meas_rise = cut(meas_ramp,0,meas_ramp_samples/2-1);
wave meas_fall = cut(meas_ramp,meas_ramp_samples/2,meas_ramp_samples-1);
wave rise_clean = zeros(meas_ramp_samples/2);
wave fall_clean = zeros(meas_ramp_samples/2);

wave hold_high = ones(meas_samples);
wave blank = zeros(meas_ramp_samples);

const ramp_off = meas_ramp[0];
cvar i; 
for(i=0; i<meas_ramp_samples/2; i++){
  rise_clean[i] = meas_rise[i]-ramp_off;
  fall_clean[i] = meas_fall[i]-ramp_off;
}
wave measurement_tone = join(rise_clean, hold_high, fall_clean);

//Qubit pulse params
const qbit_samples = min_q_time*sampleRate;
const qbit_position = qbit_samples/2-0.5;
const qbit_amp = 1.0;

wave qbit  = gauss(qbit_samples, qbit_amp, qbit_position, qbit_sigma);
wave q_rise = cut(qbit,0,qbit_samples/2-1);
wave q_fall = cut(qbit,qbit_samples/2,qbit_samples-1);
wave q_rise_clean = zeros(qbit_samples/2);
wave q_fall_clean = zeros(qbit_samples/2);

//wave q_hold_high = ones(round(hold_samples/16)*16);
wave q_hold_high = ones(hold_samples);

const q_off = qbit[0];

for(i=0; i<qbit_samples/2; i++){
  q_rise_clean[i] = q_rise[i]-q_off;
  q_fall_clean[i] = q_fall[i]-q_off;
}
wave q_tone = join(q_rise_clean, (1-q_off)*q_hold_high, q_fall_clean);

const init_wait_cycles = round((max_time-meas_delay-min_q_time-hold_time)*sequencerRate);
const pulse_sep_cycles = round(meas_delay*sequencerRate);

while(true){
  //Wait for trigger on channel 1
  waitDigTrigger(1);
  wait(init_wait_cycles);
  //Qubit tone
  playWave(2, q_tone);
  waitWave();
  wait(pulse_sep_cycles);
  //Measurement tone
  playWave(measurement_tone);
  waitWave();
  
  playWave(blank, blank);
  waitWave();
}