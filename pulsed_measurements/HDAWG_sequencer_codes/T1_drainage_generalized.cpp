// Used .cpp extension to get nice syntax highlighting 
// Generalized T1_drainage sequencer script: if tau is negative, the code will play an overlap sequence
// (i.e. leave the qubit tone on during the measurement tone for a time tau). Otherwise it will perform
// the standard T1 drainage (keep the qubit tone on most of the time, turn it off a time tau before 
// measurement). In both cases, the qubit tone stays off for a time 1.2*measurement window after the measurement
// window to allow for background collection

//Configure 
const measTime = _meas_window_;
const tau = _tau_;
//Timing control 
const max_time   = _max_time_;

//Convert times to number of samples
const sampleRate      = 2.4e+9;
const sequencerRate   = sampleRate/8;

//Measurement tone params
const meas_samples      = measTime*sampleRate;
const meas_ramp_samples = 160;

//The -0.5 is critical to have a symmetric gaussian (so that the signal goes to 0 at both ends)
wave meas_ramp = gauss(meas_ramp_samples, meas_ramp_samples/2-0.5, meas_ramp_samples/8);
wave meas_rise = cut(meas_ramp,0,meas_ramp_samples/2-1);
wave meas_fall = cut(meas_ramp,meas_ramp_samples/2,meas_ramp_samples-1);
wave rise_clean = zeros(meas_ramp_samples/2);
wave fall_clean = zeros(meas_ramp_samples/2);
wave blank     = zeros(meas_ramp_samples/2);
wave hold_high = ones(meas_samples);

const ramp_off = meas_ramp[0];
cvar i; 
for(i=0; i<meas_ramp_samples/2; i++){
  rise_clean[i] = meas_rise[i]-ramp_off;
  fall_clean[i] = meas_fall[i]-ramp_off;
}
wave measurement_tone = join(rise_clean, hold_high, fall_clean);

//qubit cycle
const tau_samples = round(abs(tau)*sampleRate);
wave init_high = ones(tau_samples);
wave final_low = zeros(meas_samples-tau_samples+meas_ramp_samples/2);
wave qubit_cycle = join(init_high, fall_clean, final_low);

const init_wait_cycles_standard = round((max_time-abs(tau))*sequencerRate);
const init_wait_cycles_overlap  = round(max_time*sequencerRate);
const pulse_sep_cycles = round(abs(tau)*sequencerRate);
// Keep the qubit pulse low for 1.2 measurement windows after the measurement pulse to allow for background collection
const background_cycles = round(measTime*sequencerRate*1.2);

const tau_pos = round((sign(tau)+1)/2);
while(true){
  //Wait for trigger on channel 1
  waitDigTrigger(1);
  if (tau_pos){
    wait(init_wait_cycles_standard);
    //Drop the qubit tone to 0
    playWave(blank, fall_clean);
    waitWave();
    wait(pulse_sep_cycles);
    //Measurement tone
    playWave(measurement_tone);
    waitWave();
  } else {
    wait(init_wait_cycles_overlap);
    playWave(measurement_tone, qubit_cycle);
    waitWave();
  }
  wait(background_cycles);
  playWave(blank, rise_clean);
  waitWave();
}