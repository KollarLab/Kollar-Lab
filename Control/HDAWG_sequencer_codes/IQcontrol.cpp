const Time      = _time_;
const NumPoints = _npoints_;
const Phi       = _angimp_;
const AlphaI    = _ampimpI_;
const AlphaQ    = _ampimpQ_;

const sampleRate      = 2.4e+9;
const sequencerRate   = sampleRate/8;
const Time_samples    = round(Time*sampleRate/16)*16;

wave w1   = ones(Time_samples);
wave sin  = AlphaQ*sine(NumPoints,1.0,Phi,2.);
wave cos  = AlphaI*cosine(NumPoints,1.0,0.,2.);

while(true){
  for (cvar i = 0; i < NumPoints; i++) {
    waitDigTrigger(1);
    playWave(w1*cos[i],w1*sin[i]);
  }
}