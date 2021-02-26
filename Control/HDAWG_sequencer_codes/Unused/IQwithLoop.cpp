const NumPoints = _npoints_;
const Time      = _time_;

const a = _a_;
const b = _b_;
const c = _c_;
const d = _d_;

const sampleRate      = 2.4e+9;
const Time_samples    = round(Time*sampleRate);

wave out  = ones(Time_samples);
wave sin  = sine(NumPoints,1.0,0.,1.);
wave cos  = cosine(NumPoints,1.0,0.,1.);

for (cvar i = 0; i < NumPoints; i++) {
    waitDigTrigger(1);
    playWave((a*cos[i] + b*sin[i])*out,(c*cos[i] + d*sin[i])*out);
    waitWave();
}