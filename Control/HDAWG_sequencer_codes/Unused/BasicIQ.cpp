const Time = _time_;
const I    = _Iamp_;
const Q    = _Qamp_;

const sampleRate      = 2.4e+9;
const Time_samples    = round(Time*sampleRate);

wave out  = ones(Time_samples);

while (true) {
    playWave(out*I, out*Q);
}