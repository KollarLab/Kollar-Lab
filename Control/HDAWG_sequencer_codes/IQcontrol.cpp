const Angle     = _angle_;
const Time      = _time_;

const sampleRate      = 2.4e+9;
const sequencerRate   = sampleRate/8;
const Time_samples    = round(Time*sampleRate/16)*16;

wave w1 = ones(Time_samples);

while(true){
    playWave(w1*cos(Angle), w1*sin(Angle));
}