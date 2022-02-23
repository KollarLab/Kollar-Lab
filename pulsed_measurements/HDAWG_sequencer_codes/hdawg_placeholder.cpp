const samples = _samples_;
wave ch1 = zeros(samples);
wave ch2 = zeros(samples);

//Note if we don't set the marker to anything it won't actually create a WF
//for it. Setting to 3 for half the time so that we have a placeholder
//This need to talk to BOTH marker channels!! Otherwise it just ignores
//the upload
wave m = marker(samples/2, 3);

ch1 = ch1+m;

while(true){
  waitDigTrigger(1);
  playWave(ch1, ch2);
  waitWave();
}