def SGS_coupling(logen, rfgen, ext_ref='HDAWG', ext_ref_freq=10, coupling='Ref', SGS_ref_freq=1000):
    if ext_ref == 'HDAWG':
        rfgen.set_External_Reference(freq = ext_ref_freq)
    else:
        rfgen.set_Internal_Reference()
    if coupling == 'LO':
        rfgen.set_RefLO_output(output='LO')
    else:
        rfgen.set_RefLO_output(output='Ref', freq = SGS_ref_freq)
    
    if coupling == 'Ref':
        logen.set_External_Reference(freq = SGS_ref_freq)
        logen.set_Internal_LO()
    if coupling == 'LO':
        logen.set_External_LO()

def HDAWG_clock(hdawg, freqs=[], channels=[], amps=[]):
    for freq, channel, amp in zip(freqs, channels, amps):
        # Configure artificial clock using HDAWG
        hdawg.OSCs[1].configure_sine(channel%2,freq)
        hdawg.Channels[channel].analog_outs = [(channel-1)%2,channel%2]
        hdawg.Channels[channel].configureChannel(amp=amp)