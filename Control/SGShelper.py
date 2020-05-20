
def SGS_coupling(ext_ref='HDAWG', ext_ref_freq=10, coupling='Ref', SGS_ref_freq=1000):
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
