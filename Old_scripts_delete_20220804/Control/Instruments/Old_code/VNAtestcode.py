# -*- coding: utf-8 -*-
"""
Created on Fri Sep 11 14:25:33 2020

@author: Kollarlab
"""

#add = 'TCPIP0::192.168.1.15::inst0::INSTR'
#inst = rm.open_resource(add)
#inst.query('*IDN?')
#inst.query('SYST:ERR?')
#inst.write('*RST')
#
#inst.write('INIT:CONT OFF')
#inst.write('CALC:PAR:DEL:ALL')
#
#inst.write("CALC:PAR:DEF:SGR 1,2")
#inst.write("DISP:WIND:TRAC2:FEED 'Ch1_SG_S11'")
#inst.write("DISP:WIND:TRAC3:FEED 'Ch1_SG_S12'")
#inst.write("DISP:WIND:TRAC4:FEED 'Ch1_SG_S21'")
#inst.write("DISP:WIND:TRAC5:FEED 'Ch1_SG_S22'")


def getAvgTrace(numAvg):
    inst.write('SENS:AVER ON')
    inst.write('SENS:AVER:CLE')
    inst.write('SENS:AVER:COUN {}'.format(numAvg))
    inst.write('SENS:AVER:MODE AUTO')
    inst.write('SWE:COUN {}'.format(numAvg)) 
      
    inst.write('INIT:IMM; *OPC')
    opc = 0
    while not (opc&1):
        print('waiting for command to complete')
        opc = eval(inst.query('*ESR?').rstrip())
        sleep(1)
        
    sdat = inst.query_ascii_values("CALC:DATA:SGR? FDAT")
    pylab.figure()
    for i in range(4):
        l = int(len(sdat)/4)
        pylab.plot(sdat[i*l:(i+1)*l])