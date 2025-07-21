
# -*- coding: utf-8 -*-
"""
Modified by jhyang, 7/17/2025
---
Added upon by Bridgette McAllister, SignalCore 7/8/2022
---
Created on Fri Mar 19 14:11:31 2021
Adapted from Chao Zhou
---
A simple driver for SignalCore SC5511A, transferred from the one written by Erick Brindock
"""

import ctypes
import logging
from typing import Any, Dict, Optional
#from pprint import pprint

#definitions for the register dictionary
INITIALIZE = 0x01
SET_SYS_ACTIVE = 0x02
SYNTH_MODE = 0x03
RF_MODE = 0x04
LIST_MODE_CONFIG = 0x05
LIST_START_FREQ = 0x06
LIST_STOP_FREQ = 0x07
LIST_STEP_FREQ = 0x08
LIST_DWELL_TIME = 0x09
LIST_CYCLE_COUNT = 0x0A
RESERVED0 = 0x0B
LIST_BUFFER_POINTS = 0x0C
LIST_BUFFER_WRITE = 0x0D
LIST_BUF_MEM_XFER = 0x0E
LIST_SOFT_TRIGGER = 0x0F

RF_FREQUENCY = 0x10
RF_LEVEL = 0x11
RF_ENABLE = 0x12
RF_PHASE = 0x13
AUTO_LEVEL_DISABLE = 0x14
RF_ALC_MODE = 0x15
RF_STANDBY = 0x16
REFERENCE_MODE = 0x17
REFERENCE_DAC_VALUE = 0x18
ALC_DAC_VALUE = 0x19
RESERVED2 = 0x1A
STORE_DEFAULT_STATE = 0x1B
RESERVED3 = 0x1C
RESERVED4 = 0x1D
RF2_STANDBY = 0x1E
RF2_FREQUENCY = 0x1F

GET_RF_PARAMETERS = 0x20
GET_TEMPERATURE = 0x21
GET_DEVICE_STATUS = 0x22
GET_DEVICE_INFO = 0x23
GET_LIST_BUFFER = 0x24
GET_ALC_DAC_VALUE = 0x25
GET_SERIAL_OUT_BUFFER = 0x26
GET_EEPROM_VALUE = 0x27

SYNTH_SELF_CAL = 0x47

#note that all registers 0x27 to 0x50 are reserved for factory use. Writing to them accidentally may
#cause the device to functionally fail

class Device_rf_params_t(ctypes.Structure):
    _fields_ = [("rf1_freq", ctypes.c_ulonglong),   #current ch#1 rf frequency
                ("start_freq", ctypes.c_ulonglong), #list start frequency
                ("stop_freq", ctypes.c_ulonglong),  #list stop frequency (>start_freq)
                ("step_freq", ctypes.c_ulonglong),  #list step frequency
                ("sweep_dwell_time", ctypes.c_uint),    #dwell time at each frequency
                ("sweep_cycles", ctypes.c_uint),    #number of cycles to sweep/list
                ("buffer_points", ctypes.c_uint),   #current number of list buffer points
                ("rf_level", ctypes.c_float),   #current ch#1 power level
                ("rf2_freq", ctypes.c_short)    #current ch#2 rf frequency
                ]

class RF2_t(ctypes.Structure):
    _fields_ = [("rf2_freq", ctypes.c_short)]
    
class Level_t(ctypes.Structure):
    _fields_ = [("level", ctypes.c_float)]

class Phase_t(ctypes.Structure):
    _fields_ = [("phase", ctypes.c_float)]

class Device_temperature_t(ctypes.Structure):
    _fields_ = [("device_temp", ctypes.c_float)]

class Dac_value_t(ctypes.Structure):
    _fields_ = [("dac_value", ctypes.c_ushort)]

class List_buffer_t(ctypes.Structure):
    _fields_ = [("transfer_mode", ctypes.c_ubyte)]

class List_buffer_set_t(ctypes.Structure):
    _fields_ = [("device_freq", ctypes.c_ulonglong)]

class Reg_read_t(ctypes.Structure):
    _fields_ = [("reg_read", ctypes.c_ulonglong)]


class Operate_status_t(ctypes.Structure):
    _fields_ = [("rf1_lock_mode", ctypes.c_ubyte),  #synthesizer lock mode for chn#1: 0 = use harmonic curcuit, 1 = fracN circuit
                ("rf1_loop_gain", ctypes.c_ubyte),  #changing the loop gain of the sum pll. 0 = normal, 1 = low. Low gain helps suppress spurs and far out phase noise, but increases the close in phase
                ("device_access", ctypes.c_ubyte),  #if a session has been open for the device
                ("rf2_standby", ctypes.c_ubyte),    #indicates chn#2 standby and output disable
                ("rf1_standby", ctypes.c_ubyte),    #indicates chn#1 standby
                ("auto_pwr_disable", ctypes.c_ubyte),   #indicates power adjustment is performed when frequency is changed
                ("alc_mode", ctypes.c_ubyte),   #indicates alc behavior: 0 is closed, 1 is open
                ("rf1_out_enable", ctypes.c_ubyte), #indicates chn#1 RF output
                ("ext_ref_lock_enable", ctypes.c_ubyte),    #indicates that 100 MHz VCXO is set to lock to an external source
                ("ext_ref_detect", ctypes.c_ubyte), #indicates external source detected
                ("ref_out_select", ctypes.c_ubyte), #indicates the reference output select: 0 = 10 MHz, 1 = 100 MHz
                ("list_mode_running", ctypes.c_ubyte),  #indicates list/sweep is triggered and currently running
                ("rf1_mode", ctypes.c_ubyte),   #indicates chn#1 rf mode set: 0=fixed tone state, 1 = list/sweep mode state
                ("harmonic_ss", ctypes.c_ubyte),    #harmonic spur suppression state
                ("over_temp", ctypes.c_ubyte)   #indicates if the temperature of the devices has exceeded ~75 deg C internally
                ]


class Pll_status_t(ctypes.Structure):
    _fields_ = [("sum_pll_ld", ctypes.c_ubyte), #lock status of main pll loop
                ("crs_pll_ld", ctypes.c_ubyte), #lock status of coarse offset pll loop (used only for harmonic mode)
                ("fine_pll_ld", ctypes.c_ubyte),    #lock status of the dds tuned fine pll loop
                ("crs_ref_pll_ld", ctypes.c_ubyte), #lock status of the coarse reference pll loop
                ("crs_aux_pll_ld", ctypes.c_ubyte), #lock status of the auxiliary coarse pll loop (used only for IntN or FracN mode)
                ("ref_100_pll_ld", ctypes.c_ubyte), #lock status of the 100 MHz VCXO pll loop
                ("ref_10_pll_ld", ctypes.c_ubyte),  #lock status of the master 10 MHz TCXO pll loop
                ("rf2_pll_ld", ctypes.c_ubyte)] #lock status of the chn#2 pll loop


class List_mode_t(ctypes.Structure):
    _fields_ = [("sss_mode", ctypes.c_ubyte),   #0 uses list for buffer, 1 calculates using stop-start-step 
                ("sweep_dir", ctypes.c_ubyte),  #0 start/beginning to stop/end, 1 stop/end to start/beginning
                ("tri_waveform", ctypes.c_ubyte),   #0 sawtooth, 1 triangular
                ("hw_trigger", ctypes.c_ubyte),     #0 soft trigger expected, 1 hard trigger expected
                ("step_on_hw_trig", ctypes.c_ubyte),    #0 trigger to sweep through list, 1 stepping on every trigger (on hard trigger only)
                ("return_to_start", ctypes.c_ubyte),    #if 1, frequency returns to start frequency after end of cycle(s)
                ("trig_out_enable", ctypes.c_ubyte),    #1 enables a trigger pulse at the trigger on pin
                ("trig_out_on_cycle", ctypes.c_ubyte)]  #0 trigger out on every frequency change, 1 trigger on cycle complete


class Device_status_t(ctypes.Structure):
    _fields_ = [("list_mode", List_mode_t), #list mode parameters
                ("operate_status_t", Operate_status_t), #operating parameters
                ("pll_status_t", Pll_status_t)] #pll status

class Manufacture_date_t(ctypes.Structure):
    _fields_ = [("year", ctypes.c_ubyte),
                ("month", ctypes.c_ubyte),
                ("day", ctypes.c_ubyte),
                ("hour", ctypes.c_ubyte)]

class Device_info_t(ctypes.Structure):
    _fields_ = [("serial_number", ctypes.c_uint32),
                ("hardware_revision", ctypes.c_float),
                ("firmware_revision", ctypes.c_float),
                ("manufacture_date", Manufacture_date_t)
                ]

class Clock_config_t(ctypes.Structure):
    _fields_ = [("ext_ref_lock_enable", ctypes.c_ubyte),    #select to lock to either external reference source
                ("ref_out_select", ctypes.c_ubyte),     #select either 10 or 100 MHz output
                ("pxi_clock_enable", ctypes.c_ubyte),   #
                ("ext_direct_clocking", ctypes.c_ubyte),    #enable direct 100 MHz clocking of the synthesizer, bypassing its internal 100 MHz clock
                ("ext_ref_freq", ctypes.c_ubyte)]   #selects the input frequency

# End of Structures------------------------------------------------------------     
 
 
def getdict(struct):
    """
    This is copied from online: 
    https://stackoverflow.com/questions/3789372/python-can-we-convert-a-ctypes-structure-to-a-dictionary
    """
    result = {}
    for field, _ in struct.fields_:
         value = getattr(struct, field)
         # if the type is not a primitive and it evaluates to False ...
         if (type(value) not in [int, float, bool]) and not bool(value):
             # it's a null pointer
             value = None
         elif hasattr(value, "_length_") and hasattr(value, "_type_"):
             # Probably an array
             value = list(value)
         elif hasattr(value, "_fields_"):
             # Probably another struct
             value = getdict(value)
         result[field] = value
    return result


class SignalCore_SC5511A(object):
    
    def __init__(self, name: str, serial_number: str, dll_path: str = 'C:\\Program Files\\SignalCore\\SC5511A\\api\\c\\x64\\sc5511a.dll', debug = False, **kwargs: Any):

        # variables appear with an underscore prefix, while their corresponding properties do not.        

        super().__init__()

        self.init = True
        self.verbose = False

        logging.info(f'Initializing instrument SignalCore generator {name}, #{serial_number}')
        self.dll = ctypes.WinDLL(dll_path)
        if debug:
            print(self.dll)
        self.dll.sc5511a_open_device.restype = ctypes.c_uint64
        self.handle = ctypes.c_void_p(self.dll.sc5511a_open_device(ctypes.c_char_p(bytes(serial_number, 'utf-8'))))
        self.serial_number = ctypes.c_char_p(bytes(serial_number, 'utf-8'))
        
        # Cannot use self.open = 0, must use the function
        # SC.set_open(0) or SC._close()
        # property setup for self.open will self reference and causes
        # stack overflow or kernel death
        # self.open_state = False
        self.open = True
        
        self._standby = False
        
        self._rf2_standby = True
        
        self._output = self.output = False
        
        self.freq = 1e9
        
        self._start_freq = 777e6
        self._stop_freq = 789e6
        self._step_freq = 1e6
        
        self._cycle_count = 0
        
        self._rf2_freq = 0
        
        self._level = self.level = 3.0
        self._auto_level_disable = self.auto_level_disable = False # setting this to 1 will lead to unstable output power
        
        self._rf_mode = self.rf_mode = "Single"
        
        self.phase_v = Phase_t(0) # stores ctype dict for processing
        self._phase = 0.0 # stores phase value, same for temperature
        
        self.temperature = Device_temperature_t(0)
        self._temp = 0
        
        self.status = Operate_status_t(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
        self.pll_status = Pll_status_t()
        self.dac_value = Dac_value_t(0)
        self.reference_dac = Dac_value_t(0)
        self.rf_params = Device_rf_params_t(0, 0, 0, 0, 0, 0, 0, 0, 0)
        self.reg_read = Reg_read_t(0)
        self.list_mode = List_mode_t()
        self.clock_config = Clock_config_t(0, 0, 0, 0, 0)
        self.manufacture_date = Manufacture_date_t(0, 0, 0, 0)
        self.device_info = Device_info_t(0, 0, 0, self.manufacture_date)
        self.device_status = Device_status_t(self.list_mode, self.status, self.pll_status)
        
        if debug:
            print(serial_number, self.handle)
            self.dll.sc5511a_get_device_status(self.handle, ctypes.byref(self.device_status))
            status = self.device_status.operate_status_t.rf1_out_enable
            print('check status', status)
        
        self.init = False
        
    
    '''
    @property
    def open(self):
        if not self.init:
            if self.open_state == True:
                print("Device session is active.")
            else:
                print("There are no active device sessions.")
        return self.open
    
    @open.setter
    def open(self, new_open):
        self.set_open(new_open)
        if not self.init:
            print(f"Set output boolean to {new_open}.")
    '''

    def set_open(self, opened):
        """opens a device session"""
        if opened and not self.open:
            self.handle = ctypes.c_void_p(self.dll.sc5511a_open_device(self.serial_number))
            self.open = True
            print("Opened device session.")
        elif not opened and self.open:
            self.dll.sc5511a_close_device(self.handle)
            self.open = False
            print("Closed device session.")
        return True

    def _close(self):
        """closes the device session"""
        self.set_open(0)
        
    @property
    def phase(self):
        error_code, p = self.get_signal_phase()
        if not self.init:
            print("Phase is {status} {error}".format(
                status = p,
                error = f"WITH ERROR CODE: :  {error_code}" if error_code else "(no errors)"))
    
    @phase.setter
    def phase(self, new_output):
        error_code, p = self.set_signal_phase(new_output)
        if not self.init:
            print("Set phase to {status} {error}".format(
                status = p, 
                error = f"WITH ERROR CODE: :  {error_code}" if error_code else "(no errors)"))
        
    def get_signal_phase(self):
        error_code = self.dll.sc5511a_get_signal_phase(self.handle, ctypes.byref(self.phase_v))
        signal_phase = self.phase_v.phase
        return error_code, signal_phase

    def set_signal_phase(self, phase = 0.0):
        """sets the signal phase of the device
        input: float"""
        c_phase= ctypes.c_float(phase)
        error_code = self.dll.sc5511a_set_signal_phase(self.handle, c_phase)
        self._phase = phase
        return error_code, phase

    @property
    def standby(self):
        print("Standby status is {status}.".format(
            status = "active" if self._standby else "not active"))
    
    @standby.setter
    def standby(self, new_standby):
        error_code, enable = self.set_standby(new_standby)
        self._standby = enable
        if not self.init:
            print("Set standby status to {status} {error}".format(
                status = "active" if self._standby else "not active", 
                error = f"WITH ERROR CODE: :  {error_code}" if error_code else "(no errors)"))

    def set_standby(self, enable = False):
        """sets the standby of the device.
        input: 
        standby (bool) true = standby enabled, false  = standby disabled."""
        error_code = self.dll.sc5511a_set_standby(self.handle, enable)
        return error_code, enable
    
    @property
    def rf2_standby(self):
        print("RF2 standby status is {status}.".format(
            status = "active" if self._standby else "not active"))
    
    @rf2_standby.setter
    def rf2_standby(self, new_standby):
        error_code, enable = self.set_rf2_standby(new_standby)
        self._rf2_standby = enable
        if not self.init:
            print("Set standby status to {status} {error}".format(
                status = "active" if self._rf2_standby else "not active", 
                error = f"WITH ERROR CODE: :  {error_code}" if error_code else "(no errors)"))

    def set_rf2_standby(self, rf2_enable=False):
        """Sets the standby of the device for rf2.
        input:
        standby (bool) true = standby enabled, false = standby disabled."""
        error_code = self.dll.sc5511a_set_rf2_standby(self.handle, rf2_enable)
        return error_code, rf2_enable

    @property
    def output(self):
        error_code, status = self.get_output()
        if not self.init:
            print("Device output status is {status} {error}".format(
                status = "active" if status else "not active", 
                error = f"WITH ERROR CODE: :  {error_code}" if error_code else "(no errors)"))
    
    @output.setter
    def output(self, new_output):
        error_code, output_enable = self.set_output(new_output)
        if not self.init:
            self._output = output_enable
            print("Set device output to {status} {error}".format(
                status = "active" if output_enable else "not active", 
                error = f"WITH ERROR CODE: :  {error_code}" if error_code else "(no errors)"))

    def set_output(self, output_enable):
        """
        Turns the output of RF1 on or off.
            Input:
                enable (int) = OFF = 0 ; ON = 1
        """
        c_enable = ctypes.c_ubyte(output_enable)
        error_code = self.dll.sc5511a_set_output(self.handle, c_enable)
        return error_code, output_enable
    

    def get_output(self):
        '''
        Reads the output status of RF1
            Output:
                status (int) : OFF = 0 ; ON = 1
        '''
        error_code = self.dll.sc5511a_get_device_status(self.handle, ctypes.byref(self.device_status))
        status = self.device_status.operate_status_t.rf1_out_enable
        return error_code, status

    @property
    def rf_mode(self):
        error_code, status = self.get_rf_mode()
        if not self.init:
            print("RF Mode is in {status} Mode, {error}".format(
                status = "Sweep" if status else "Single", 
                error = f"WITH ERROR CODE: :  {error_code}" if error_code else "(no errors)"))
    
    @rf_mode.setter
    def rf_mode(self, new_mode):
        if "single" in new_mode.lower():
            mode = 0
        elif "sweep" in new_mode.lower():
            mode = 1
        else:
            print("Mode not set. Available modes are Single and Sweep.")
            return
            
        error_code, output_enable = self.set_rf_mode(mode)
        if not self.init:
            self._rf_mode = new_mode
            print("Set RF Mode to {status} {error}".format(
                status = "Sweep Mode" if output_enable else "Single Mode", 
                error = f"WITH ERROR CODE: :  {error_code}" if error_code else "(no errors)"))

    def get_rf_mode(self):
        """
        Gets RF mode
        (int) = 0 = single tone, 1 = List/Sweep"""
        error_code = self.dll.sc5511a_get_device_status(self.handle, ctypes.byref(self.device_status))
        rf_mode = self.device_status.operate_status_t.rf1_mode
        return error_code, rf_mode
    
    def set_rf_mode(self, rf_mode = 0):
        """
        sets the rf mode of the device
        """
        error_code = self.dll.sc5511a_set_rf_mode(self.handle, rf_mode)
        return error_code, rf_mode
     
    @property
    def start_freq(self):
        print(f"Sweep starting frequency is {self._start_freq} Hz.") 
        
    @start_freq.setter
    def start_freq(self, start):
        error_code, start_freq = self.set_list_start_freq(int(start))
        if not self.init:
            print("Set starting frequency to {status} {error}".format(status = start,
                error = f"WITH ERROR CODE: :  {error_code}" if error_code else "(no errors)"))
     
    def set_list_start_freq(self, start_freq):
        """frequency in Hz"""
        frequency=ctypes.c_ulonglong(start_freq)
        error_code = self.dll.sc5511a_list_start_freq(self.handle, frequency)
        return error_code, start_freq
    
    @property
    def stop_freq(self):
        print(f"Sweep stop frequency is {self._stop_freq} Hz.") 
        
    @stop_freq.setter
    def stop_freq(self, stop):
        error_code, stop_freq = self.set_list_stop_freq(int(stop))
        if not self.init:
            print("Set stopping frequency to {status} {error}".format(status = stop,
                error = f"WITH ERROR CODE:  {error_code}" if error_code else "(no errors)"))

    def set_list_stop_freq(self, stop_freq):
        """Sets the list stop frequency
        input: int"""
        frequency = ctypes.c_ulonglong(stop_freq)
        error_code = self.dll.sc5511a_list_stop_freq(self.handle, frequency)
        return error_code, stop_freq
    
    @property
    def step_freq(self):
        print(f"Sweep frequency step is {self._step_freq} Hz.") 
        
    @step_freq.setter
    def step_freq(self, step):
        error_code, step_freq = self.set_list_step_freq(int(step))
        if not self.init:
            print("Set frequency step to {status} {error}".format(status = step,
                error = f"WITH ERROR CODE:  {error_code}" if error_code else "(no errors)"))

    def set_list_step_freq(self, step_freq):
        """Sets the list step frequency
        input: int"""
        frequency = ctypes.c_ulonglong(step_freq)
        error_code = self.dll.sc5511a_list_step_freq(self.handle, frequency)
        return error_code, step_freq
    
    @property
    def cycle_count(self):
        print(f"Sweep cycle count is {self._cycle_count}.") 
        
    @cycle_count.setter
    def cycle_count(self, count):
        error_code, cycles = self.set_list_cycle_count(int(count))
        if not self.init:
            print("Set number of sweep cycles to {status} {error}".format(status = cycles,
                error = f"WITH ERROR CODE:  {error_code}" if error_code else "(no errors)"))
    
    def set_list_cycle_count(self, cycle_num = 0):
        """Sets the list cycle count value
        input: integer
        """
        cycle_count = ctypes.c_uint(cycle_num)
        error_code = self.dll.sc5511a_list_cycle_count(self.handle, cycle_count)
        self.list_cycle_count = cycle_num
        return error_code, cycle_num

    def set_list_soft_trigger(self):
        """Sets the list soft trigger
        """
        error_code = self.dll.sc5511a_list_soft_trigger(self.handle)
        return error_code

    @property
    def freq(self):
        error_code, f = self.get_frequency()
        if not self.init:
            print("Frequency is {status} {error}".format(
                status = f, 
                error = f"WITH ERROR CODE:  {error_code}" if error_code else "(no errors)"))
    
    @freq.setter
    def freq(self, new_output):
        error_code, output_freq = self.set_frequency(new_output)
        if not self.init:
            print("Set frequency to {status} {error}".format(
                status = output_freq,
                error = f"WITH ERROR CODE:  {error_code}" if error_code else "(no errors)"))

    def get_frequency(self):
        """Gets the frequency value
        """
        error_code = self.dll.sc5511a_get_rf_parameters(self.handle, ctypes.byref(self.rf_params))
        frequency = self.rf_params.rf1_freq
        return error_code, frequency
    
    def set_frequency(self, frequency):
        """
        Sets RF1 frequency. Valid between 100MHz and 20GHz
            Args:
                frequency (int) = frequency in Hz
        """

        c_freq = ctypes.c_ulonglong(int(frequency))
        error_code = self.dll.sc5511a_set_freq(self.handle, c_freq)
        return error_code, frequency

    def set_clock_reference(self, ref, direct_lock, high_, lock_to_external):
        """Sets the clock reference
        input: (int) ..............................
        """
        ext_ref = ctypes.c_ubyte(ref)
        ext_direct_lock = ctypes.c_ubyte(direct_lock)
        high = ctypes.c_ubyte(high_)
        lock = ctypes.c_ubyte(lock_to_external)
        error_code = self.dll.sc5511a_set_clock_reference(self.handle, ext_ref, ext_direct_lock, high, lock)
        return error_code, ref, ext_direct_lock, high_, lock_to_external

    def get_reference_source(self):
        """ Gets the reference source
        ..........................
        """
        error_code =  self.device_status.operate_status_t.ext_ref_lock_enable
        return error_code

    @property
    def level(self):
        error_code, l = self.get_level()
        if not self.init:
            print("Power level is {status} {error}".format(
                status = l, 
                error = f"WITH ERROR CODE:  {error_code}" if error_code else "(no errors)"))
    
    @level.setter
    def level(self, new_output):
        error_code, output_level = self.set_level(new_output)
        if not self.init:
            print("Set power level to {status} {error}".format(
                status = output_level,
                error = f"WITH ERROR CODE: :  {error_code}" if error_code else "(no errors)"))

    def set_level(self, power):
        """Sets the power level
        input: float
        """
        c_power = ctypes.c_float(power)
        error_code = self.dll.sc5511a_set_level(self.handle, c_power)
        self._level = c_power
        return error_code, power

    def get_level(self):
        """Gets the power level of the device
        """
        error_code = self.dll.sc5511a_get_rf_parameters(self.handle, ctypes.byref(self.rf_params))
        rf_level = self.rf_params.rf_level
        return error_code, rf_level
    
    @property
    def auto_level(self):
        error_code, enable = self.get_auto_level_disable()
        if not self.init:
            print("Auto leveling is {status} {error}".format(
                status = "enabled" if enable else "disabled", 
                error = f"WITH ERROR CODE:  {error_code}" if error_code else "(no errors)"))
    
    @auto_level.setter
    def auto_level(self, new_output):
        error_code, enable = self.set_auto_level_disable(new_output)
        if not self.init:
            print("Auto leveling has been {status} {error}".format(
                status = "enabled" if enable else "disabled", 
                error = f"WITH ERROR CODE: :  {error_code}" if error_code else "(no errors)"))

    def set_auto_level_disable(self, enable):
        """Sets the auto level to either enable or disable
        """
        if enable == 1:
            enable = 0
        elif enable == 0:
            enable = 1
        c_enable = ctypes.c_ubyte(enable)
        error_code = self.dll.sc5511a_set_auto_level_disable(self.handle, c_enable)
        self._auto_level_disable = c_enable
        return error_code, enable

    def get_auto_level_disable(self):
        """Gets the current status of the auto level
        """
        error_code = self.dll.sc5511a_get_device_status(self.handle, ctypes.byref(self.device_status))
        enabled = self.device_status.operate_status_t.auto_pwr_disable
        if enabled == 1:
            enabled = 0
        elif enabled == 0:
            enabled = 1
        return error_code, enabled

    def get_ext_ref_lock_enable(self):
        #This function works only for devices with firmware >= ver3.6 and hardware > ver16.0
        error_code = self.dll.sc5511a_get_clock_config(self.handle, ctypes.byref(self.clock_config))
        ext_ref_lock = self.clock_config.ext_ref_lock_enable
        return error_code, ext_ref_lock

    def get_ref_out_select(self):
        #This function works only for devices with firmware >= ver3.6 and hardware > ver16.0
        error_code = self.dll.sc5511a_get_clock_config(self.handle, ctypes.byref(self.clock_config))
        ref_out_select = self.clock_config.ref_out_select
        return error_code, ref_out_select

    def get_pxi_clock_enable(self):
        #This function works only for devices with firmware >= ver3.6 and hardware > ver16.0
        error_code = self.dll.sc5511a_get_clock_config(self.handle, ctypes.byref(self.clock_config))
        pxi_clock_enable = self.clock_config.pxi_clock_enable
        return error_code, pxi_clock_enable

    def get_ext_direct_clock(self):
        #This function works only for devices with firmware >= ver3.6 and hardware > ver16.0
        error_code = self.dll.sc5511a_get_clock_config(self.handle, ctypes.byref(self.clock_config))
        ext_direct_clock = self.clock_config.ext_direct_clocking
        return error_code, ext_direct_clock

    def get_ext_ref_freq(self):
        #This function works only for devices with firmware >= ver3.6 and hardware > ver16.0
        error_code = self.dll.sc5511a_get_clock_config(self.handle, ctypes.byref(self.clock_config))
        ext_ref_freq = self.clock_config.ext_ref_freq
        return error_code, ext_ref_freq

    def get_alc_dac(self) -> int:
        """Gets the alc dac value
        """
        error_code = self.dll.sc5511a_get_alc_dac(self.handle, ctypes.byref(self.dac_value))
        dac_value = self.dac_value.dac_value
        return error_code, dac_value

    def store_default_state(self):
        """Stores the default state of the device
        """
        error_code = self.dll.sc5511a_store_default_state(self.handle)
        return error_code

    @property
    def temp(self):
        error_code, device_temp = self.get_temperature()
        self._temp = device_temp
        print("Device temp is {status} degrees Celsius {error}".format(
            status = device_temp,
            error = f"WITH ERROR CODE: :  {error_code}" if error_code else "(no errors)")) 
        return device_temp
        

    def get_temperature(self):
        """Gets the current temperature of the device
        """
        error_code = self.dll.sc5511a_get_temperature(self.handle, ctypes.byref(self.temperature))
        device_temp = self.temperature.device_temp
        return error_code, device_temp

    def set_list_mode(self, sss_mode=1, sweep_dir=0, tri_waveform=0, hw_trigger=0, step_on_hw_trig=0, return_to_start=0, trig_out_enable=1, trig_out_on_cycle=0) -> Dict[str, Optional[str]]:
        """
        Configures list mode
        sss_mode= 0: List mode, 1:sweep mode
        sweep_dir= 0:Forward, 1: reverse
        tri_waveform= 0:Sawtooth waveform, 1: Triangular waveform
        hw_trigger= 0:Software trigger, 1: Hardware trigger
        step_on_hw_trig= 0:Start/stop behavior, 1: Step on trigger, see manual for more details
        return_to_start = 0:stop at end of sweep/list, 1:return to start
        trig_out_enable= 0:No output trigger, 1: Output trigger enabled on trigger pin
        trig_out_on_cycle= 0: puts out a trigger pulse at each frequency change, 1: trigger pulse at the completion of each sweep/list cycle
        """
        lm = List_mode_t(sss_mode=sss_mode, sweep_dir=sweep_dir, tri_waveform=tri_waveform, hw_trigger=hw_trigger, step_on_hw_trig=step_on_hw_trig, return_to_start=return_to_start, trig_out_enable=trig_out_enable, trig_out_on_cycle=trig_out_on_cycle)
        error_code = self.dll.sc5511a_list_mode_config(self.handle, ctypes.byref(lm))
        self.list_mode = lm
        return error_code, lm

    def get_idn(self) -> Dict[str, Optional[str]]:
        """Gets the device identity
        """
        error_code = self.dll.sc5511a_get_device_info(self.handle, ctypes.byref(self.device_info))
        #device_info = self.device_info

        IDN: Dict[str, Optional[str]] = {
            'vendor': "SignalCore",
            'model': "SC5511A",
            'serial_number': self.serial_number.value.decode("utf-8"),
            'firmware_revision': self.device_info.firmware_revision,
            'hardware_revision': self.device_info.hardware_revision,
            'manufacture_date': '20{}={}-{}, {}:00'.format(self.device_info.manufacture_date.year, self.device_info.manufacture_date.month, self.device_info.manufacture_date.day, self.device_info.manufacture_date.hour)
            }
        self.device_info = IDN
        return error_code, IDN
    
    def settings(self):
        """Returns all of the parameters"""
        error_code = self.dll.sc5511a_get_rf_parameters(self.handle, ctypes.byref(self.rf_params))
        if not error_code:
            print("Returned parameters successfully:")
        else:
            print(f"Failed to return parameters with error code {error_code}.")
        params = self.rf_params

        RF: Dict[str, Optional[str]] = {
            'output' : str(self._output),
            'rf_mode' : self._rf_mode,
            'freq' : params.rf1_freq,
            'phase' : str(self._phase),
            'start_freq' : params.start_freq,
            'stop_freq' : params.stop_freq,
            'step_freq' : params.step_freq,
            #'sweep_dwell_time' : params.sweep_dwell_time,
            'cycle_count' : params.sweep_cycles,
            #'buffer_points' : params.buffer_points,
            'level' : params.rf_level,
            'auto_level' : str(self._auto_level_disable),
            'standby' : str(self._standby),
            'rf2_standby' : str(self._rf2_standby),
            'rf2_freq' : params.rf2_freq,
            'temp' : str(self.temp),
            }

        return RF

    def get_clock_config(self)->Dict[str, Optional[str]]:
        #This function works only for devices with firmware >= ver3.6 and hardware > ver16.0
        error_code = self.dll.sc5511a_get_clock_config(self.handle, ctypes.byref(self.clock_config))
        clock_config = self.clock_config

        CLOCK: Dict[str, Optional[str]] = {
            'ext_ref_lock_enable' : clock_config.clock_config,
            'ref_out_select' : clock_config.clock_config,
            'pxi_clock_enable' : clock_config.clock_config,
            'ext_direct_clocking' : clock_config.clock_config,
            'ext_ref_freq' : clock_config.clock_config}

        return error_code, CLOCK

    def get_list_mode(self)-> Dict[str, Optional[str]]:
        """Gets the current list mode values
        """
        error_code = self.dll.sc5511a_get_device_status(self.handle, ctypes.byref(self.device_status))
        device_status = self.device_status

        LIST: Dict[str, Optional[str]] = {
            'vendor': "SignalCore",
            'model': "SC5511A",
            'sss_mode' : device_status.list_mode.sss_mode,
            'sweep_dir': device_status.list_mode.sweep_dir,
            'tri_waveform': device_status.list_mode.tri_waveform,
            'hw_trigger': device_status.list_mode.hw_trigger,
            'step_on_hw_trig': device_status.list_mode.step_on_hw_trig,
            'return_to_start': device_status.list_mode.return_to_start,
            'trig_out_enable': device_status.list_mode.trig_out_enable,
            'trig_out_on_cycle': device_status.list_mode.trig_out_on_cycle
            }
        return error_code, LIST

    def get_operate_status(self)-> Dict[str, Optional[str]]:
        """Gets the current operate status values
        """
        error_code = self.dll.sc5511a_get_device_status(self.handle, ctypes.byref(self.device_status))
        device_status = self.device_status

        OPERATE: Dict[str, Optional[str]] = {
            'vendor' : "SignalCore",
            'model' : "SC5511A",
            'rf1_lock_mode' : device_status.operate_status_t.rf1_lock_mode,
            'rf1_loop_gain' : device_status.operate_status_t.rf1_loop_gain,
            'device_access' : device_status.operate_status_t.device_access,
            'rf2_standby' : device_status.operate_status_t.rf2_standby,
            'rf1_standby' : device_status.operate_status_t.rf1_standby,
            'auto_pwr_disable' : device_status.operate_status_t.auto_pwr_disable,
            'alc_mode' : device_status.operate_status_t.alc_mode,
            'rf1_out_enable' : device_status.operate_status_t.rf1_out_enable,
            'ext_ref_lock_enable' : device_status.operate_status_t.ext_ref_lock_enable,
            'ext_ref_detect' : device_status.operate_status_t.ext_ref_detect,
            'ref_out_select' : device_status.operate_status_t.ref_out_select,
            'list_mode_running' : device_status.operate_status_t.list_mode_running,
            'rf1_mode' : device_status.operate_status_t.rf1_mode,
            'over_temp' : device_status.operate_status_t.over_temp,
            'harmonic_ss' : device_status.operate_status_t.harmonic_ss}
        return error_code, OPERATE

    def get_pll_status(self)-> Dict[str, Optional[str]]:
        """Gets the current PLL status values
        """
        error_code = self.dll.sc5511a_get_device_status(self.handle, ctypes.byref(self.device_status))
        device_status = self.device_status
        #_pll_status_t = self.pll_status
        PLL: Dict[str, Optional[str]] = {
            'vendor' : "SignalCore",
            'model' : "SC5511A",
            'sum_pll_ld' : device_status.pll_status_t.sum_pll_ld,
            'crs_pll_ld' : device_status.pll_status_t.crs_pll_ld,
            'fine_pll_ld' : device_status.pll_status_t.fine_pll_ld,
            'crs_ref_pll_ld' : device_status.pll_status_t.crs_ref_pll_ld,
            'crs_aux_pll_ld' : device_status.pll_status_t.crs_aux_pll_ld,
            'ref_100_pll_ld' : device_status.pll_status_t.ref_100_pll_ld,
            'ref_10_pll_ld' : device_status.pll_status_t.ref_10_pll_ld,
            'rf2_pll_ld' : device_status.pll_status_t.rf2_pll_ld}

        return error_code, PLL 


    def reg_write(self, reg_byte, ins_word = int):
        """Writes the specific value to the register, see documentation for more info
        """
        instruct_word = ctypes.c_uint(ins_word)
        error_code = self.dll.sc5511a_reg_write(self.handle, reg_byte, instruct_word)
        # not sure if I should set reg_read variable to reg_byte or ins_word or... both?
        # don't think this function will be used in normal operation
        return error_code, ins_word

    def reg_read(self, reg_byte, ins_word):
        """Reads back the specific value back, see documentation for more info
        """
        reg_byte = ctypes.c_ubyte(reg_byte)
        instruct_word = ctypes.c_ulonglong(ins_word)
        error_code = self.dll.sc5511a_reg_read(self.handle, reg_byte, instruct_word, ctypes.byref(self.reg_read))
        rec_word = self.reg_read.reg_read

        return error_code, rec_word

    def set_synth_mode(self, disable_spur_suppress = 0, low_loop_gain = 0, lock_mode = 0):
        """
        sets the rf mode of the device
        Disable spur suppress: (only takes effect when lock mode is harmonic)
            input: integer = 0 = Normal gain, 1 = spur suppress by lowering loop gain and/or ping ponging between lock modes automatically
        Low loop gain:
            input: integer = 0 = Normal gain, 1 = low gain
        Lock mode:
            input: integer = 0 = harmonic, 1 = fractional
        """
        error_code = self.dll.sc5511a_set_synth_mode(self.handle, disable_spur_suppress, low_loop_gain, lock_mode)
        return error_code, disable_spur_suppress, low_loop_gain, lock_mode

    def set_list_dwell_time(self, tunit):
        """ 1: 500 us, 2: 1ms"""
        error_code = self.dll.sc5511a_list_dwell_time(self.handle, tunit)
        #time = tunit*500
        return error_code, tunit

    def set_alc_mode(self, alc_mode = 0):
        """Sets the ALC to close(0) or open (1) mode operation for channel RF1"""
        error_code = self.dll.sc5511a_set_alc_mode(self.handle, alc_mode)
        self.alc_mode = alc_mode
        return error_code, alc_mode

    def set_reference_dac(self, d_value = 0):
        """Sets a value to the reference dac
        input: integer
        """
        dac_value = ctypes.c_uint(d_value)
        error_code = self.dll.sc5511a_set_reference_dac(self.handle, dac_value)
        self.reference_dac.dac_value = d_value
        return error_code, d_value

    def set_alc_dac(self, alc_value = 0):
        """Sets the ALC DAC value
        input: integer
        """
        dac_value = ctypes.c_uint(alc_value)
        error_code = self.dll.sc5511a_set_alc_dac(self.handle, dac_value)
        self.dac_value.dac_value = alc_value
        return error_code, alc_value

    def list_buffer_points(self, l_points = 0):
        """Sets the list buffer points
        input: integer
        """
        list_points = ctypes.c_uint(l_points)
        error_code = self.dll.sc5511a_list_buffer_points(self.handle, list_points)
        return error_code, l_points

    def list_buffer_write(self, command):
        """Writes to the list buffer
        input: HEX or float values for frequency
        """
        command1 = ctypes.c_ulonglong(command)
        error_code = self.dll.sc5511a_list_buffer_write(self.handle, command1)
        return error_code, command
        

    def list_buffer_transfer(self, transfer_mode = 0):
        """Transfers values to and from the EEPROM/buffer
        """
        transfer = ctypes.c_ubyte(transfer_mode)
        error_code = self.dll.sc5511a_list_buffer_transfer(self.handle, transfer)
        return error_code, transfer_mode


    def list_buffer_read(self, address=0) -> int:
        """Reads back values from the list buffer
        """
        address = ctypes.c_uint(address)
        error_code = self.dll.sc5511a_list_buffer_read(self.handle, address, ctypes.byref(self.freq))
        freq = self.freq.device_freq
        return error_code, address, freq

    @property
    def rf2_freq(self):
        print(f"RF2 frequency is {self._rf2_freq} Hz.") 
        
    @rf2_freq.setter
    def rf2_freq(self, f):
        error_code, output_freq = self.set_rf2_frequency(int(f))
        if not self.init:
            print("Set RF2 frequency to {status} {error}".format(status = f,
                error = f"WITH ERROR CODE:  {error_code}" if error_code else "(no errors)"))

    def set_rf2_frequency(self, freq):
        """Sets the RF2 frequency value
        input: integer
        """
        error_code = self.dll.sc5511a_set_rf2_freq(self.handle, freq)
        return error_code, freq

    def synth_self_cal(self):
        """Calibrates the synthesizer of the device
        """
        error_code = self.dll.sc5511a_synth_self_cal(self.handle)
        return error_code