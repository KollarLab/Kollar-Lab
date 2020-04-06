# -*- coding: utf-8 -*-
import ctypes

class AqMd3(object):
    def __init__(self, library, ResourceName):        
        self._prefix = "AqMD3"
        self._lib = ctypes.cdll.LoadLibrary(library)
              
        self.visession = ctypes.c_int()
        
        self.call('InitWithOptions',
            ctypes.create_string_buffer(ResourceName.encode('utf-8')),
            False,
            False, 
            '',
            ctypes.byref(self.visession))

    def SetAttribute(self, RepCapIdentifier, Attribute, Value, AttributeType) :

        # Encode the RepCapIdentifier if needed
        if (RepCapIdentifier != None):
            RepCapIdentifier = ctypes.c_char_p(RepCapIdentifier.encode('utf-8'))
                   
        if (AttributeType == 'ViBoolean') :            
            self.call('SetAttributeViBoolean',
                  self.visession,
                  RepCapIdentifier,
                  Attribute,
                  ctypes.c_bool(Value))
            
        if (AttributeType == 'ViInt32') :
            self.call('SetAttributeViInt32',
                  self.visession,
                  RepCapIdentifier,
                  Attribute,
                 ctypes.c_int32(Value))
        
        if (AttributeType == 'ViInt64') :            
            self.call('SetAttributeViInt64',
                  self.visession,
                  RepCapIdentifier,
                  Attribute,
                  ctypes.c_int64(Value))

        if (AttributeType == 'ViReal64') :
            self.call('SetAttributeViReal64',
                  self.visession,
                  RepCapIdentifier,
                  Attribute,
                  ctypes.c_double(Value))
            
        if (AttributeType == 'ViString') :
            buffer = ctypes.create_string_buffer(256)
            self.call('SetAttributeViString',
                  self.visession,
                  RepCapIdentifier,
                  Attribute,
                  ctypes.c_char_p(Value.encode('utf-8')))

    def GetAttribute(self, RepCapIdentifier, Attribute, AttributeType) :
        
        # Encode the RepCapIdentifier if needed
        if (RepCapIdentifier != None):
            RepCapIdentifier = ctypes.c_char_p(RepCapIdentifier.encode('utf-8'))
                   
        if (AttributeType == 'ViBoolean') :
            buffer = ctypes.c_bool();
            self.call('GetAttributeViBoolean',
                  self.visession,
                  RepCapIdentifier,
                  Attribute, 
                  ctypes.byref(buffer))
            return buffer.value      
            
        if (AttributeType == 'ViInt32') :
            buffer = ctypes.c_int32()
            self.call('GetAttributeViInt32',
                  self.visession,
                  RepCapIdentifier,
                  Attribute, 
                  ctypes.byref(buffer))
            return buffer.value
        
        if (AttributeType == 'ViInt64') :
            buffer = ctypes.c_int64()
            self.call('GetAttributeViInt64',
                  self.visession,
                  RepCapIdentifier,
                  Attribute, 
                  ctypes.byref(buffer))
            return buffer.value

        if (AttributeType == 'ViReal64') :
            buffer = ctypes.c_double()
            self.call('GetAttributeViReal64',
                  self.visession,
                  RepCapIdentifier,
                  Attribute, 
                  ctypes.byref(buffer))
            return buffer.value
            
        if (AttributeType == 'ViString') :
            buffer = ctypes.create_string_buffer(256)
            self.call('GetAttributeViString',
                  self.visession,
                  RepCapIdentifier,
                  Attribute, 
                  256,
                  ctypes.byref(buffer))
            return buffer.value.decode('ascii')        
           
        return ''

    def call(self, funcname, *args):
        method_to_call = getattr(self._lib, self._prefix + '_' + funcname)
        return_value = method_to_call(*args)        
        if return_value:
            error_message = ctypes.create_string_buffer(256)
            self.call('error_message',
                      self.visession,
                      return_value, 
                      ctypes.byref(error_message))
            raise IviCError(error_message.value)        
            
            
    def close(self):    
        self.call('close', self.visession) 
        
    def SelfCalibrate(self):    
        self.call('SelfCalibrate', self.visession) 

if __name__ == '__main__':
    
    ResourceName = "PXI2::0::0::INSTR"
    
    print("Initialize instrument @", ResourceName)
    drv = AqMd3("AqMd3_64.dll", ResourceName)
    
    # Declaration of some attributes (values from AqMd3.h)
    IVI_ATTR_BASE          = 1000000
    IVI_INHERENT_ATTR_BASE = IVI_ATTR_BASE +  50000
    IVI_CLASS_ATTR_BASE    = IVI_ATTR_BASE +  250000
    IVI_LXISYNC_ATTR_BASE  = IVI_ATTR_BASE +  950000
    IVI_SPECIFIC_ATTR_BASE = IVI_ATTR_BASE +  150000
        
    AQMD3_ATTR_SPECIFIC_DRIVER_DESCRIPTION          = IVI_INHERENT_ATTR_BASE + 514  # ViString, read-only 
    AQMD3_ATTR_SPECIFIC_DRIVER_REVISION             = IVI_INHERENT_ATTR_BASE + 551  # ViString, read-only   
    AQMD3_ATTR_INSTRUMENT_FIRMWARE_REVISION         = IVI_INHERENT_ATTR_BASE + 510  # ViString, read-only
    AQMD3_ATTR_INSTRUMENT_MANUFACTURER              = IVI_INHERENT_ATTR_BASE + 511  # ViString, read-only
    AQMD3_ATTR_INSTRUMENT_MODEL                     = IVI_INHERENT_ATTR_BASE + 512  # ViString, read-only
    AQMD3_ATTR_INSTRUMENT_INFO_SERIAL_NUMBER_STRING = IVI_SPECIFIC_ATTR_BASE + 8    # ViString, read-only 
    
    AQMD3_ATTR_VERTICAL_OFFSET                      = IVI_CLASS_ATTR_BASE + 25      # ViReal64, read-write 
    
    AQMD3_ATTR_ACTIVE_TRIGGER_SOURCE                = IVI_CLASS_ATTR_BASE + 1       # ViString, read-write 
    AQMD3_ATTR_ACQUISITION_MODE                     = IVI_SPECIFIC_ATTR_BASE + 11   # ViInt32, read-write 
   
    # Read some attributes
    print('Manufacturer  : {}'.format(drv.GetAttribute(None, AQMD3_ATTR_INSTRUMENT_MANUFACTURER, 'ViString')))
    print('Description   : {}'.format(drv.GetAttribute(None, AQMD3_ATTR_SPECIFIC_DRIVER_DESCRIPTION, 'ViString')))
    print('Revision      : {}'.format(drv.GetAttribute(None, AQMD3_ATTR_SPECIFIC_DRIVER_REVISION, 'ViString')))
    print('Model         : {}'.format(drv.GetAttribute(None, AQMD3_ATTR_INSTRUMENT_MODEL, 'ViString')))
    print('Firmware Rev. : {}'.format(drv.GetAttribute(None, AQMD3_ATTR_INSTRUMENT_FIRMWARE_REVISION, 'ViString')))
    print('Serial #      : {}'.format(drv.GetAttribute(None, AQMD3_ATTR_INSTRUMENT_INFO_SERIAL_NUMBER_STRING, 'ViString')))
    
    # Read Channel1 offset
    print('\nRead Channel1 Offset: {}'.format(drv.GetAttribute('Channel1', AQMD3_ATTR_VERTICAL_OFFSET, 'ViReal64')))
    
    # Set Channel1 offset to 0.5V
    print('\nChange Channel1 offset to 0.05 V')
    drv.SetAttribute('Channel1', AQMD3_ATTR_VERTICAL_OFFSET, 0.05, 'ViReal64')

    # Read Channel1 offset
    print('\nRead Channel1 Offset: {}'.format(drv.GetAttribute('Channel1', AQMD3_ATTR_VERTICAL_OFFSET, 'ViReal64')))    
    
    # Read Channel1 offset
    print('\nRead Channel1 Offset: {}'.format(drv.GetAttribute('Channel1', AQMD3_ATTR_VERTICAL_OFFSET, 'ViReal64')))    
    
    # Read active trigger source
    print('\nRead active trigger source : {}'.format(drv.GetAttribute(None, AQMD3_ATTR_ACTIVE_TRIGGER_SOURCE, 'ViString')))    

    # Set active trigger source to External1
    print('\nSet active trigger source to External1')
    drv.SetAttribute(None, AQMD3_ATTR_ACTIVE_TRIGGER_SOURCE, 'External1', 'ViString')
    
    # Read active trigger source
    print('\nRead active trigger source : {}'.format(drv.GetAttribute(None, AQMD3_ATTR_ACTIVE_TRIGGER_SOURCE, 'ViString')))  

    # Set acquisition mode to averager
    print('\nSet acquisition mode to averager')
    drv.SetAttribute(None, AQMD3_ATTR_ACQUISITION_MODE, 1, 'ViInt32')   
 
    # Read acquisition mode
    print('\nRead acquisition mode: {}'.format(drv.GetAttribute(None, AQMD3_ATTR_ACQUISITION_MODE, 'ViInt32')))
    
     
    # Perform a self-calibration
    print("\nPerform internal calibrations")    
    drv.SelfCalibrate()
    
    # Close the instrument
    print("\nClose the instrument")    
    drv.close()