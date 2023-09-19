"""
Support functions.
"""
from typing import Union, List
import numpy as np
import json
import base64
from collections import OrderedDict


def gauss(mu=0, si=25, length=100, maxv=30000):
    """
    Create a numpy array containing a Gaussian function

    :param mu: Mu (peak offset) of Gaussian
    :type mu: float
    :param sigma: Sigma (standard deviation) of Gaussian
    :type sigma: float
    :param length: Length of array
    :type length: int
    :param maxv: Maximum amplitude of Gaussian
    :type maxv: float
    :return: Numpy array containing a Gaussian function
    :rtype: array
    """
    x = np.arange(0, length)
    y = maxv * np.exp(-(x-mu)**2/si**2)
    return y


def DRAG(mu, si, length, maxv, delta, alpha):
    """
    Create I and Q arrays for a DRAG pulse.
    Based on QubiC and Qiskit-Pulse implementations.

    :param mu: Mu (peak offset) of Gaussian
    :type mu: float
    :param si: Sigma (standard deviation) of Gaussian
    :type si: float
    :param length: Length of array
    :type length: int
    :param maxv: Maximum amplitude of Gaussian
    :type maxv: float
    :param delta: anharmonicity of the qubit (units of 1/sample time)
    :type delta: float
    :param alpha: alpha parameter of DRAG (order-1 scale factor)
    :type alpha: float
    :return: Numpy array with I and Q components of the DRAG pulse
    :rtype: array, array
    """
    x = np.arange(0, length)
    gaus = maxv * np.exp(-(x-mu)**2/si**2)
    # derivative of the gaussian
    dgaus = -(x-mu)/(si**2)*gaus
    idata = gaus
    qdata = -1 * alpha * dgaus / delta
    return idata, qdata


def triang(length=100, maxv=30000):
    """
    Create a numpy array containing a triangle function

    :param length: Length of array
    :type length: int
    :param maxv: Maximum amplitude of triangle function
    :type maxv: float
    :return: Numpy array containing a triangle function
    :rtype: array
    """
    y = np.zeros(length)

    # if length is even, there are length//2 samples in the ramp
    # if length is odd, there are length//2 + 1 samples in the ramp
    halflength = (length + 1) // 2

    y1 = np.linspace(0, maxv, halflength)
    y[:halflength] = y1
    y[length//2:length] = np.flip(y1)
    return y

class NpEncoder(json.JSONEncoder):
    """
    JSON encoder with support for numpy objects.
    Taken from https://stackoverflow.com/questions/50916422/python-typeerror-object-of-type-int64-is-not-json-serializable
    """
    def default(self, obj):
        if isinstance(obj, np.integer):
            return int(obj)
        if isinstance(obj, np.floating):
            return float(obj)
        if isinstance(obj, np.ndarray):
            # base64 is considerably more compact and faster to pack/unpack
            # return obj.tolist()
            return (base64.b64encode(obj.tobytes()).decode(), obj.shape, obj.dtype.str)
        return super().default(obj)

def progs2json(proglist):
    return json.dumps(proglist, cls=NpEncoder)

def json2progs(s):
    if hasattr(s, 'read'):
        # input is file-like, we should use json.load()
        # be sure to read dicts back in order (only matters for Python <3.7)
        proglist = json.load(s, object_pairs_hook=OrderedDict)
    else:
        # input is string or bytes
        # be sure to read dicts back in order (only matters for Python <3.7)
        proglist = json.loads(s, object_pairs_hook=OrderedDict)

    for progdict in proglist:
        # tweak data structures that got screwed up by JSON:
        # in JSON, dict keys are always strings, so we must cast back to int
        progdict['gen_chs'] = OrderedDict([(int(k),v) for k,v in progdict['gen_chs'].items()])
        progdict['ro_chs'] = OrderedDict([(int(k),v) for k,v in progdict['ro_chs'].items()])
        # the envelope arrays need to be restored as numpy arrays with the proper type
        for iCh, pulsedict in enumerate(progdict['pulses']):
            for name, pulse in pulsedict.items():
                #pulse['data'] = np.array(pulse['data'], dtype=self._gen_mgrs[iCh].env_dtype)
                data, shape, dtype = pulse['data']
                pulse['data'] = np.frombuffer(base64.b64decode(data), dtype=np.dtype(dtype)).reshape(shape)
    return proglist

class QickMetadata:
    """
    Provides information about the connections between IP blocks, extracted from the HWH file.
    The HWH parser is very different between PYNQ 2.6/2.7 and 3.0+, so this class serves as a common interface.
    """
    def __init__(self, soc):
        # We will use the HWH parser to extract information about signal connections between blocks.
        self.sigparser = None
        self.busparser = None
        self.systemgraph = None
        self.xml = None

        if hasattr(soc, 'systemgraph'):
            # PYNQ 3.0 and higher have a "system graph"
            self.systemgraph = soc.systemgraph
            # TODO: We shouldn't need to use BusParser, but we think there's a bug in how pynqmetadata handles axis_switch.
            self.busparser = BusParser(self.systemgraph._root)
            self.xml = soc.systemgraph._element_tree
        else:
            self.sigparser = soc.parser
            # Since the HWH parser doesn't parse buses, we also make our own BusParser.
            self.busparser = BusParser(self.sigparser.root)
            self.xml = soc.parser.root

    def trace_sig(self, blockname, portname):
        if self.systemgraph is not None:
            dests = self.systemgraph.blocks[blockname].ports[portname].destinations()
            result = []
            for port, block in dests.items():
                blockname = block.parent().name
                if blockname==self.systemgraph.name:
                    result.append([port])
                else:
                    result.append([blockname, port])
            return result

        return self._trace_net(self.sigparser, blockname, portname)

    def trace_bus(self, blockname, portname):
        return self._trace_net(self.busparser, blockname, portname)

    def _trace_net(self, parser, blockname, portname):
        """
        Find the block and port that connect to this block and port.
        If you expect to only get one block+port as a result, you can assign the result to ((block, port),)

        :param parser: HWH parser object (from Overlay.parser, or BusParser)
        :param blockname: the IP block of interest
        :type blockname: string
        :param portname: the port we want to trace
        :type portname: string

        :return: a list of [block, port] pairs, or just [port] for ports of the top-level design
        :rtype: list
        """
        fullport = blockname+"/"+portname
        # the net connected to this port
        netname = parser.pins[fullport]
        if netname == '__NOC__':
            return []
        # get the list of other ports on this net, discard the port we started at and ILA ports
        return [x.split('/') for x in parser.nets[netname] if x != fullport and 'system_ila_' not in x]

    def get_fclk(self, blockname, portname):
        """
        Find the frequency of a clock port.

        :param parser: HWH parser object (from Overlay.parser, or BusParser)
        :param blockname: the IP block of interest
        :type blockname: string
        :param portname: the port we want to trace
        :type portname: string

        :return: frequency in MHz
        :rtype: float
        """
        xmlpath = "./MODULES/MODULE[@FULLNAME='/{0}']/PORTS/PORT[@NAME='{1}']".format(
            blockname, portname)
        port = self.xml.find(xmlpath)
        return float(port.get('CLKFREQUENCY'))/1e6

    def get_param(self, blockname, parname):
        """
        Find the value of an IP parameter. This works for all IPs, including those that do not show up in ip_dict because they're not addressable.

        :param parser: HWH parser object (from Overlay.parser, or BusParser)
        :param blockname: the IP block of interest
        :type blockname: string
        :param parname: the parameter of interest
        :type parname: string

        :return: parameter value
        :rtype: string
        """
        xmlpath = "./MODULES/MODULE[@FULLNAME='/{0}']/PARAMETERS/PARAMETER[@NAME='{1}']".format(
            blockname, parname)
        param = self.xml.find(xmlpath)
        return param.get('VALUE')

    def mod2type(self, blockname):
        if self.systemgraph is not None:
            return self.systemgraph.blocks[blockname].vlnv.name
        return self.busparser.mod2type[blockname]

class BusParser:
    def __init__(self, root):
        """
        Matching all the buses in the modules from the HWH file.
        This is essentially a copy of the HWH parser's match_nets() and match_pins(),
        but working on buses instead of signals.

        In addition, there's a map from module names to module types.

        :param root: HWH XML tree (from Overlay.parser.root)
        """
        self.nets = {}
        self.pins = {}
        self.mod2type = {}
        for module in root.findall('./MODULES/MODULE'):
            fullpath = module.get('FULLNAME').lstrip('/')
            self.mod2type[fullpath] = module.get('MODTYPE')
            for bus in module.findall('./BUSINTERFACES/BUSINTERFACE'):
                port = fullpath + '/' + bus.get('NAME')
                busname = bus.get('BUSNAME')
                self.pins[port] = busname
                if busname in self.nets:
                    self.nets[busname] |= set([port])
                else:
                    self.nets[busname] = set([port])


def ch2list(ch: Union[List[int], int]) -> List[int]:
    """
    convert a channel number or a list of ch numbers to list of integers

    :param ch: channel number or list of channel numbers
    :return: list of channel number(s)
    """
    if ch is None:
        return []
    try:
        ch_list = [int(ch)]
    except TypeError:
        ch_list = ch
    return ch_list
