import sdxf
from math import sin,cos,pi,floor,asin,acos,tan,atan,sqrt
from alphanum import alphanum_dict
from random import randrange
from MaskMakerPro import *
from scipy.optimize import fsolve
import numpy
from sympy import nsolve, Symbol
from scipy import optimize
   
def angleDifference(a,b):
    return (a-b + 180) % 360 - 180

class getDrawDictionaries:
    """ Returns the angles of all the vertices and which vertices each resonator starts and stops at """
    def __init__(self,vertexDict,vertexPositions,resonators,roundDepth = 3):
        """ """
        self.vertexAngles = {} 
        self.capacitorAngles = {}
        self.verboseVertexDict = {}
        self.resonatorStartStopVertices = {}
        self.resonatorInputOuputAngles = {}

        self.vertexDict = vertexDict
        self.vertexPositions = vertexPositions
        self.resonators = resonators
    

        for vtex in range(len(vertexPositions)):
        # connectedResonators = vertexDict[vtex]
            connected = self.vertexDict[vtex]
            angles = []
            for res in connected:

                # find resonator positions
                resPos = resonators[res]
                startPos = [numpy.round(resPos[0],roundDepth), numpy.round(resPos[1],roundDepth)]
                stopPos = [numpy.round(resPos[2],roundDepth), numpy.round(resPos[3],roundDepth)]

                # find which vertex a resonator starts from
                if startPos[0] == vertexPositions[vtex][0] and startPos[1] == vertexPositions[vtex][1]:
                    deltaY = stopPos[1] - startPos[1]
                    deltaX = stopPos[0] - startPos[0]

                    if res in self.resonatorStartStopVertices.keys():
                        resSS = self.resonatorStartStopVertices[res]
                        self.resonatorStartStopVertices[res] = [vtex, resSS[1]]
                    else:
                        self.resonatorStartStopVertices[res] = [vtex, -1]

                # find which vertex a resonator ends at
                elif stopPos[0] == vertexPositions[vtex][0] and stopPos[1] == vertexPositions[vtex][1]:
                    deltaY = startPos[1] - stopPos[1]
                    deltaX = startPos[0] - stopPos[0]

                    if res in self.resonatorStartStopVertices.keys():
                        resSS = self.resonatorStartStopVertices[res]
                        self.resonatorStartStopVertices[res] = [resSS[0],vtex]
                    else:
                        self.resonatorStartStopVertices[res] = [-1,vtex]

                # find the angle of a given resonator from that vertex
                tempAngle = numpy.arctan2(deltaY,deltaX)*180/numpy.pi
                if tempAngle<0:
                    tempAngle = 360.+tempAngle
                angles.append(tempAngle)
                
                # vertexAngles now has a 1-1 correspondence with res    
                self.vertexAngles[vtex] = angles


                # now try some orientations and figure out the minimal offset orientation
                threewayInteriorAngle = 120
                minAngleOffset = 3300000. # some really big number
                if len(angles) == 1:
                    minThreewayAngles = [angles[0]-120.,angles[0],angles[0]+120.] 
                    self.verboseVertexDict[vtex] = [-1,res,-1]

                elif len(angles) == 2:
                    # symmetrically place the capacitors
                    meanAngle = ((angleDifference(angles[0],angles[1])/2 + angles[1])+180.)%360
                    possibleAngles = [(meanAngle-120)%360,(meanAngle+120)%360]
                    
                    # find the possible angle that is closest to angles[0]
                    minCapAngleIdx0 = numpy.argmin(abs(angles[0]-possibleAngles))
                    minCapAngle0 = possibleAngles[minCapAngleIdx0]

                    # find the possible angle that is closest to angles[1]
                    minCapAngleIdx1 = numpy.argmin(abs(angles[1]-possibleAngles))

                    # check to see if that angle has been used already
                    if minCapAngleIdx1 == minCapAngleIdx0:
                        possibleAngles.remove(possibleAngles[minCapAngleIdx1])
                        minCapAngleIdx1 = numpy.argmin(abs(angles[1]-possibleAngles))
                    minCapAngle1 = possibleAngles[minCapAngleIdx1]

                    # assign the appropriate angles to a temporary array
                    minThreewayAngles = [minCapAngle0, minCapAngle1, meanAngle] 

                    # create a verbose dictionary which has -1 for stub resonators
                    self.verboseVertexDict[vtex] = [connected[0],connected[1],-1]

                    # print('vtex',vtex,'meanAngle',meanAngle,'angles',angles,'minCapAngle0',minCapAngle0,'minCapAngle1',minCapAngle1,'minThreewayAngles',minThreewayAngles)
                
                else:
                    # try all orientations and figure out the optimal way to lay out capacitors
                    for odx in range(len(angles)):
                        orientation = angles[odx]
                        # possibleAngles = [orientation-240.,orientation-120,orientation+120,orientation+240]
                        possibleAngles = [(orientation-120)%360,(orientation+120)%360]
                        

                        if odx == 0:
                            # find the possible angle that is closest to angles[1],angles[2]
                            minCapAngleIdx0 = numpy.argmin(abs(angleDifference(angles[1],possibleAngles)))
                            minCapAngle0 = possibleAngles[minCapAngleIdx0]

                            possibleAngles.remove(minCapAngle0)
                            minCapAngle1 = possibleAngles[0]

                            # assign the appropriate angles to a temporary array
                            threewayAngles = [orientation, minCapAngle0, minCapAngle1]  
                        
                        elif odx == 1:
                            # find the possible angle that is closest to angles[0],angles[2]
                            minCapAngleIdx0 = numpy.argmin(abs(angleDifference(angles[0],possibleAngles)))
                            minCapAngle0 = possibleAngles[minCapAngleIdx0]

                            possibleAngles.remove(minCapAngle0)
                            minCapAngle1 = possibleAngles[0]

                            # assign the appropriate angles to a temporary array
                            threewayAngles = [minCapAngle0, orientation, minCapAngle1]  
                        
                        elif odx == 2:
                            # find the possible angle that is closest to angles[0],angles[1]
                            minCapAngleIdx0 = numpy.argmin(abs(angleDifference(angles[0],possibleAngles)))
                            minCapAngle0 = possibleAngles[minCapAngleIdx0]
                            possibleAngles.remove(minCapAngle0)
                            minCapAngle1 = possibleAngles[0]

                            # assign the appropriate angles for the vertex
                            threewayAngles = [minCapAngle0, minCapAngle1, orientation]  
                            # print('checking','vtex',vtex,'res',res,'orientation',angles[odx],'minCapAngle1',minCapAngle1)

                        # find the total offset angle for a given orienation
                        angleOffset = 0
                        for adx in range(len(angles)):
                            angleOffset += abs(angleDifference(threewayAngles[adx],angles[adx]))
                        angleOffset = angleOffset
                        
                        # check to see if it is the new smallest
                        if angleOffset < minAngleOffset:
                            # assign those threewayAngles to minThreewayAngles with appropriate mod
                            minThreewayAngles = [threewayAngles[0], threewayAngles[1], threewayAngles[2]] 
                            minOrientation = orientation
                            minAngleOffset = angleOffset
                            self.verboseVertexDict[vtex] = connected

                        # print('vtex',vtex,'res',res,'orientation',angles[odx],'angles',angles,'threewayAngles',threewayAngles,'minThreewayAngles',minThreewayAngles,'angleOffset',angleOffset)
                
                # assign the minThreewayAngles to capacitorAngles
                self.capacitorAngles[vtex] = minThreewayAngles
            
        for vtex in range(len(vertexPositions)):
            # now find the input and output angles for all resonators            
            connected = self.vertexDict[vtex]
            for rdx in range(len(connected)):
                
                res = connected[rdx]
                capacitorAngles = self.capacitorAngles[vtex]
                angles = self.vertexAngles[vtex]

                resPos = self.resonators[res]
                startPos = [numpy.round(resPos[0],roundDepth), numpy.round(resPos[1],roundDepth)]
                stopPos = [numpy.round(resPos[2],roundDepth), numpy.round(resPos[3],roundDepth)]



                if len(angles)==1:
                    res = connected[0]
                    if res not in self.resonatorInputOuputAngles.keys():
                        self.resonatorInputOuputAngles[res] = [0.0, 0.0]
                else:
                    # find difference between vertex angle and the capacitor angle
                    angleDiff = angleDifference(capacitorAngles[rdx],angles[rdx])
                    # print('vtex',vtex,'rdx',rdx,'capacitorAngles',capacitorAngles,'angles',angles,'angleDiff',angleDiff,'resonatorInputOuputAngles',self.resonatorInputOuputAngles)

                    # assign angleDiff to the input and output
                    if startPos[0] == self.vertexPositions[vtex][0] and startPos[1] == self.vertexPositions[vtex][1]:
                        # print('start','angleDiff',angleDiff,'resonatorInputOuputAngles',self.resonatorInputOuputAngles)
                        if res in self.resonatorInputOuputAngles.keys():
                            # print('3')
                            resonatorInputOuput = self.resonatorInputOuputAngles[res]
                            self.resonatorInputOuputAngles[res] = [-angleDiff,resonatorInputOuput[1]]
                        else:
                            # print('4')
                            self.resonatorInputOuputAngles[res] = [-angleDiff,0.0]
                        # print('start','angleDiff',angleDiff,'resonatorInputOuputAngles',self.resonatorInputOuputAngles)
                    elif stopPos[0] == self.vertexPositions[vtex][0] and stopPos[1] == self.vertexPositions[vtex][1]:
                        # print('stop','angleDiff',angleDiff,'resonatorInputOuputAngles',self.resonatorInputOuputAngles)
                        if res in self.resonatorInputOuputAngles.keys():
                            # print('1')
                            resonatorInputOuput = self.resonatorInputOuputAngles[res]
                            self.resonatorInputOuputAngles[res] = [resonatorInputOuput[0],-angleDiff]
                        else:
                            # print('2')
                            self.resonatorInputOuputAngles[res] = [0.0,-angleDiff]
            #             print('stop','angleDiff',angleDiff,'resonatorInputOuputAngles',self.resonatorInputOuputAngles)
            # print('vtex',vtex,'resonatorInputOuputAngles',self.resonatorInputOuputAngles)

class placeThreewayCap:
    """Draw a section of lattice with well-defined radius"""
    def __init__(self,s,connectedResonators,centerPosition,capacitorAngles,cap_length,cap_gap,gapw_buffer,stub_length,pinw,stop_pinw,gapw,stop_gapw,pinwLinear,gapwLinear,excludedResonatorIndices,bufferDistance):

        # if pinw is None: pinw=structure.defaults['pinw']
        # if gapw is None: gapw=structure.defaults['gapw']
        # if stub_length is None: stub_length = 0.0

        # self.capEndPoints = []
        # self.capEndDirections = []

        for cdx in range(len(connectedResonators)):
            # print('capacitorAngles[cdx]',capacitorAngles[cdx])
            capAngle = capacitorAngles[cdx]

            s.last_direction = capAngle+180      

            if connectedResonators[cdx] == -1 and (cdx not in excludedResonatorIndices):
                s.last = (centerPosition[0]+numpy.cos(capAngle*numpy.pi/180)*(cap_length+stub_length),centerPosition[1]+numpy.sin(capAngle*numpy.pi/180)*(cap_length+stub_length))
                CPWStraight(s,stub_length,pinwLinear,gapwLinear)
                Inner_end_cap_bondpad_buffer(s,cap_length,cap_gap,pinw,stop_pinw,gapw,stop_gapw,gapw_buffer,pinwLinear,gapwLinear,bufferDistance)

            elif cdx in excludedResonatorIndices:
                s.last = (centerPosition[0]+numpy.cos(capAngle*numpy.pi/180)*cap_length,centerPosition[1]+numpy.sin(capAngle*numpy.pi/180)*cap_length)
                Inner_end_cap_bondpad_buffer(s,cap_length,cap_gap,pinw,stop_pinw,gapw,stop_gapw,gapw_buffer,pinwLinear,gapwLinear,bufferDistance)

            else:
                s.last = (centerPosition[0]+numpy.cos(capAngle*numpy.pi/180)*cap_length,centerPosition[1]+numpy.sin(capAngle*numpy.pi/180)*cap_length)
                Inner_end_cap_buffer(s,cap_length,cap_gap,pinw,stop_pinw,gapw,stop_gapw,gapw_buffer,bufferDistance)
                

        
def drawThreewayCaps(structure,drawDictionaries,cap_length,cap_gap,gapw_buffer,stub_length,pinw,stop_pinw,gapw,stop_gapw,bondPadConnectionParameters,bufferDistance):
    

    excludedVertices = bondPadConnectionParameters['startVertices']
    allExcludedResonatorsIndices = bondPadConnectionParameters['resonatorIndices']
    pinwLinear = bondPadConnectionParameters['pinwLinear']
    gapwLinear = bondPadConnectionParameters['gapwLinear']
    
    for vtex in range(len(drawDictionaries.vertexPositions)):
        
        excludedResonatorIndices = []
        if vtex in excludedVertices:

            for idx in range(len(excludedVertices)):
                if excludedVertices[idx]==vtex:
                    excludedResonatorIndices.append(allExcludedResonatorsIndices[idx])
            # print('vtex',vtex,'excludedResonatorIndices',excludedResonatorIndices,'excludedVertices',excludedVertices,'allExcludedResonatorsIndices',allExcludedResonatorsIndices)

        placeThreewayCap(structure,drawDictionaries.verboseVertexDict[vtex],drawDictionaries.vertexPositions[vtex],drawDictionaries.capacitorAngles[vtex],cap_length,cap_gap,gapw_buffer,stub_length,pinw,stop_pinw,gapw,stop_gapw,pinwLinear,gapwLinear,excludedResonatorIndices,bufferDistance)
       

def drawResonators(structure,drawDictionaries,reslength,x,r,cap_length,pinw,gapw,num_wiggles,up,roundDepth=3):
    
    for res in range(len(drawDictionaries.resonators)):
        
        resonator = drawDictionaries.resonators[res]
        plusEnd = numpy.round([resonator[0],resonator[1]],roundDepth)
        minusEnd = numpy.round([resonator[2],resonator[3]],roundDepth)
        startVertex = drawDictionaries.resonatorStartStopVertices[res][0]
        stopVertex = drawDictionaries.resonatorStartStopVertices[res][1]

        whichResonatorsStart = drawDictionaries.verboseVertexDict[startVertex]
        startCapAngles = drawDictionaries.capacitorAngles[startVertex]
        for idx in range(len(whichResonatorsStart)):
            if res == whichResonatorsStart[idx]:
                whichIndex = idx
        startCapAngle = startCapAngles[whichIndex]

        whichResonatorsStop = drawDictionaries.verboseVertexDict[stopVertex]
        stopCapAngles = drawDictionaries.capacitorAngles[stopVertex]
        for idx in range(len(whichResonatorsStop)):
            if res == whichResonatorsStop[idx]:
                whichIndex = idx
        stopCapAngle = stopCapAngles[whichIndex]

        startPosition = [plusEnd[0]+cap_length*cos(startCapAngle*pi/180),plusEnd[1]+cap_length*sin(startCapAngle*pi/180)]
        stopPosition = [minusEnd[0]+cap_length*cos(stopCapAngle*pi/180),minusEnd[1]+cap_length*sin(stopCapAngle*pi/180)]

        x = numpy.sqrt((plusEnd[0]-minusEnd[0])**2 + (plusEnd[1]-minusEnd[1])**2)

        # print('res',res,'startCapAngle',startCapAngle,'stopCapAngle',stopCapAngle,'resonatorInputOuputAngles[res]',drawDictionaries.resonatorInputOuputAngles[res])
        connectionResonator(structure,reslength,x,r,res,resonator,startCapAngle,stopCapAngle,startPosition,drawDictionaries.resonatorInputOuputAngles[res],cap_length,pinw,gapw,num_wiggles,up)


class connectionResonator:
    """Draw a section of lattice with well-defined radius"""
    def __init__(self,s,total_length,x,r,res,resonator,startCapAngle,stopCapAngle,startPosition,resonatorInputOutputAngles,cap_length,pinw=None,gapw=None,num_wiggles=None,up=None):
        
        if total_length == 0:
            self.vlength=0
            return
        
        if pinw is None: pinw=structure.defaults['pinw']
        if gapw is None: gapw=structure.defaults['gapw']
        if num_wiggles is None: num_wiggles=1
        if up is None: up = 1

        # move to start position and set proper orientation
        s.last = startPosition
        s.last_direction=startCapAngle       

        # get the input and output angles
        input_angle = resonatorInputOutputAngles[0]
        output_angle = resonatorInputOutputAngles[1]
        input_angle_radians = input_angle*numpy.pi/180.
        output_angle_radians = output_angle*numpy.pi/180.

        # set the proper signs
        if input_angle>0:
            inSign = 1
        elif input_angle<0:
            inSign = -1
        else:
            inSign = 0

        if output_angle>0:
            outSign = 1
        elif output_angle<0:
            outSign = -1
        else:
            outSign = 0

        # compute all the lengths for the input section
        inputCapDeltaX = abs((cap_length) * cos(input_angle_radians))
        inputTurn1DeltaX = abs(2*r*sin(input_angle_radians))
        inputTurn2DeltaX = abs(r*sin(input_angle_radians))
        inputTurn2DeltaY = abs(r*(1-cos(input_angle_radians)))

        # compute how much further you need to go to meet the stright line resonator path
        inputYNeeded = abs(cap_length*sin(input_angle_radians))
        inputExtensionY = abs(inputTurn2DeltaY- inputYNeeded)

        # figure out what the extension length needs to be
        if numpy.round(abs(input_angle_radians),3)!=0:
            inputExtension = inputExtensionY/abs(sin(input_angle_radians))
        else:
            inputExtension = 0.0
        inputExtensionDeltaX = inputExtension*cos(input_angle_radians)
        
        # compute the total X distance and the total distance
        inputAngleCorrectionX = inputCapDeltaX + inputTurn1DeltaX + inputExtensionDeltaX + inputTurn2DeltaX
        inputAngleCorrectionLength = cap_length + r*2*abs(input_angle_radians) + inputExtension + r*abs(input_angle_radians)

        # compute all the lengths for the output section
        outputCapDeltaX = abs((cap_length) * cos(output_angle_radians))
        outputTurn2DeltaX = abs(2*r*sin(output_angle_radians))
        outputTurn1DeltaX = abs(r*sin(output_angle_radians))
        outputTurn1DeltaY = abs(r*(1-cos(output_angle_radians)))

        # compute how much further you need to go to meet the stright line resonator path
        outputYNeeded = abs(cap_length*sin(output_angle_radians))
        outputExtensionY = abs(outputTurn1DeltaY-outputYNeeded)

        # figure out what the extension length needs to be
        if numpy.round(abs(output_angle_radians),3)!=0:
            outputExtension = outputExtensionY/abs(sin(output_angle_radians))
        else:
            outputExtension = 0.0
        outputExtensionDeltaX = outputExtension*cos(output_angle_radians)

        # compute the total X distance and the total distance
        outputAngleCorrectionX = outputCapDeltaX + outputTurn1DeltaX + outputExtensionDeltaX + outputTurn2DeltaX
        outputAngleCorrectionLength = cap_length + r*2*abs(output_angle_radians) + outputExtension + r*abs(output_angle_radians)

        # print('res',res,'x',x,'inputAngleCorrectionX',inputAngleCorrectionX,'outputAngleCorrectionX',outputAngleCorrectionX)

        # compute how long the extension length along the straight line needs to be to make it to the other capacitor
        extension_length = (x-inputAngleCorrectionX-outputAngleCorrectionX-2*r-2*r*num_wiggles)/2

        if num_wiggles % 2 ==1:
            vlength = (((total_length-inputAngleCorrectionLength-outputAngleCorrectionLength-2*extension_length-pi*r)/num_wiggles)-pi*r)/2
        else:
            vlength = (((total_length-inputAngleCorrectionLength-outputAngleCorrectionLength-2*extension_length-pi*r+2*r)/num_wiggles)-pi*r-2*r)/2
        
        talliedLength = cap_length
        # print('added', cap_length-cap_gap1/2)
        # print('talliedLength',talliedLength)

        

        # draw the input section
        CPWBend(s,2*input_angle,pinw,gapw,r)
        CPWStraight(s,inputExtension,pinw,gapw)
        CPWBend(s,-input_angle,pinw,gapw,r)
        talliedLength = talliedLength + abs(3*r*input_angle_radians)+inputExtension
        # print('added', abs(r*input_angle))
        # print('talliedLength',talliedLength)

        # extend
        CPWStraight(s,extension_length,pinw,gapw)
        talliedLength = talliedLength + extension_length

        # draw squiggles
        if num_wiggles == 1:
            CPWBend(s,up*90,pinw,gapw,r)
            talliedLength = talliedLength + r*pi/2
            # print('added',r*pi/2)
            # print('talliedLength',talliedLength)    
              
            CPWStraight(s,vlength,pinw,gapw)
            talliedLength = talliedLength + vlength
            # print('added',vlength)
            # print('talliedLength',talliedLength)    

            CPWBend(s,-up*180,pinw,gapw,r)
            talliedLength = talliedLength + r*pi
            # print('added',r*pi)
            # print('talliedLength',talliedLength)     

            CPWStraight(s,vlength,pinw,gapw)
            talliedLength = talliedLength + vlength

            CPWBend(s,up*90,pinw,gapw,r)
            talliedLength = talliedLength + r*pi/2

        else:        
            CPWBend(s,-90,pinw,gapw,r)
            talliedLength = talliedLength + r*pi/2
                # print('added',r*pi/2)
                # print('talliedLength',talliedLength)    

            for ii in range(num_wiggles):
                isign=(2*(ii%2)-1)                
                CPWStraight(s,vlength,pinw,gapw)
                talliedLength = talliedLength + vlength
                # print('added',vlength)
                # print('talliedLength',talliedLength)    

                CPWBend(s,-isign*180,pinw,gapw,r)
                talliedLength = talliedLength + r*pi
                # print('added',r*pi)
                # print('talliedLength',talliedLength)     

                CPWStraight(s,vlength,pinw,gapw)
                talliedLength = talliedLength + vlength

                if num_wiggles % 2 == 0 and ii != (num_wiggles-1):
                    CPWStraight(s,2*r,pinw,gapw)
                    talliedLength = talliedLength + 2*r
                # print('added',vlength)
                # print('talliedLength',talliedLength)     

            CPWBend(s,isign*90,pinw,gapw,r)
            talliedLength = talliedLength + r*pi/2
            # print('added',r*pi/2)
            # print('talliedLength',talliedLength)  

        # extend
        CPWStraight(s,extension_length,pinw,gapw)
        talliedLength = talliedLength + extension_length
        # print('added',extension_length)
        # print('talliedLength',talliedLength)  

        CPWBend(s,output_angle,pinw,gapw,r)
        CPWStraight(s,outputExtension,pinw,gapw)
        CPWBend(s,-2*output_angle,pinw,gapw,r)

        # CPWBend(s,output_angle,pinw,gapw,r)
        talliedLength = talliedLength + abs(3*r*output_angle_radians)+outputExtension
        # print('added',abs(r*output_angle))
        # print('talliedLength',talliedLength)     

        talliedLength = talliedLength + cap_length
        # print('added',cap_length)
        # print('res',res,'talliedLength',talliedLength,'abs(3*r*output_angle)+outputExtension', abs(3*r*output_angle)+outputExtension)
        print('res',res,'talliedLength',talliedLength)
           

class bondPadConnection:
    """Makes a connection between the edge of a ring and a bond pad"""
    def __init__(self,structure,drawDictionaries,extension_length,lowerLeftChipCorner,chipSize,edgeBufferDistance,directLine,r,bondPadPos,bondPadAng,startVertex,startRes,cap_length,pinw,gapw):

        s = structure
        verboseVertex = drawDictionaries.verboseVertexDict[startVertex]
        capacitorAngles = drawDictionaries.capacitorAngles[startVertex]
        vertexPosition = drawDictionaries.vertexPositions[startVertex]


        resonatorAngle = capacitorAngles[startRes]



        xMax = lowerLeftChipCorner[0] + chipSize - edgeBufferDistance
        xMin = lowerLeftChipCorner[0] + edgeBufferDistance

        yMax = lowerLeftChipCorner[1] + chipSize - edgeBufferDistance
        yMin = lowerLeftChipCorner[1] + edgeBufferDistance

        s.last = [vertexPosition[0]+(cap_length)*numpy.cos(resonatorAngle*numpy.pi/180.),vertexPosition[1]+(cap_length)*numpy.sin(resonatorAngle*numpy.pi/180.)]
        s.last_direction = resonatorAngle
        CPWStraight(s,extension_length,pinw,gapw)
        
        s.last_direction = s.last_direction % 360


        if 0<=s.last_direction<=90 or 270<s.last_direction<=360:
            rotationAngle1 = angleDifference(0,s.last_direction)
            CPWBend(s,rotationAngle1,pinw,gapw,r)

            if directLine==1:
                xDist1 = bondPadPos[0]-s.last[0]
                yDist1 = bondPadPos[1]-s.last[1]
                if xDist1-r>0:
                    if yDist1 > 0:
                        CPWStraight(s,xDist1-r,pinw,gapw)
                        CPWBend(s,90,pinw,gapw,r)
                        CPWStraight(s,yDist1-r,pinw,gapw)
                    else:
                        CPWStraight(s,xDist1-r,pinw,gapw)
                        CPWBend(s,-90,pinw,gapw,r)
                        CPWStraight(s,abs(yDist1)-r,pinw,gapw)
                else:
                    if yDist1 > 0:
                        CPWStraight(s,xDist1+2*r,pinw,gapw)
                        CPWBend(s,90,pinw,gapw,r)
                        CPWStraight(s,yDist1-3*r,pinw,gapw)
                        CPWBend(s,90,pinw,gapw,r)
                        xDist2 = s.last[0]-bondPadPos[0]-r
                        CPWStraight(s,xDist2,pinw,gapw)
                        CPWBend(s,-90,pinw,gapw,r)
                    else:
                        CPWStraight(s,xDist1+2*r,pinw,gapw)
                        CPWBend(s,-90,pinw,gapw,r)
                        CPWStraight(s,abs(yDist1)-3*r,pinw,gapw)
                        CPWBend(s,-90,pinw,gapw,r)
                        xDist2 = s.last[0]-bondPadPos[0]-r
                        CPWStraight(s,xDist2,pinw,gapw)
                        CPWBend(s,90,pinw,gapw,r)

            else:
                xDist1 = xMax-s.last[0]-r
                if s.last[1] > bondPadPos[1]:
                    yDist1 = s.last[1]-3*r-bondPadPos[1]
                    CPWStraight(s,xDist1,pinw,gapw)
                    CPWBend(s,-90,pinw,gapw,r)
                    CPWStraight(s,yDist1,pinw,gapw)
                    CPWBend(s,-90,pinw,gapw,r)
                else:
                    yDist1 = bondPadPos[1]-3*r-s.last[1]
                    CPWStraight(s,xDist1,pinw,gapw)
                    CPWBend(s,90,pinw,gapw,r)
                    CPWStraight(s,yDist1,pinw,gapw)
                    CPWBend(s,90,pinw,gapw,r)


                xDist2 = s.last[0]-bondPadPos[0]-r
                CPWStraight(s,xDist2,pinw,gapw)
                if s.last[1]>bondPadPos[1]:
                    CPWBend(s,90,pinw,gapw,r)
                else:
                    CPWBend(s,-90,pinw,gapw,r)

        else:
            rotationAngle1 = angleDifference(180,s.last_direction)
            CPWBend(s,rotationAngle1,pinw,gapw,r)
            
            if directLine==1:
                xDist1 = s.last[0]-bondPadPos[0]
                yDist1 = s.last[1]-bondPadPos[1]

                if xDist1-r>0:
                    if yDist1 > 0:
                        CPWStraight(s,xDist1-r,pinw,gapw)
                        CPWBend(s,90,pinw,gapw,r)
                        CPWStraight(s,yDist1-r,pinw,gapw)
                    else:
                        CPWStraight(s,xDist1-r,pinw,gapw)
                        CPWBend(s,-90,pinw,gapw,r)
                        CPWStraight(s,abs(yDist1)-r,pinw,gapw)
                else:
                    if yDist1 > 0:
                        CPWStraight(s,xDist1+2*r,pinw,gapw)
                        CPWBend(s,90,pinw,gapw,r)
                        CPWStraight(s,yDist1-3*r,pinw,gapw)
                        CPWBend(s,90,pinw,gapw,r)
                        xDist2 = bondPadPos[0]-s.last[0]-r
                        CPWStraight(s,xDist2,pinw,gapw)
                        CPWBend(s,-90,pinw,gapw,r)
                    else:
                        CPWStraight(s,xDist1+2*r,pinw,gapw)
                        CPWBend(s,-90,pinw,gapw,r)
                        CPWStraight(s,abs(yDist1)-3*r,pinw,gapw)
                        CPWBend(s,-90,pinw,gapw,r)
                        xDist2 = bondPadPos[0]-s.last[0]-r
                        CPWStraight(s,xDist2,pinw,gapw)
                        CPWBend(s,90,pinw,gapw,r)

            else:
                xDist1 = s.last[0]-xMin-r

                if s.last[1] > bondPadPos[1]:
                    yDist1 = s.last[1]-3*r-bondPadPos[1]
                    CPWStraight(s,xDist1,pinw,gapw)
                    CPWBend(s,90,pinw,gapw,r)
                    CPWStraight(s,yDist1,pinw,gapw)
                    CPWBend(s,90,pinw,gapw,r)

                else:
                    yDist1 = bondPadPos[1]-3*r-s.last[1]
                    CPWStraight(s,xDist1,pinw,gapw)
                    CPWBend(s,-90,pinw,gapw,r)
                    CPWStraight(s,yDist1,pinw,gapw)
                    CPWBend(s,-90,pinw,gapw,r)


                xDist2 = bondPadPos[0]-s.last[0]-r
                CPWStraight(s,xDist2,pinw,gapw)

                if s.last[1]>bondPadPos[1]:
                    CPWBend(s,-90,pinw,gapw,r)
                else:
                    CPWBend(s,90,pinw,gapw,r)



        

class bondPadConnectionHyperbolic:
    """Makes a connection between the edge of a ring and a bond pad"""
    def __init__(self,structure,radialExtension,r,startPos,startAngle,stopPos,stopAngle,cap_length,cap_gap,cap_gap_original,pinw,stop_pinw,gapw,stop_gapw):

        s = structure
#        cap_gap = cap_gap_original
        cap = cap_length + cap_gap_original/2
        
        initX = (cap+radialExtension)*cos(startAngle*pi/180)
        initY = (cap+radialExtension)*sin(startAngle*pi/180)
        
        deltaX = stopPos[0] - startPos[0]
        deltaY = stopPos[1] - startPos[1]
        
        if deltaX > 0 and deltaY > 0:
            
            initialRotationAngle = (stopAngle-90)-startAngle
            finalRotationAngle = 90
        
        elif deltaX > 0 and deltaY < 0:
            
            initialRotationAngle = ((stopAngle+90)-startAngle)%360
            finalRotationAngle = -90
            
        elif deltaX < 0 and deltaY < 0:
        
            initialRotationAngle = -((startAngle+(stopAngle-90))%360)
            finalRotationAngle = 90
        
        else:
        
            initialRotationAngle = ((stopAngle+90)-startAngle)%360
            finalRotationAngle = -90
            
        
        initialRotationDeltaX = r*(abs(sin(initialRotationAngle*pi/180)))
        initialRotationDeltaY = r*(1-abs(cos(initialRotationAngle*pi/180)))
        finalRotationDeltaX = r*(1-cos(finalRotationAngle*pi/180))
        finalRotationDeltaY = r*abs(sin(finalRotationAngle*pi/180))
        
        if deltaX < 0:    
            xExcursion = -deltaX + initX - initialRotationDeltaX - finalRotationDeltaX
        else:
            xExcursion = deltaX - initX - initialRotationDeltaX - finalRotationDeltaX
        
        if deltaY < 0:
            yExcursion = -deltaY + initY - initialRotationDeltaY - finalRotationDeltaY
        else:
            yExcursion = deltaY - initY - initialRotationDeltaY - finalRotationDeltaY
            
            
        # print('deltaX',deltaX,'deltaY',deltaY,'initialRotationAngle',initialRotationAngle,'finalRotationAngle',finalRotationAngle,
        #       'initialRotationDeltaX',initialRotationDeltaX,'initialRotationDeltaY',initialRotationDeltaY,
        #       'finalRotationDeltaX',finalRotationDeltaX,'finalRotationDeltaY',finalRotationDeltaY,'startAngle',startAngle,
        #       'stopAngle',stopAngle)
        
        s.last_direction = startAngle
        s.last = startPos
        
        s.last = (startPos[0] + cap * cos(s.last_direction*pi/180), startPos[1] + cap * sin(s.last_direction*pi/180)) 
        s.last_direction = s.last_direction - 180;
#        Inner_end_cap(s,cap_length,cap_gap_original,pinw,stop_pinw,gapw,stop_gapw)  
        outerHalfPac(s,cap_length,cap_gap,cap_gap_original,pinw,stop_pinw,gapw,stop_gapw)   
        s.last_direction = s.last_direction + 180;
        s.last = (startPos[0] + cap * cos(s.last_direction*pi/180), startPos[1] + cap * sin(s.last_direction*pi/180))  
        
        
        CPWStraight(s,radialExtension,pinw,gapw)
        CPWBend(s,initialRotationAngle,pinw,gapw,r)
        CPWStraight(s,xExcursion,pinw,gapw)
        CPWBend(s,finalRotationAngle,pinw,gapw,r)
        CPWStraight(s,yExcursion,pinw,gapw)
        
        
class outerHalfPac:
     def __init__(self,structure,cap_length,cap_gap,cap_gap_original,start_pinw,stop_pinw,start_gapw,stop_gapw):
        if cap_length==0: return 
        if cap_length==0: return 

        s=structure
        start=s.last
        start_taperX= cap_length-(stop_pinw/2)/tan(60*pi/180) # X-pos of where the centerpin taper starts
        start_taperY=((cap_length-start_taperX)+cap_gap/2)*tan(60*pi/180)
    
        "Intersection point calculations"
        x_intersect=(start_pinw/2 + start_gapw - sqrt(3)*(cap_length+cap_gap_original/2))/(-sqrt(3)-(stop_gapw/start_taperX))        
        y_intersect=-sqrt(3)*(x_intersect-(cap_length+cap_gap_original/2))
        "draw points that form the end cap geometry"
        EndCap=[
            (start[0],start[1]+start_pinw/2),
            (start[0],start[1]+start_pinw/2+start_gapw),
            (start[0]+x_intersect,start[1]+y_intersect),
            (start[0]+cap_length+cap_gap_original/2,start[1]),
            (start[0]+x_intersect,start[1]-y_intersect),
            (start[0],start[1]-start_pinw/2-start_gapw),
            (start[0],start[1]-start_pinw/2),
            (start[0]+start_taperX,start[1]-stop_pinw/2),
            (start[0]+cap_length,start[1]),
            # (start[0]+start_taperX,start[1]+stop_pinw/2),
            (start[0]+(cap_length+start_taperX)/2,start[1]+(stop_pinw/2)/2),
            (start[0],start[1]+start_pinw/2)
            ]

        "rotate structure to proper orientation"
        EndCap=rotate_pts(EndCap,s.last_direction,start)

        "create polylines and append to drawing /connect the dots"
        s.append(sdxf.PolyLine(EndCap))
            
        "update last anchor position"
        stop=rotate_pt((start[0]+cap_length+cap_gap/2,start[1]),s.last_direction,start)
        s.last=stop
        
class drawSingleAlignmentMark:
    """ Make a single alignment mark"""
    def __init__(self,structure,linewidth,size,startPt):
        lw=linewidth/2.
        w=size[0]/2.
        h=size[1]/2.
        pts=[ (-lw,-h), (lw,-h), (lw,-lw),(w,-lw),(w,lw),(lw,lw),(lw,h),(-lw,h),(-lw,lw),(-w,lw),(-w,-lw),(-lw,-lw),(-lw,-h)]
        
        pts=translate_pts(pts,startPt)
        
        structure.append(sdxf.PolyLine(pts))
    
class drawAlignmentMarks:
    """ Draw a set of alignment marks"""
    def __init__(self,structure,linewidth,size,edgeBuffer,chipSize,pos):
        drawSingleAlignmentMark(structure,linewidth,size,[pos[0]+edgeBuffer,pos[1]-edgeBuffer])
        drawSingleAlignmentMark(structure,linewidth,size,[pos[0]+chipSize-edgeBuffer,pos[1]-edgeBuffer])
        drawSingleAlignmentMark(structure,linewidth,size,[pos[0]+chipSize-edgeBuffer,pos[1]-chipSize+edgeBuffer])
        drawSingleAlignmentMark(structure,linewidth,size,[pos[0]+edgeBuffer,pos[1]-chipSize+edgeBuffer])

class drawCycloneResonator:
    """Draw a section of lattice with well-defined radius"""
    def __init__(self,structure,total_length,x,radius,pinw,gapw,num_wiggles=None):
        """ """
        
        if total_length == 0:
            self.vlength=0
            return

        s=structure
        start=structure.last    

        if num_wiggles is None: num_wiggles=3     
        
        extension_length = (x - radius*(2+2*num_wiggles))/2
       
        vlength = (((total_length - 2*extension_length - 2*pi/2*radius)/num_wiggles)-pi*radius)/2

        s.last = (start[0], start[1])
        tempLast = s.last_direction
        s.last_direction = s.last_direction - 180;        

        #     # Outer_Pacman_cap(s,cap_length,cap_gap1,pinw,stop_pinw,gapw,stop_gapw)
        # drawCycloneCap(s,cap_length,cap_separation,bendradius1,bendradius2,pinw,gapw,indentGapw) 

        # s.last_direction = tempLast;
        # s.last = (start[0] + cap_length * cos(s.last_direction*pi/180), start[1] + cap_length * sin(s.last_direction*pi/180))  

        # talliedLength = cap_length
        # print('added', abs(radius*input_angle))
        # print('talliedLength',talliedLength)

        CPWStraight(s,extension_length,pinw,gapw)
        talliedLength = extension_length
        # print('added', extension_length)
        # print('talliedLength',talliedLength)           
        if num_wiggles == 1:
            CPWBend(s,90,pinw,gapw,radius)  
            talliedLength = talliedLength + radius*pi/2
            # print('added',radius*pi/2)
            # print('talliedLength',talliedLength)  

            CPWStraight(s,vlength,pinw,gapw)
            talliedLength = talliedLength + vlength
            # print('added',vlength)
            # print('talliedLength',talliedLength) 

            CPWBend(s,-180,pinw,gapw,radius)
            talliedLength = talliedLength + pi*radius
            # print('added',pi*radius)
            # print('talliedLength',talliedLength) 
            
            CPWStraight(s,vlength,pinw,gapw)
            talliedLength = talliedLength + vlength
            # print('added',vlength)
            # print('talliedLength',talliedLength) 
            
            CPWBend(s,90,pinw,gapw,radius)
            talliedLength = talliedLength + radius*pi/2
            # print('added',radius*pi/2)
            # print('talliedLength',talliedLength) 

            CPWStraight(s,extension_length,pinw,gapw)
            talliedLength = talliedLength + extension_length
            # print('added',extension_length)
            # print('talliedLength',talliedLength)

            # print('added',abs(output_angle*radius))
            # print('talliedLength',talliedLength) 
        else:
            CPWBend(s,90,pinw,gapw,radius)
            talliedLength = talliedLength + radius*pi/2
            # print('added',radius*pi/2)
            # print('talliedLength',talliedLength)    

            for ii in range(num_wiggles):
                isign=2*(ii%2)-1                
                CPWStraight(s,vlength,pinw,gapw)
                talliedLength = talliedLength + vlength
                # print('added',vlength)
                # print('talliedLength',talliedLength)    

                CPWBend(s,isign*180,pinw,gapw,radius)
                talliedLength = talliedLength + radius*pi
                # print('added',radius*pi)
                # print('talliedLength',talliedLength)     

                CPWStraight(s,vlength,pinw,gapw)
                talliedLength = talliedLength + vlength
                # print('added',vlength)
                # print('talliedLength',talliedLength)     

            CPWBend(s,90,pinw,gapw,radius)
            talliedLength = talliedLength + radius*pi/2
            # print('added',radius*pi/2)
            # print('talliedLength',talliedLength)  

            CPWStraight(s,extension_length,pinw,gapw)
            talliedLength = talliedLength + extension_length
            # print('added',extension_length)
            # print('talliedLength',talliedLength)   
   
class cycloneInterior:
    """A CPW bend"""
    def __init__(self,structure,turn_angle,pinw=None,gapw=None,radius=None,segments=60):
        """creates a CPW bend with pinw/gapw/radius
            @param turn_angle: turn_angle is in degrees, positive is CCW, negative is CW
        """
        #load default values if necessary
        
        if turn_angle==0: return
        
        s=structure
#        print('radius',radius)

        if radius is None: radius=s.defaults['radius']
        if pinw is None:   pinw=s.defaults['pinw']
        if gapw is None:   gapw=s.defaults['gapw']
        
        self.structure=structure
        self.turn_angle=turn_angle
        self.pinw=pinw
        self.gapw=gapw
        self.radius=radius
        self.segments=segments

        self.start=s.last
        self.start_angle=s.last_direction
        self.stop_angle=self.start_angle+self.turn_angle
        
        if turn_angle>0: self.asign=1
        else:            self.asign=-1
       
        #DXF uses the angle of the radial vector for its start and stop angles
        #so we have to rotate our angles by 90 degrees to get them right
        #also it only knows about arcs with CCW sense to them, so we have to rotate our angles appropriately
        self.astart_angle=self.start_angle-self.asign*90
        self.astop_angle=self.stop_angle-self.asign*90
        #calculate location of Arc center
        self.center=rotate_pt( (self.start[0],self.start[1]+self.asign*self.radius),self.start_angle,self.start)
        
        self.poly_arc_bend(structure)
        
        self.structure.last=rotate_pt(self.start,self.stop_angle-self.start_angle,self.center)
        self.structure.last_direction=self.stop_angle

    def poly_arc_bend(self,structure):
    
        # #lower gap
        # pts=arc_pts(self.astart_angle,self.astop_angle,self.radius+self.pinw/2.+self.gapw,self.segments)
        # pts.extend(arc_pts(self.astop_angle,self.astart_angle,self.radius+self.pinw/2.,self.segments))
        # pts.append(pts[0])

        pts=arc_pts(self.astart_angle,self.astop_angle,self.radius-self.pinw/2,self.segments)
        pts.extend(arc_pts(self.astop_angle,self.astart_angle,self.radius+self.pinw/2,self.segments))
        pts.append(pts[0])
        self.structure.append(sdxf.PolyLine(translate_pts(pts,self.center),layer=structure.layer,color=structure.color))
        # self.structure.append(sdxf.PolyLine(translate_pts(pts2,self.center)))
   
class cycloneBoundingCircle:
    """A CPW bend"""
    def __init__(self,structure,center,radius,segments=200):

        
        s=structure
    
        self.structure=structure
        self.radius=radius
        self.segments=segments

        self.start=center
        self.start_angle=0
        self.stop_angle=360
        
       
        #DXF uses the angle of the radial vector for its start and stop angles
        #so we have to rotate our angles by 90 degrees to get them right
        #also it only knows about arcs with CCW sense to them, so we have to rotate our angles appropriately
        self.astart_angle=self.start_angle
        self.astop_angle=self.stop_angle
        #calculate location of Arc center
        self.center=rotate_pt( (self.start[0],self.start[1]+self.radius),self.start_angle,self.start)
        
        self.poly_arc_bend(structure)
        
        self.structure.last=rotate_pt(self.start,self.stop_angle-self.start_angle,self.center)
        self.structure.last_direction=self.stop_angle

    def poly_arc_bend(self,structure):
    
        # #lower gap
        # pts=arc_pts(self.astart_angle,self.astop_angle,self.radius+self.pinw/2.+self.gapw,self.segments)
        # pts.extend(arc_pts(self.astop_angle,self.astart_angle,self.radius+self.pinw/2.,self.segments))
        # pts.append(pts[0])

        pts=arc_pts(self.astart_angle,self.astop_angle,self.radius,self.segments)
        pts.extend(arc_pts(self.astop_angle,self.astart_angle,self.radius,self.segments))
        pts.append(pts[0])
        self.structure.append(sdxf.PolyLine(translate_pts(pts,self.center),layer=structure.layer,color = structure.color))

class drawCycloneCap:
     def __init__(self,interiorStructure,boundingCircleStructure,cap_length,cap_separation,bendradius,outerRadius,pinw,gapw,indentGapw):
        """
        Class that draws a singlehexagonal endcap for one part of the 3way coupling capacitor
        variables:
        cap_length= linear length of end cap
        cap_gap= width of capacitive gap b/n end caps
        start_pinw= beginning width of end cap
        stop_pinw= width of end cap before taper
        """
        
        """
        The issue with this code is that for small cap_gap, the space between the capacitor and ground plane isn't big enough.
        """
        "Load attributes"    

        s1=interiorStructure
        s2=boundingCircleStructure


        start=s1.last 
        startAngle = s1.last_direction
        
        overlapAngle = 30
        rotationAngle = 45
              
        deltaX = 2*(bendradius+gapw+pinw)*sin(rotationAngle*pi/180)-(bendradius)
        deltaY = 2*(bendradius+gapw+pinw)*(1-cos(rotationAngle*pi/180))
        innerRadius = outerRadius-deltaY
        print('innerRadius',innerRadius)
        beta = atan(bendradius/(2.0*innerRadius + bendradius)) * 180/pi 
        # input_angle = (bendradius/outerRadius)*180/pi


        for i in range(4):
            cycloneInterior(s1,-90,pinw,gapw,bendradius)
            cycloneInterior(s1,90,pinw,gapw,outerRadius)
            referencePosition = s1.last #remember where this sits
            referenceAngle = s1.last_direction
            cycloneInterior(s1,rotationAngle,pinw,gapw,bendradius+gapw+pinw)
            cycloneInterior(s1,-rotationAngle,pinw,gapw,bendradius+gapw+pinw)
            # print('last_direction',s1.last_direction)
            # s1.last_direction = s1.last_direction + 180
            # cycloneInterior(s1,overlapAngle,pinw,gapw,innerRadius)
        
        
            print(deltaY)

            if i == 0:
                s1.last = [referencePosition[0]+deltaX , referencePosition[1]+(deltaY-gapw-pinw-bendradius)]
                s1.last_direction = 90
            elif i==1:
                s1.last = [referencePosition[0]-(deltaY-gapw-pinw-bendradius),referencePosition[1]+deltaX ]
                s1.last_direction = 180
            elif i==2:
                s1.last = [referencePosition[0]-deltaX ,referencePosition[1]-(deltaY-gapw-pinw-bendradius)]
                s1.last_direction = 270




        cycloneBoundingCircle(s2,start,40)


class drawRingResonator:
    """Draw a section of lattice with well-defined radius"""
    def __init__(self,s,total_length,x,r,input_angle,output_angle,cap_length,isInputPac,isOutputPac,cap_gap1,cap_gap2,pinw=None,stop_pinw=None,gapw=None,stop_gapw=None,num_wiggles=None,isOuterRing=0,up=None):
        if total_length == 0:
            self.vlength=0
            return

        start=s.last
        # print('start is ', start)        


        if pinw is None: pinw=structure.defaults['pinw']
        if gapw is None: gapw=structure.defaults['gapw']
        if num_wiggles is None: num_wiggles=3
        if isOuterRing is None: isOuterRing=0
        if up is None: up=1
        
        input_angle = input_angle*pi/180
        output_angle = output_angle*pi/180

        print('input_angle', input_angle, 'output_angle',output_angle)


        if input_angle>0:
            inSign = 1
        elif input_angle<0:
            inSign = -1
        else:
            inSign = 0

        if output_angle>0:
            outSign = 1
        elif output_angle<0:
            outSign = -1
        else:
            outSign = 0

        capDeltaX1 = (cap_length) * cos(input_angle)
        capDeltaY1 = -(cap_length) * sin(input_angle)
        deltaX1 = abs(r*sin(input_angle))
        deltaY1 = -r*(1-cos(input_angle))*inSign

        capDeltaX2 = (cap_length) * cos(output_angle)
        capDeltaY2 = (cap_length) * sin(output_angle)
        deltaX2 = abs(r*sin(output_angle))
        deltaY2 = r*(1-cos(output_angle))*outSign



        # print('input_angle', input_angle, 'output_angle',output_angle,'capDeltaX1',capDeltaX1,'capDeltaY1',capDeltaY1,'deltaX1',deltaX1,'deltaY1',deltaY1,'capDeltaX2',capDeltaX2,'capDeltaY2',capDeltaY2,'deltaX2',deltaX2,'deltaY2',deltaY2)
    
        extension_length = (x-capDeltaX1-deltaX1-2*r-2*r*num_wiggles-capDeltaX2-deltaX2)/2

        if num_wiggles % 2 ==1:
            vlength = (((total_length-2*cap_length-r*(abs(input_angle)+abs(output_angle))-2*extension_length-pi*r)/num_wiggles)-pi*r)/2
        else:
            vlength = (((total_length-2*cap_length-r*(abs(input_angle)+abs(output_angle))-2*extension_length-pi*r-2*r)/num_wiggles)-pi*r-r)/2

        # print('extension_length',extension_length,'vlength',vlength)
        s.last = (start[0] + (cap_length) * cos(s.last_direction*pi/180), start[1] + (cap_length) * sin(s.last_direction*pi/180))
        s.last_direction = s.last_direction - 180;        
        if isInputPac:  
            Outer_Pacman_cap(s,cap_length,cap_gap1,pinw,stop_pinw,gapw,stop_gapw)
        else:
            Inner_end_cap(s,cap_length,cap_gap1,pinw,stop_pinw,gapw,stop_gapw)   
        s.last_direction = s.last_direction+180;
        s.last = (start[0] + (cap_length) * cos(s.last_direction*pi/180), start[1] + (cap_length) * sin(s.last_direction*pi/180))  

        talliedLength = cap_length
        # print('added', cap_length-cap_gap1/2)
        # print('talliedLength',talliedLength)

        CPWBend(s,input_angle*180/pi,pinw,gapw,r)
        talliedLength = talliedLength + abs(r*input_angle)
        # print('added', abs(r*input_angle))
        # print('talliedLength',talliedLength)

        CPWStraight(s,extension_length,pinw,gapw)
        talliedLength = talliedLength + extension_length

        if num_wiggles == 1:
            CPWBend(s,up*90,pinw,gapw,r)
            talliedLength = talliedLength + r*pi/2
            # print('added',r*pi/2)
            # print('talliedLength',talliedLength)    
              
            CPWStraight(s,vlength,pinw,gapw)
            talliedLength = talliedLength + vlength
            # print('added',vlength)
            # print('talliedLength',talliedLength)    

            CPWBend(s,-up*180,pinw,gapw,r)
            talliedLength = talliedLength + r*pi
            # print('added',r*pi)
            # print('talliedLength',talliedLength)     

            CPWStraight(s,vlength,pinw,gapw)
            talliedLength = talliedLength + vlength

            CPWBend(s,up*90,pinw,gapw,r)
            talliedLength = talliedLength + r*pi/2

        else:        
            CPWBend(s,-90,pinw,gapw,r)
            talliedLength = talliedLength + r*pi/2
                # print('added',r*pi/2)
                # print('talliedLength',talliedLength)    

            print('range', range(4), 'num_wiggles', num_wiggles)
            for ii in range(num_wiggles):
                isign=(2*(ii%2)-1)                
                CPWStraight(s,vlength,pinw,gapw)
                talliedLength = talliedLength + vlength
                # print('added',vlength)
                # print('talliedLength',talliedLength)    

                CPWBend(s,-isign*180,pinw,gapw,r)
                talliedLength = talliedLength + r*pi
                # print('added',r*pi)
                # print('talliedLength',talliedLength)     

                CPWStraight(s,vlength,pinw,gapw)
                talliedLength = talliedLength + vlength

                if num_wiggles % 2 == 0 and ii != (num_wiggles-1):
                    CPWStraight(s,2*r,pinw,gapw)
                    talliedLength = talliedLength + 2*r

                # print('added',vlength)
                # print('talliedLength',talliedLength)     

            CPWBend(s,isign*90,pinw,gapw,r)
            talliedLength = talliedLength + r*pi/2



        # print('added',r*pi/2)
        # print('talliedLength',talliedLength)  

        CPWStraight(s,extension_length,pinw,gapw)
        talliedLength = talliedLength + extension_length
        # print('added',extension_length)
        # print('talliedLength',talliedLength)   

        CPWBend(s,output_angle*180/pi,pinw,gapw,r)
        talliedLength = talliedLength + abs(output_angle*r)
        # print('added',abs(r*output_angle))
        # print('talliedLength',talliedLength)     
            
        if isOutputPac:
            Outer_Pacman_cap(s,cap_length,cap_gap2,pinw,stop_pinw,gapw,stop_gapw) 
        else:
            Inner_end_cap(s,cap_length,cap_gap2,pinw,stop_pinw,gapw,stop_gapw) 

        talliedLength = talliedLength + cap_length
        # print('added',cap_length)
        # print('talliedLength',talliedLength)




class drawUnitCell:
    """Draw a section of lattice with well-defined radius"""
    def __init__(self,s,reslength,x,r,cap_length,cap_gap,pinw=None,stop_pinw=None,gapw=None,stop_gapw=None,excludeLast=None):
      
        if excludeLast is None: excludeLast=0
        if pinw is None: pinw=structure.defaults['pinw']
        if gapw is None: gapw=structure.defaults['gapw']

        # create unit cell outline
        s.last_direction = 240
        drawRingResonator(s,reslength,x,r,0,0,cap_length,0,0,cap_gap,cap_gap,pinw,stop_pinw,gapw,stop_gapw,5)
        s.last_direction = s.last_direction + 60 
        drawRingResonator(s,reslength,x,r,0,0,cap_length,0,0,cap_gap,cap_gap,pinw,stop_pinw,gapw,stop_gapw,5)
        s.last_direction = s.last_direction - 60 
        breakoffPt1 = s.last
        drawRingResonator(s,reslength,x,r,0,0,cap_length,0,0,cap_gap,cap_gap,pinw,stop_pinw,gapw,stop_gapw,5)
        s.last_direction = s.last_direction + 60
        breakoffPt2 = s.last
        drawRingResonator(s,reslength,x,r,0,0,cap_length,0,0,cap_gap,cap_gap,pinw,stop_pinw,gapw,stop_gapw,5)
        s.last_direction = s.last_direction - 60
        breakoffPt3 = s.last
        drawRingResonator(s,reslength,x,r,0,0,cap_length,0,0,cap_gap,cap_gap,pinw,stop_pinw,gapw,stop_gapw,5)
        s.last_direction = s.last_direction + 60
        breakoffPt4 = s.last
        drawRingResonator(s,reslength,x,r,0,0,cap_length,0,0,cap_gap,cap_gap,pinw,stop_pinw,gapw,stop_gapw,5)
        s.last_direction = s.last_direction + 60
        drawRingResonator(s,reslength,x,r,0,0,cap_length,0,0,cap_gap,cap_gap,pinw,stop_pinw,gapw,stop_gapw,5)

        if excludeLast==0:
            s.last_direction = s.last_direction + 60
            drawRingResonator(s,reslength,x,r,0,0,cap_length,0,0,cap_gap,cap_gap,pinw,stop_pinw,gapw,stop_gapw,5)

        #fill in the other resonators
        s.last = breakoffPt1
        s.last_direction = 0
        drawRingResonator(s,reslength,breakoffPt1[1]-breakoffPt3[1],r,-90,-90,cap_length,0,0,cap_gap,cap_gap,pinw,stop_pinw,gapw,stop_gapw,5)

        s.last = breakoffPt2
        s.last_direction = 180
        drawRingResonator(s,reslength,breakoffPt2[1]-breakoffPt4[1],r,90,90,cap_length,0,0,cap_gap,cap_gap,pinw,stop_pinw,gapw,stop_gapw,5)

class radialDraw:
    """A connection between a point of radius r and and point of radius R subtending an angle theta"""
    def __init__(self,structure,total_length,rStart,r,R,theta,extension_length,cap_length,cap_gap=None,pinw=None,stop_pinw=None,gapw=None,stop_gapw=None):
        """ """
        
        if total_length == 0:
            self.vlength=0
            return

        s=structure
        start=structure.last
                
        if pinw is None: pinw=structure.defaults['pinw']
        if gapw is None: gapw=structure.defaults['gapw']
        
        num_wiggles=4

        # Determines the direction of the turns
        if theta>0:
            asign=1
        else:
            asign=-1
            
        # Various parameters to get the lengths right
        a = rStart+cap_length+extension_length+(1+2*num_wiggles)*r
        L = abs(a*tan(theta))
        alpha = pi/2 + abs(theta)
        d = abs(r/(tan(alpha/2)))
        s1 = L+r-d
        s2 = R-d-a/cos(theta)-cap_length
        gamma = pi-alpha
        if gamma < -3*pi/2:
            gamma = gamma+2*pi
        sarc = gamma*r
        
        # Calculate vlength to make up the rest of the cavity length
        # vlength = (total_length-(cap_length-cap_gap/2)-(cap_length-cap_gap/2)-extension_length-pi*r*num_wiggles+2*r-s1-sarc-s2)/(2*num_wiggles)
        vlength = (total_length-2*(cap_length-cap_gap/2)-extension_length-pi*r/2-(num_wiggles-1)*pi*r-pi*r+2*r-s1-sarc-s2)/(2*(num_wiggles-1)+2)

        s.last = (start[0] + cap_length * cos(s.last_direction*pi/180), start[1] + cap_length * sin(s.last_direction*pi/180)) 
        s.last_direction = s.last_direction - 180;
        Inner_end_cap(s,cap_length,cap_gap,pinw,stop_pinw,gapw,stop_gapw)   
        s.last_direction = s.last_direction + 180;
        s.last = (start[0] + cap_length * cos(s.last_direction*pi/180), start[1] + cap_length * sin(s.last_direction*pi/180))  
        
        talliedLength = cap_length-cap_gap/2
        # print('added',cap_length-cap_gap/2)
        # print('talliedLength',talliedLength)   

        CPWStraight(s,extension_length,pinw,gapw)
        talliedLength = talliedLength + extension_length
        # print('added',extension_length)
        # print('talliedLength',talliedLength)

        CPWBend(s,asign*90,pinw,gapw,r)
        talliedLength = talliedLength + r*pi/2
        # print('added', r*pi/2)
        # print('talliedLength',talliedLength)
        
        for ii in range(num_wiggles-1):
            isign=2*(ii%2)-1
            CPWStraight(s,vlength,pinw,gapw)
            talliedLength = talliedLength + vlength
            # print('added', vlength)
            # print('talliedLength',talliedLength)

            CPWBend(s,isign*asign*180,pinw,gapw,r)
            talliedLength = talliedLength + r*pi
            # print('added', r*pi)
            # print('talliedLength',talliedLength)

            CPWStraight(s,vlength,pinw,gapw)
            talliedLength = talliedLength + vlength
            # print('added', vlength)
            # print('talliedLength',talliedLength)
        
        CPWStraight(s,vlength,pinw,gapw)
        talliedLength = talliedLength + vlength
        # print('added', vlength)
        # print('talliedLength',talliedLength)

        CPWBend(s,asign*180,pinw,gapw,r)   
        talliedLength = talliedLength + r*pi
        # print('added', r*pi)
        # print('talliedLength',talliedLength)

        CPWStraight(s,vlength-2*r,pinw,gapw)
        talliedLength = talliedLength + vlength-2*r
        # print('added', vlength-2*r)
        # print('talliedLength',talliedLength)
    
        CPWStraight(s,s1,pinw,gapw)
        talliedLength = talliedLength + s1
        # print('added', s1)
        # print('talliedLength',talliedLength)

        CPWBend(s,-asign*gamma*180/pi,pinw,gapw,r)
        talliedLength = talliedLength + r*gamma
        # print('added', pi*r)
        # print('talliedLength',talliedLength)

        CPWStraight(s,s2,pinw,gapw)
        talliedLength = talliedLength + s2
        # print('added', s2)
        # print('talliedLength',talliedLength)

        Inner_end_cap(s,cap_length,cap_gap,pinw,stop_pinw,gapw,stop_gapw) 
        talliedLength = talliedLength + cap_length-cap_gap/2
        # print('added', cap_length-cap_gap/2)
        # print('talliedLength',talliedLength)

