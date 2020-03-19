# -*- coding: utf-8 -*-
"""
Created on Sept 5, 2018

@author: Pranav Mundada, Mattias Fitzpatrick

"""

import os
import sys

pkgDir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.realpath(__file__))))

if not pkgDir in sys.path:
    sys.path.append(pkgDir)

import context

import DrawCodes.sdxf as sdxf
import random

#import ezdxf
from DrawCodes.MaskMakerPro import *
from math import sin,cos,pi,floor,asin,acos,tan,atan,sqrt

import re
from scipy import *
import pylab
import numpy
import time

from DrawCodes.alphanum import alphanum_dict
from random import randrange

import pickle
import datetime

from matplotlib.pyplot import *

# Code to define lattice info
from GraphCodes.EuclideanLayoutGenerator2 import EuclideanLayout

from GraphCodes.GeneralLayoutGenerator import *

from GraphCodes.LayoutGenerator5 import PlanarLayout

# Helper functions for drawing the graph
from DrawCodes.generalLayoutDrawHelpers import *


class layout(Chip):
    
    # chipSize = 24592.
    chipSize = 22650.
    #chipSize = 7000.
    
    def __init__(self,name,size=(chipSize,chipSize),mask_id_loc=(0,0),chip_id_loc=(0,0)):
        
        Chip.__init__(self,name,size,mask_id_loc,chip_id_loc)
        
        #####################################################################################################################
        "Make bounding box and define chip parameters"        
        #####################################################################################################################
        
        smallChipSize = 7000.
        smallChipGap = 350.
        innerEdgeBuffer = 175.
        outerEdgeBuffer = 300.

        "Creates the maskborder box"
        maskborder=350. #How much larger than the chip border        
        border=Structure(self,start=self.bottomleft_corner,color=5,layer="maskborder")
        box=[   (self.bottomleft_corner[0],self.bottomleft_corner[1]),
                (self.bottomright_corner[0],self.bottomright_corner[1]),
                (self.topright_corner[0],self.topright_corner[1]),
                (self.topleft_corner[0],self.topleft_corner[1]),
                (self.bottomleft_corner[0],self.bottomleft_corner[0])
                ]
        border.append(sdxf.PolyLine(box,layer=border.layer,color=border.color))   

        chipSpacing = smallChipSize + smallChipGap

        chipStarts = [(self.topleft_corner[0]+outerEdgeBuffer+innerEdgeBuffer,self.bottomleft_corner[1]+outerEdgeBuffer+innerEdgeBuffer),
                    # (self.topleft_corner[0]+outerEdgeBuffer+innerEdgeBuffer,self.bottomleft_corner[1]+outerEdgeBuffer+innerEdgeBuffer),
                    (self.bottomleft_corner[0]+outerEdgeBuffer+innerEdgeBuffer+chipSpacing,self.bottomleft_corner[1]+outerEdgeBuffer+innerEdgeBuffer),
                    # (self.bottomleft_corner[0]+outerEdgeBuffer+innerEdgeBuffer+chipSpacing,self.bottomleft_corner[1]+outerEdgeBuffer+innerEdgeBuffer),
                    (self.bottomleft_corner[0]+outerEdgeBuffer+innerEdgeBuffer+2*chipSpacing,self.bottomleft_corner[1]+outerEdgeBuffer+innerEdgeBuffer),
                    # (self.bottomleft_corner[0]+outerEdgeBuffer+innerEdgeBuffer+2*chipSpacing,self.bottomleft_corner[1]+outerEdgeBuffer+innerEdgeBuffer),
                    (self.bottomleft_corner[0]+outerEdgeBuffer+innerEdgeBuffer,self.bottomleft_corner[1]+outerEdgeBuffer+innerEdgeBuffer+chipSpacing),
                    (self.bottomleft_corner[0]+outerEdgeBuffer+innerEdgeBuffer+chipSpacing,self.bottomleft_corner[1]+outerEdgeBuffer+innerEdgeBuffer+chipSpacing),
                    (self.bottomleft_corner[0]+outerEdgeBuffer+innerEdgeBuffer+2*chipSpacing,self.bottomleft_corner[1]+outerEdgeBuffer+innerEdgeBuffer+chipSpacing),
                    (self.bottomleft_corner[0]+outerEdgeBuffer+innerEdgeBuffer,self.bottomleft_corner[1]+outerEdgeBuffer+innerEdgeBuffer+2*chipSpacing),
                    (self.bottomleft_corner[0]+outerEdgeBuffer+innerEdgeBuffer+chipSpacing,self.bottomleft_corner[1]+outerEdgeBuffer+innerEdgeBuffer+2*chipSpacing),
                    (self.bottomleft_corner[0]+outerEdgeBuffer+innerEdgeBuffer+2*chipSpacing,self.bottomleft_corner[1]+outerEdgeBuffer+innerEdgeBuffer+2*chipSpacing)]

        for startPos in chipStarts:
            border=Structure(self,start=startPos,color=3,layer="chipborder")
            box=[   (startPos[0],startPos[1]),
                    (startPos[0]+smallChipSize,startPos[1]),
                    (startPos[0]+smallChipSize,startPos[1]+smallChipSize),
                    (startPos[0],startPos[1]+smallChipSize),
                    (startPos[0],startPos[1])]
            border.append(sdxf.PolyLine(box,layer=border.layer,color=border.color))
        




        pos=[0,0]    #initial position of device
        chipborder=24592  
        chipX=chipborder #Length of device
        chipY=chipborder #Width of device
        cut_border=chipborder + 350
        edgeBufferDistance = 850

    
        """
        CPW Parameters: defines the basic geometry of the resonators 
        """  
        pinw_Z0 = 11.3
        gapw_Z0 = 5.54
        scalefactor=1.
        pinwLinear=scalefactor*pinw_Z0
        gapwLinear=scalefactor*gapw_Z0
        x = 3000.
        r = 90.
        num_wiggles = 4
        up = 1
        

        "Cap Parameters"     
        cap_length = 80.0
        cap_gap = 2.0

        scalefactor=1.5
        pinw=scalefactor*pinw_Z0
        gapw=scalefactor*gapw_Z0
        stop_gapw=10*gapw
        stop_pinw=8*pinw
        gapw_buffer = 0.0
        cap_gap_out = 100        
        cap_gap_ext = 0 #This values increases the gap for the outer capacitors
        
        "Bondpad Parameters"        
        bond_pad_length=350. #Length of Rectangular portion of bond pad
        launcher_scalefactor = 30
        launcher_pinw=launcher_scalefactor*pinw_Z0  
        launcher_gapw=launcher_scalefactor*gapw_Z0      
        taper_length= 300.
        launcher_padding = launcher_gapw

        stub_length = 100.
        extension_length = 30.
        chipCornerToBondPad = [1200,1200]
        
        defaultRadius = 100. #the name of this variable is important. Do not change it. It functions as a global.

        def drawHanger(structure, start, startDirection, leftRight, pinw, gapw, couplingGap, couplingStraight, resonatorStraight1, resonatorStraight2, resonatorStraight3, targetLength, upDown = 1):
            s = structure
            s.last = start
            s.last_direction = startDirection

            CPWStraight(s,couplingGap,0,gapw+pinw/2)
            CPWStraight(s,couplingStraight, pinw,gapw)
            CPWBend(s,-leftRight*startDirection,pinw,gapw,radius)
            br_base = 30
            br_width = 50
            CPWStraight(mask,resonatorStraight1,pinw,gapw) 
            accumulatedLength = couplingStraight + (pi/2)*radius + resonatorStraight1
            isign = -1
            if upDown == 1:
                usign = 1
            if upDown == -1:
                usign = -1;
            while accumulatedLength<targetLength:
                if (accumulatedLength+pi*radius/2)>targetLength:
                    angle = usign*isign*((targetLength-accumulatedLength)/radius)* 180/pi
                else:
                    angle = usign*isign*90
                CPWBend(mask,angle,pinw,gapw,radius)
                accumulatedLength = accumulatedLength + radius*abs(angle*pi/180)
                if (accumulatedLength+resonatorStraight3)>targetLength:
                    delta = targetLength-accumulatedLength
                else:
                    delta = resonatorStraight3
                CPWStraight(mask,delta,pinw,gapw)
                accumulatedLength = accumulatedLength + delta
                if (accumulatedLength+pi*radius/2)>targetLength:
                    angle = usign*isign*((targetLength-accumulatedLength)/radius)* 180/pi
                else:
                    angle = usign*isign*90
                CPWBend(mask,angle,pinw,gapw,radius)
                accumulatedLength = accumulatedLength + radius*abs(angle*pi/180)
                print(('after angle:', angle, 'accumulatedLength',accumulatedLength))

                isign = isign*-1
                if (accumulatedLength+resonatorStraight2)>targetLength:
                    delta = targetLength-accumulatedLength
                else:
                    delta = resonatorStraight2
                CPWStraight(mask,delta,pinw,gapw) 
                accumulatedLength = accumulatedLength + delta
                print(('after delta:', delta, 'accumulatedLength',accumulatedLength))
            print('             ')
            
        def drawMixedRadiusHanger(structure, start, startDirection, leftRight, pinw, gapw, couplingGap, couplingStraight, resonatorStraight1, resonatorStraight2, resonatorStraight3, targetLength, secondRadius = defaultRadius, upDown = 1):
            s = structure
            s.last = start
            s.last_direction = startDirection

            CPWStraight(s,couplingGap,0,gapw+pinw/2)
            CPWStraight(s,couplingStraight, pinw,gapw)
            CPWBend(s,-leftRight*startDirection,pinw,gapw,defaultRadius)
            br_base = 30
            br_width = 50
            CPWStraight(mask,resonatorStraight1,pinw,gapw) 
            accumulatedLength = couplingStraight + (pi/2)*defaultRadius + resonatorStraight1
            isign = -1
            if upDown == 1:
                usign = 1
            if upDown == -1:
                usign = -1;
            
            #swtich to new radius
            while accumulatedLength<targetLength:
                if (accumulatedLength+pi*secondRadius/2)>targetLength:
                    angle = usign*isign*((targetLength-accumulatedLength)/secondRadius)* 180/pi
                else:
                    angle = usign*isign*90
                CPWBend(mask,angle,pinw,gapw,secondRadius)
                accumulatedLength = accumulatedLength + secondRadius*abs(angle*pi/180)
                if (accumulatedLength+resonatorStraight3)>targetLength:
                    delta = targetLength-accumulatedLength
                else:
                    delta = resonatorStraight3
                CPWStraight(mask,delta,pinw,gapw)
                accumulatedLength = accumulatedLength + delta
                if (accumulatedLength+pi*secondRadius/2)>targetLength:
                    angle = usign*isign*((targetLength-accumulatedLength)/secondRadius)* 180/pi
                else:
                    angle = usign*isign*90
                CPWBend(mask,angle,pinw,gapw,secondRadius)
                accumulatedLength = accumulatedLength + secondRadius*abs(angle*pi/180)
                print(('after angle:', angle, 'accumulatedLength',accumulatedLength))

                isign = isign*-1
                if (accumulatedLength+resonatorStraight2)>targetLength:
                    delta = targetLength-accumulatedLength
                else:
                    delta = resonatorStraight2
                CPWStraight(mask,delta,pinw,gapw) 
                accumulatedLength = accumulatedLength + delta
                print(('after delta:', delta, 'accumulatedLength',accumulatedLength))
            print('             ')
            
        def drawTightHanger(structure, start, startDirection, leftRight, pinw, gapw, couplingGap, couplingStraight, resonatorStraight1, meanderSize, turnEdge, targetLength, hangerRadius = 100):
            s = structure
            s.last = start
            s.last_direction = startDirection

            CPWStraight(s,couplingGap,0,gapw+pinw/2)
            CPWStraight(s,couplingStraight, pinw,gapw)
            CPWBend(s,-leftRight*startDirection,pinw,gapw,hangerRadius)
            br_base = 30
            br_width = 50
            CPWStraight(mask,resonatorStraight1,pinw,gapw) 
            accumulatedLength = couplingStraight + (pi/2)*hangerRadius + resonatorStraight1
            isign = -1
            
            ind = 1
            while accumulatedLength<targetLength:
                print(ind)
                #make the first turn to parallel to the feedline
                if (accumulatedLength+pi*hangerRadius/2)>targetLength:
                    angle = isign*((targetLength-accumulatedLength)/hangerRadius)* 180/pi
                else:
                    angle = isign*90
                CPWBend(mask,angle,pinw,gapw,hangerRadius)
                accumulatedLength = accumulatedLength + hangerRadius*abs(angle*pi/180)
                
                #start the meander
                if (accumulatedLength+meanderSize)>targetLength:
                    delta = targetLength-accumulatedLength
                else:
                    delta = meanderSize
                CPWStraight(mask,delta,pinw,gapw)
                accumulatedLength = accumulatedLength + delta
                
                #turn at the end of the first meander
                if (accumulatedLength+pi*hangerRadius/2)>targetLength:
                    angle = -isign*((targetLength-accumulatedLength)/hangerRadius)* 180/pi
                else:
                    angle = -isign*90
                CPWBend(mask,angle,pinw,gapw,hangerRadius)
                accumulatedLength = accumulatedLength + hangerRadius*abs(angle*pi/180)
                print(('after angle:', angle, 'accumulatedLength',accumulatedLength))

                #do the straight in the meander
                isign = isign*-1
                if (accumulatedLength+turnEdge)>targetLength:
                    delta = targetLength-accumulatedLength
                else:
                    delta = turnEdge
                CPWStraight(mask,delta,pinw,gapw) 
                accumulatedLength = accumulatedLength + delta
                print(('after delta:', delta, 'accumulatedLength',accumulatedLength))
                
                ind = ind +1
            print('             ')
            
        def drawTTestTightHanger(structure, start, startDirection, leftRight, pinw, gapw, couplingGap, couplingStraight, resonatorStraight1, meanderSize, turnEdge, targetLength, hangerRadius = 100, finalStraight = 300, upDown = 1):
            s = structure
            s.last = start
            s.last_direction = startDirection

            CPWStraight(s,couplingGap,0,gapw+pinw/2)
            CPWStraight(s,couplingStraight, pinw,gapw)
            CPWBend(s,-leftRight*startDirection,pinw,gapw,hangerRadius)
            br_base = 30
            br_width = 50
            CPWStraight(mask,resonatorStraight1,pinw,gapw) 
            accumulatedLength = couplingStraight + (pi/2)*hangerRadius + resonatorStraight1
            isign = -1
            
            if upDown == 1:
                usign = 1
            if upDown == -1:
                usign = -1;
            
            ind = 1
            
            #leave room for the final turn
            finalAngle = 30
            targetLength = targetLength - hangerRadius*abs(finalAngle*pi/180) - finalStraight
            
            while accumulatedLength<targetLength:
                print(ind)
                #make the first turn to parallel to the feedline
                if (accumulatedLength+pi*hangerRadius/2)>targetLength:
                    angle = usign*isign*((targetLength-accumulatedLength)/hangerRadius)* 180/pi
                else:
                    angle = usign*isign*90
                CPWBend(mask,angle,pinw,gapw,hangerRadius)
                accumulatedLength = accumulatedLength + hangerRadius*abs(angle*pi/180)
                
                #start the meander
                if (accumulatedLength+meanderSize)>targetLength:
                    delta = targetLength-accumulatedLength
                else:
                    delta = meanderSize
                CPWStraight(mask,delta,pinw,gapw)
                accumulatedLength = accumulatedLength + delta
                
                #turn at the end of the first meander
                if (accumulatedLength+pi*hangerRadius/2)>targetLength:
                    angle = -usign*isign*((targetLength-accumulatedLength)/hangerRadius)* 180/pi
                else:
                    angle = -usign*isign*90
                CPWBend(mask,angle,pinw,gapw,hangerRadius)
                accumulatedLength = accumulatedLength + hangerRadius*abs(angle*pi/180)
                print(('after angle:', angle, 'accumulatedLength',accumulatedLength))

                #do the straight in the meander
                isign = isign*-1
                if (accumulatedLength+turnEdge)>targetLength:
                    delta = targetLength-accumulatedLength
                else:
                    delta = turnEdge
                CPWStraight(mask,delta,pinw,gapw) 
                accumulatedLength = accumulatedLength + delta
                print(('after delta:', delta, 'accumulatedLength',accumulatedLength))
                
                ind = ind +1
                
            #draw the last turn and funal straight
            angle = usign*finalAngle*isign
            CPWBend(mask,angle,pinw,gapw,hangerRadius)
            accumulatedLength = accumulatedLength + hangerRadius*abs(angle*pi/180)
            
            delta = finalStraight
            CPWStraight(mask,delta,pinw,gapw)
            accumulatedLength = accumulatedLength + delta
            
            print('             ')

        def drawHanger_with_bridges(structure, start, startDirection, leftRight, pinw, gapw, couplingGap, \
            couplingStraight, resonatorStraight1, resonatorStraight2, resonatorStraight3, targetLength, br_width, br_base,bridge_gap=200):
            s = structure
            s.last = start
            s.last_direction = startDirection

            CPWStraight(s,couplingGap,0,gapw+pinw/2)
            CPWStraight(s,couplingStraight, pinw,gapw)
            CPWBend(s,-leftRight*startDirection,pinw,gapw,radius)
            CPWStraight_Bridges_Layer1(mask,resonatorStraight1, br_base, br_width, pinw,gapw,bridge_gap=bridge_gap) 
            accumulatedLength = couplingStraight + (pi/2)*radius + resonatorStraight1
            isign = -1
            while accumulatedLength<targetLength:
                if (accumulatedLength+pi*radius/2)>targetLength:
                    angle = isign*((targetLength-accumulatedLength)/radius)* 180/pi
                else:
                    angle = isign*90
                CPWBend(mask,angle,pinw,gapw,radius)
                accumulatedLength = accumulatedLength + radius*abs(angle*pi/180)
                if (accumulatedLength+resonatorStraight3)>targetLength:
                    delta = targetLength-accumulatedLength
                else:
                    delta = resonatorStraight3
                CPWStraight_Bridges_Layer1(mask,delta, br_base, br_width, pinw,gapw,bridge_gap=bridge_gap) 
                accumulatedLength = accumulatedLength + delta
                if (accumulatedLength+pi*radius/2)>targetLength:
                    angle = isign*((targetLength-accumulatedLength)/radius)* 180/pi
                else:
                    angle = isign*90
                CPWBend(mask,angle,pinw,gapw,radius)
                accumulatedLength = accumulatedLength + radius*abs(angle*pi/180)
                print(('after angle:', angle, 'accumulatedLength',accumulatedLength))

                isign = isign*-1
                if (accumulatedLength+resonatorStraight2)>targetLength:
                    delta = targetLength-accumulatedLength
                else:
                    delta = resonatorStraight2
                CPWStraight_Bridges_Layer1(mask,delta, br_base, br_width, pinw,gapw,bridge_gap=bridge_gap)  
                accumulatedLength = accumulatedLength + delta
                print(('after delta:', delta, 'accumulatedLength',accumulatedLength))
            print('             ')

        #####################################################################################################################
        "Start laying out resonators"
        #####################################################################################################################
        
        mask=Structure(self,start=0)
        mask.defaults = {'pinw':pinw,'gapw':gapw,'bendradius':r} #Set these as defaults for simplified inherited input
   
        # Draw some alignment marks
        edgeBuffer = 530
        # drawAlignmentMarks(mask,20,[520,520],edgeBuffer,chipborder-300,pos)
       
        mid=(1250,-7000+1400-200) 

        mask.last = mid

        # Parameters for Bond Pads
        edgeOffset = 1255
        sideShift = 2158
        centerShift = 11866
        launcher_padding = launcher_gapw
        

        couplingSpacing = 2.8 #smaller than CW gap
#        couplingSpacing = gapw_Z0 #CPW gap I think. But this makes the gap double
        couplingSpacing = 0.
        couplingGap = 40.
        couplingStraight = 80.
        resonatorSpacing = 500.
        resonatorStraight1 = 1700.
        resonatorStraight2 = 1350.
        targetLength = 4000.
        radius = defaultRadius #the name of this variable is important. Do not change it. It functions as a global.

        showLattice=0

        scalefactorSmall=2.0/10
        pinwSmall=scalefactorSmall*pinw_Z0
        gapwSmall=scalefactorSmall*gapw_Z0
        
        feedlineRadius = 120.
        feedlinePinw = pinwLinear
        feedlineGapw = gapwLinear

        scalefactorLarge=3.0/10
        pinwLarge=scalefactorLarge*pinw_Z0
        gapwLarge=scalefactorLarge*gapw_Z0        

        for idx in range(9):
            # startVertices = bondPadConnectionParameters['startVertices']
            
            bufferDistance = 1;
            
            # define potential bond pad positions
            bondPadConnectionParameters = {'startVertices':[0,1],'resonatorIndices':[0,2],'stopBondPad':[0,1],'directLines':[1,1],'pinwLinear':pinwLinear,'gapwLinear':gapwLinear}
            stopBondPads = bondPadConnectionParameters['stopBondPad']
            directLines = bondPadConnectionParameters['directLines']
            startResonators = bondPadConnectionParameters['resonatorIndices']
            lowerLeftChipCorner = chipStarts[idx]
            bondPadPositions = [[lowerLeftChipCorner[0]+chipCornerToBondPad[0],lowerLeftChipCorner[1]+chipCornerToBondPad[1]],
                                [lowerLeftChipCorner[0]+smallChipSize-chipCornerToBondPad[0],lowerLeftChipCorner[1]+chipCornerToBondPad[1]],
                                [lowerLeftChipCorner[0]+smallChipSize-chipCornerToBondPad[0],lowerLeftChipCorner[1]+smallChipSize-chipCornerToBondPad[1]],
                                [lowerLeftChipCorner[0]+chipCornerToBondPad[0],lowerLeftChipCorner[1]+smallChipSize-chipCornerToBondPad[1]]]
            bondPadAngles = [90,90,270,270]
            
            for idy in range(len(bondPadPositions)):

                # find where that bond pad is
                bondPadPos = bondPadPositions[idy]
                bondPadAng = bondPadAngles[idy]
                bondPadEndPos = [bondPadPos[0]+(bond_pad_length+launcher_padding+taper_length)*numpy.cos(bondPadAng*numpy.pi/180.),bondPadPos[1]+(bond_pad_length+launcher_padding+taper_length)*numpy.sin(bondPadAng*numpy.pi/180.)]

                # figure out if that bond pad is connected
                if idy in stopBondPads:

                    drawBondPad(mask,bondPadPos,bondPadAng,pinwLinear,gapwLinear,bond_pad_length,launcher_pinw,launcher_gapw,taper_length,launcher_padding)
                    if idy == stopBondPads[0]:
                        current_pos = mask.last
                    index = stopBondPads.index(idy)
                    # print('idx',idx,'stopBondPads',stopBondPads[index],'startVertex',startVertex)
                    
                    
            #######################        
            #do each individual 7x7
            #######################
            
            def draw_ParanvOriginal():
                #Pranav's original
                
                #draw the feedline
                startPtBondPad = bondPadPositions[stopBondPads[-1]]
                mask.last = current_pos
                endPtBondPad = mask.last
                i = stopBondPads[0]
                j = stopBondPads[1]
                CPWStraight(mask,4000,feedlinePinw,feedlineGapw)
                CPWBend(mask,-90,feedlinePinw,feedlineGapw,feedlineRadius)  
                CPWStraight(mask,bondPadPositions[j][0] - bondPadPositions[i][0]-2*feedlineRadius,feedlinePinw,feedlineGapw)
                CPWBend(mask,-90,feedlinePinw,feedlineGapw,feedlineRadius)  
                CPWStraight(mask,4000,feedlinePinw,feedlineGapw)
                mask.last = (mask.last[0], mask.last[1]-(endPtBondPad[1]-startPtBondPad[1]))
                
                
                #setup for the hangers
                leftRight = 1
                scalefactor= 1.
                pinw=scalefactor*pinw_Z0
                gapw=scalefactor*gapw_Z0
                startDirection = -90
                start=  (endPtBondPad[0]+(feedlinePinw/2+feedlineGapw+pinw/2+couplingSpacing),endPtBondPad[1]+resonatorSpacing+1100)
                
                
                #resonator 1 : lower left
                resonatorStraight3 = 400.
                resonatorStraight1 = 900.
                resonatorStraight2 = 550.
                resLength = 3400.
                drawHanger(mask, start, startDirection, leftRight, pinw, gapw, couplingGap, couplingStraight, resonatorStraight1, resonatorStraight2, resonatorStraight3, resLength)

                #resonator2: clockwise from 1
                resLength = 3150.
                resonatorStraight3 = 0.
                resonatorStraight1 = 1050.
                resonatorStraight2 = 650.
                start=  (endPtBondPad[0]+(feedlinePinw/2+feedlineGapw+pinw/2+couplingSpacing),endPtBondPad[1]+resonatorSpacing+2100)
                drawHanger(mask, start, startDirection, leftRight, pinw, gapw, couplingGap, couplingStraight, resonatorStraight1, resonatorStraight2, resonatorStraight3, resLength)
                
                #resonator3: clockwise from 2 
                resLength = 3200.
                resonatorStraight3 = 0.
                resonatorStraight1 = 1700.
                resonatorStraight2 = 1350.
                start=  (endPtBondPad[0]+(feedlinePinw/2+feedlineGapw+pinw/2+couplingSpacing),endPtBondPad[1]+resonatorSpacing+3100)
                drawHanger(mask, start, startDirection, leftRight, pinw, gapw, couplingGap, couplingStraight, resonatorStraight1, resonatorStraight2, resonatorStraight3, resLength)

                #resonator4: upper right
                resLength = 3250.
                leftRight = -1
                resonatorStraight3 = 400.
                resonatorStraight1 = 1350.
                resonatorStraight2 = 1000.
                start=  (endPtBondPad[0]-(feedlinePinw/2+feedlineGapw+pinw/2)+bondPadPositions[j][0] - bondPadPositions[i][0]-couplingSpacing,endPtBondPad[1]+resonatorSpacing+3*900+100)
                drawHanger(mask, start, startDirection, leftRight, pinw, gapw, couplingGap, couplingStraight, resonatorStraight1, resonatorStraight2, resonatorStraight3, resLength)

                #resonator5: clockwise from 4
                resonatorStraight3 = 600
                resLength = 3300.
                resonatorStraight3 = 400.
                resonatorStraight1 = 1700.
                resonatorStraight2 = 1350.
                start=  (endPtBondPad[0]-(feedlinePinw/2+feedlineGapw+pinw/2)+bondPadPositions[j][0] - bondPadPositions[i][0]-couplingSpacing,endPtBondPad[1]+resonatorSpacing+2*900)
                drawHanger(mask, start, startDirection, leftRight, pinw, gapw, couplingGap, couplingStraight, resonatorStraight1, resonatorStraight2, resonatorStraight3, resLength)

                #resonator6:
                resonatorStraight3 = 400
                resLength = 3350.
                resonatorStraight3 = 400.
                resonatorStraight1 = 2600.
                resonatorStraight2 = 1350.
                start=  (endPtBondPad[0]-(feedlinePinw/2+feedlineGapw+pinw/2)+bondPadPositions[j][0] - bondPadPositions[i][0]-couplingSpacing,endPtBondPad[1]+resonatorSpacing+800)
                drawHanger(mask, start, startDirection, leftRight, pinw, gapw, couplingGap, couplingStraight, resonatorStraight1, resonatorStraight2, resonatorStraight3, resLength)

                #resonator7:
                resonatorStraight3 = 400
                resLength = 3100.
                resonatorStraight3 = 0.
                resonatorStraight1 = 3100.
                resonatorStraight2 = 0.
                start=  (endPtBondPad[0]-(feedlinePinw/2+feedlineGapw+pinw/2)+bondPadPositions[j][0] - bondPadPositions[i][0]-couplingSpacing,endPtBondPad[1]+resonatorSpacing+100)
                drawHanger(mask, start, startDirection, leftRight, pinw, gapw, couplingGap, couplingStraight, resonatorStraight1, resonatorStraight2, resonatorStraight3, resLength)
                
                return
            
            
            def draw_BendTestChip(step = 150):
                #bend test chip
#                if idx == 0:
#                    step = 150
#                if idx == 1:
#                    step = 300
#                step = 150
                
                #draw the feedline
                startPtBondPad = bondPadPositions[stopBondPads[-1]]
                mask.last = current_pos
                endPtBondPad = mask.last
                i = stopBondPads[0]
                j = stopBondPads[1]
                CPWStraight(mask,4000,feedlinePinw,feedlineGapw)
                CPWBend(mask,-90,feedlinePinw,feedlineGapw,feedlineRadius)  
                CPWStraight(mask,bondPadPositions[j][0] - bondPadPositions[i][0]-2*feedlineRadius,feedlinePinw,feedlineGapw)
                CPWBend(mask,-90,feedlinePinw,feedlineGapw,feedlineRadius)  
                CPWStraight(mask,4000,feedlinePinw,feedlineGapw)
                mask.last = (mask.last[0], mask.last[1]-(endPtBondPad[1]-startPtBondPad[1]))
                
                
                #setup for the hangers
                leftRight = 1
                scalefactor= 1.
                pinw=scalefactor*pinw_Z0
                gapw=scalefactor*gapw_Z0
                startDirection = -90
                start=  (endPtBondPad[0]+(feedlinePinw/2+feedlineGapw+pinw/2+couplingSpacing),endPtBondPad[1]+resonatorSpacing+1100)
                
                
                baseLength = 7000/2.
                jog = 400
                
                #resonator 1 : lower left
                resLength = baseLength +4*step
                leftRight = 1
                if step == 150:
                    resonatorStraight1 = 900.
                    resonatorStraight2 = 550.
                    resonatorStraight3 = jog
                else:
                    resonatorStraight1 = 900.
                    resonatorStraight2 = 550.
                    resonatorStraight3 = jog
                drawHanger(mask, start, startDirection, leftRight, pinw, gapw, couplingGap, couplingStraight, resonatorStraight1, resonatorStraight2, resonatorStraight3, resLength)

#                #resonator2: clockwise from 1
                resLength = baseLength + step
                leftRight = 1
                if step == 150:
                    resonatorStraight1 = 1400.
                    resonatorStraight2 = 800.
                    resonatorStraight3 = jog
                else:
                    resonatorStraight1 = 1400.
                    resonatorStraight2 = 800.
                    resonatorStraight3 = jog + 100
                start=  (endPtBondPad[0]+(feedlinePinw/2+feedlineGapw+pinw/2+couplingSpacing),endPtBondPad[1]+resonatorSpacing+2900)
                drawHanger(mask, start, startDirection, leftRight, pinw, gapw, couplingGap, couplingStraight, resonatorStraight1, resonatorStraight2, resonatorStraight3, resLength)

                #resonator5: clockwise from 4
                resLength = baseLength + 2*step
                leftRight = -1
                if step == 150:
                    resonatorStraight1 = 1600.
                    resonatorStraight2 = resLength - resonatorStraight1 - jog
                    resonatorStraight3 = jog
                else:
                    resonatorStraight1 = 1750.
                    resonatorStraight2 = resLength - resonatorStraight1 - jog
                    resonatorStraight3 = jog
                start=  (endPtBondPad[0]-(feedlinePinw/2+feedlineGapw+pinw/2)+bondPadPositions[j][0] - bondPadPositions[i][0]-couplingSpacing,endPtBondPad[1]+resonatorSpacing+2400)
                drawHanger(mask, start, startDirection, leftRight, pinw, gapw, couplingGap, couplingStraight, resonatorStraight1, resonatorStraight2, resonatorStraight3, resLength)

                #resonator6:
                resLength = baseLength + 3*step
                leftRight = -1
                resonatorStraight1 = 2250.
                resonatorStraight2 = resLength - resonatorStraight1 - jog
                resonatorStraight3 = 1800
                start=  (endPtBondPad[0]-(feedlinePinw/2+feedlineGapw+pinw/2)+bondPadPositions[j][0] - bondPadPositions[i][0]-couplingSpacing,endPtBondPad[1]+resonatorSpacing+1300)
                drawHanger(mask, start, startDirection, leftRight, pinw, gapw, couplingGap, couplingStraight, resonatorStraight1, resonatorStraight2, resonatorStraight3, resLength)

                #resonator7:
                resLength = baseLength
                leftRight = -1
                if step == 150:
                    resonatorStraight1 = resLength
                    resonatorStraight2 = 0.
                    resonatorStraight3 = 0.
                else:
                    resonatorStraight1 = resLength
                    resonatorStraight2 = 0.
                    resonatorStraight3 = 0.
                start=  (endPtBondPad[0]-(feedlinePinw/2+feedlineGapw+pinw/2)+bondPadPositions[j][0] - bondPadPositions[i][0]-couplingSpacing,endPtBondPad[1]+resonatorSpacing+100 + 100)
                drawHanger(mask, start, startDirection, leftRight, pinw, gapw, couplingGap, couplingStraight, resonatorStraight1, resonatorStraight2, resonatorStraight3, resLength)
            
                return
            
            def draw_ParallelStraighTestChip_twoRadii(step = 150):
                    #parallel straight test chip
    #                if idx == 3:
    #                    step = 300
    #                if idx == 4:
    #                    step = 250
    #                if idx == 5:
    #                    step = 200
#                    step = 200
#                    step = 150
                    
                    #draw the feedline
                    startPtBondPad = bondPadPositions[stopBondPads[-1]]
                    mask.last = current_pos
                    endPtBondPad = mask.last
                    i = stopBondPads[0]
                    j = stopBondPads[1]
                    CPWStraight(mask,4000,feedlinePinw,feedlineGapw)
                    CPWBend(mask,-90,feedlinePinw,feedlineGapw,feedlineRadius)  
                    CPWStraight(mask,bondPadPositions[j][0] - bondPadPositions[i][0]-2*feedlineRadius,feedlinePinw,feedlineGapw)
                    CPWBend(mask,-90,feedlinePinw,feedlineGapw,feedlineRadius)  
                    CPWStraight(mask,4000,feedlinePinw,feedlineGapw)
                    mask.last = (mask.last[0], mask.last[1]-(endPtBondPad[1]-startPtBondPad[1]))
                    
                    
                    #setup for the hangers
                    leftRight = 1
                    scalefactor= 1.
                    pinw=scalefactor*pinw_Z0
                    gapw=scalefactor*gapw_Z0
                    startDirection = -90
                    start=  (endPtBondPad[0]+(feedlinePinw/2+feedlineGapw+pinw/2+couplingSpacing),endPtBondPad[1]+resonatorSpacing+1100)
                    
                    
                    baseLength = 7000/2.
                    jog = 400
                    
                    #resonator 1 : lower left
                    resLength = baseLength
                    leftRight = 1
                    meanderSize = 300.
                    turnEdge = 0.
                    nbends = 5
                    resonatorStraight1 = resLength -couplingStraight - radius * numpy.pi/2 - (nbends-1)*turnEdge - nbends*meanderSize - nbends*2*radius * numpy.pi/2
                    start=  (endPtBondPad[0]+(feedlinePinw/2+feedlineGapw+pinw/2+couplingSpacing),endPtBondPad[1]+resonatorSpacing+400)
                    drawTightHanger(mask, start, startDirection, leftRight, pinw, gapw, couplingGap, couplingStraight, resonatorStraight1, meanderSize, turnEdge, resLength, hangerRadius = radius)
                    
                    #resonator 2 : middle left
                    resLength = baseLength +step
                    leftRight = 1
                    meanderSize = 300.
                    turnEdge = 12.5
                    nbends = 5
                    resonatorStraight1 = resLength -couplingStraight - radius * numpy.pi/2 - (nbends-1)*turnEdge - nbends*meanderSize - nbends*2*radius * numpy.pi/2
                    start=  (endPtBondPad[0]+(feedlinePinw/2+feedlineGapw+pinw/2+couplingSpacing),endPtBondPad[1]+resonatorSpacing+1400)
                    drawTightHanger(mask, start, startDirection, leftRight, pinw, gapw, couplingGap, couplingStraight, resonatorStraight1, meanderSize, turnEdge, resLength, hangerRadius = radius)
                    
                    #resonator 3 : upper left
                    resLength = baseLength +2*step
                    leftRight = 1
                    meanderSize = 300.
                    turnEdge = 25.
                    nbends = 5
                    resonatorStraight1 = resLength -couplingStraight - radius * numpy.pi/2 - (nbends-1)*turnEdge - nbends*meanderSize - nbends*2*radius * numpy.pi/2
                    start=  (endPtBondPad[0]+(feedlinePinw/2+feedlineGapw+pinw/2+couplingSpacing),endPtBondPad[1]+resonatorSpacing+2400)
                    drawTightHanger(mask, start, startDirection, leftRight, pinw, gapw, couplingGap, couplingStraight, resonatorStraight1, meanderSize, turnEdge, resLength, hangerRadius = radius)
                    
                    #resonator 4 : upper left
                    resLength = baseLength +3*step
                    leftRight = 1
                    meanderSize = 300.
                    turnEdge = 50.
                    nbends = 5
                    resonatorStraight1 = resLength -couplingStraight - radius * numpy.pi/2 - (nbends-1)*turnEdge - nbends*meanderSize - nbends*2*radius * numpy.pi/2
                    start=  (endPtBondPad[0]+(feedlinePinw/2+feedlineGapw+pinw/2+couplingSpacing),endPtBondPad[1]+resonatorSpacing+3400)
                    drawTightHanger(mask, start, startDirection, leftRight, pinw, gapw, couplingGap, couplingStraight, resonatorStraight1, meanderSize, turnEdge, resLength, hangerRadius = radius)
                    
                    
                    #resonator 5 : upper right
                    resLength = baseLength +4*step
                    leftRight = -1
                    meanderSize = 300.
                    nbends = 5.
                    turnEdge = 100.
                    resonatorStraight1 = resLength -couplingStraight - radius * numpy.pi/2 - (nbends-1)*turnEdge - nbends*meanderSize - nbends*2*radius * numpy.pi/2
                    if step == 150:
                        start=  (endPtBondPad[0]-(feedlinePinw/2+feedlineGapw+pinw/2)+bondPadPositions[j][0] - bondPadPositions[i][0]-couplingSpacing,endPtBondPad[1]+resonatorSpacing+2900)
                    else:
                        start=  (endPtBondPad[0]-(feedlinePinw/2+feedlineGapw+pinw/2)+bondPadPositions[j][0] - bondPadPositions[i][0]-couplingSpacing,endPtBondPad[1]+resonatorSpacing+2200-125)
                    drawTightHanger(mask, start, startDirection, leftRight, pinw, gapw, couplingGap, couplingStraight, resonatorStraight1, meanderSize, turnEdge, resLength, hangerRadius = radius)
                    
                    #resonator 6 : middle right
                    resLength = baseLength +5*step
                    leftRight = -1
                    meanderSize = 300.
                    nbends = 5
                    if step == 150:
                        turnEdge = 200.
                        resonatorStraight1 = resLength -couplingStraight - radius * numpy.pi/2 - (nbends-1)*turnEdge - nbends*meanderSize - nbends*2*radius * numpy.pi/2
                        start=  (endPtBondPad[0]-(feedlinePinw/2+feedlineGapw+pinw/2)+bondPadPositions[j][0] - bondPadPositions[i][0]-couplingSpacing,endPtBondPad[1]+resonatorSpacing+1900)
                    else:
                        turnEdge = 200.
                        resonatorStraight1 = resLength -couplingStraight - radius * numpy.pi/2 - (nbends-1)*turnEdge - nbends*meanderSize - nbends*2*radius * numpy.pi/2
                        start=  (endPtBondPad[0]-(feedlinePinw/2+feedlineGapw+pinw/2)+bondPadPositions[j][0] - bondPadPositions[i][0]-couplingSpacing,endPtBondPad[1]+resonatorSpacing+1200-125)
                    drawTightHanger(mask, start, startDirection, leftRight, pinw, gapw, couplingGap, couplingStraight, resonatorStraight1, meanderSize, turnEdge, resLength, hangerRadius = radius)
                    
                    #resonator 7 : lower right
                    resLength = baseLength +6*step
                    leftRight = -1
                    meanderSize = 300.
                    nbends = 5
                    if step == 150:
                        turnEdge = 300.
                        newRadius = 50.
                        resonatorStraight1 = resLength -couplingStraight - newRadius * numpy.pi/2 - (nbends-1)*turnEdge - nbends*meanderSize - nbends*2*newRadius * numpy.pi/2
                        start=  (endPtBondPad[0]-(feedlinePinw/2+feedlineGapw+pinw/2)+bondPadPositions[j][0] - bondPadPositions[i][0]-couplingSpacing,endPtBondPad[1]+resonatorSpacing+900)
                        drawTightHanger(mask, start, startDirection, leftRight, pinw, gapw, couplingGap, couplingStraight, resonatorStraight1, meanderSize, turnEdge, resLength, hangerRadius = newRadius)
                    else:
                        turnEdge = 300.
                        resonatorStraight1 = resLength -couplingStraight - radius * numpy.pi/2 - (nbends-1)*turnEdge - nbends*meanderSize - nbends*2*radius * numpy.pi/2
                        start=  (endPtBondPad[0]-(feedlinePinw/2+feedlineGapw+pinw/2)+bondPadPositions[j][0] - bondPadPositions[i][0]-couplingSpacing,endPtBondPad[1]+resonatorSpacing+125)
                        drawTightHanger(mask, start, startDirection, leftRight, pinw, gapw, couplingGap, couplingStraight, resonatorStraight1, meanderSize, turnEdge, resLength, hangerRadius = radius)
                    
                    if step == 150:
                        resLength = baseLength +7*step
                        leftRight = -1
                        meanderSize = 300.
                        turnEdge = 0.
                        nbends = 5
                        newRadius = 50.
                        resonatorStraight1 = resLength -couplingStraight - newRadius * numpy.pi/2 - (nbends-1)*turnEdge - nbends*meanderSize - nbends*2*newRadius * numpy.pi/2
                        if step == 150:
                            start=  (endPtBondPad[0]-(feedlinePinw/2+feedlineGapw+pinw/2)+bondPadPositions[j][0] - bondPadPositions[i][0]-couplingSpacing,endPtBondPad[1]+resonatorSpacing-100)
                        else:
                            start=  (endPtBondPad[0]-(feedlinePinw/2+feedlineGapw+pinw/2)+bondPadPositions[j][0] - bondPadPositions[i][0]-couplingSpacing,endPtBondPad[1]+resonatorSpacing+200-125)
                        drawTightHanger(mask, start, startDirection, leftRight, pinw, gapw, couplingGap, couplingStraight, resonatorStraight1, meanderSize, turnEdge, resLength, hangerRadius = newRadius)
                        
                    return
                
            def draw_ParallelStraighTestChip(step = 150):
                #parallel straight test chip
#                if idx == 3:
#                    step = 300
#                if idx == 4:
#                    step = 250
#                if idx == 5:
#                    step = 200
#                step = 200
#                step = 150
                
                #draw the feedline
                startPtBondPad = bondPadPositions[stopBondPads[-1]]
                mask.last = current_pos
                endPtBondPad = mask.last
                i = stopBondPads[0]
                j = stopBondPads[1]
                CPWStraight(mask,4000,feedlinePinw,feedlineGapw)
                CPWBend(mask,-90,feedlinePinw,feedlineGapw,feedlineRadius)  
                CPWStraight(mask,bondPadPositions[j][0] - bondPadPositions[i][0]-2*feedlineRadius,feedlinePinw,feedlineGapw)
                CPWBend(mask,-90,feedlinePinw,feedlineGapw,feedlineRadius)  
                CPWStraight(mask,4000,feedlinePinw,feedlineGapw)
                mask.last = (mask.last[0], mask.last[1]-(endPtBondPad[1]-startPtBondPad[1]))
                
                
                #setup for the hangers
                leftRight = 1
                scalefactor= 1.
                pinw=scalefactor*pinw_Z0
                gapw=scalefactor*gapw_Z0
                startDirection = -90
                start=  (endPtBondPad[0]+(feedlinePinw/2+feedlineGapw+pinw/2+couplingSpacing),endPtBondPad[1]+resonatorSpacing+1100)
                
                
                baseLength = 7000/2.
                jog = 400
                
                #resonator 1 : lower left
                resLength = baseLength
                leftRight = 1
                meanderSize = 300.
                turnEdge = 0.
                nbends = 5
                resonatorStraight1 = resLength -couplingStraight - radius * numpy.pi/2 - (nbends-1)*turnEdge - nbends*meanderSize - nbends*2*radius * numpy.pi/2
                start=  (endPtBondPad[0]+(feedlinePinw/2+feedlineGapw+pinw/2+couplingSpacing),endPtBondPad[1]+resonatorSpacing+400)
                drawTightHanger(mask, start, startDirection, leftRight, pinw, gapw, couplingGap, couplingStraight, resonatorStraight1, meanderSize, turnEdge, resLength, hangerRadius = radius)
                
                #resonator 2 : middle left
                resLength = baseLength +step
                leftRight = 1
                meanderSize = 300.
                turnEdge = 12.5
                nbends = 5
                resonatorStraight1 = resLength -couplingStraight - radius * numpy.pi/2 - (nbends-1)*turnEdge - nbends*meanderSize - nbends*2*radius * numpy.pi/2
                start=  (endPtBondPad[0]+(feedlinePinw/2+feedlineGapw+pinw/2+couplingSpacing),endPtBondPad[1]+resonatorSpacing+1400)
                drawTightHanger(mask, start, startDirection, leftRight, pinw, gapw, couplingGap, couplingStraight, resonatorStraight1, meanderSize, turnEdge, resLength, hangerRadius = radius)
                
                #resonator 3 : upper left
                resLength = baseLength +2*step
                leftRight = 1
                meanderSize = 300.
                turnEdge = 25.
                nbends = 5
                resonatorStraight1 = resLength -couplingStraight - radius * numpy.pi/2 - (nbends-1)*turnEdge - nbends*meanderSize - nbends*2*radius * numpy.pi/2
                start=  (endPtBondPad[0]+(feedlinePinw/2+feedlineGapw+pinw/2+couplingSpacing),endPtBondPad[1]+resonatorSpacing+2400)
                drawTightHanger(mask, start, startDirection, leftRight, pinw, gapw, couplingGap, couplingStraight, resonatorStraight1, meanderSize, turnEdge, resLength, hangerRadius = radius)
                
                #resonator 4 : upper left
                resLength = baseLength +3*step
                leftRight = 1
                meanderSize = 300.
                turnEdge = 50.
                nbends = 5
                resonatorStraight1 = resLength -couplingStraight - radius * numpy.pi/2 - (nbends-1)*turnEdge - nbends*meanderSize - nbends*2*radius * numpy.pi/2
                start=  (endPtBondPad[0]+(feedlinePinw/2+feedlineGapw+pinw/2+couplingSpacing),endPtBondPad[1]+resonatorSpacing+3400)
                drawTightHanger(mask, start, startDirection, leftRight, pinw, gapw, couplingGap, couplingStraight, resonatorStraight1, meanderSize, turnEdge, resLength, hangerRadius = radius)
                
                
                #resonator 5 : upper right
                resLength = baseLength +4*step
                leftRight = -1
                meanderSize = 300.
                nbends = 5.
                turnEdge = 100.
                resonatorStraight1 = resLength -couplingStraight - radius * numpy.pi/2 - (nbends-1)*turnEdge - nbends*meanderSize - nbends*2*radius * numpy.pi/2
                if step == 150:
                    start=  (endPtBondPad[0]-(feedlinePinw/2+feedlineGapw+pinw/2)+bondPadPositions[j][0] - bondPadPositions[i][0]-couplingSpacing,endPtBondPad[1]+resonatorSpacing+2900)
                else:
                    start=  (endPtBondPad[0]-(feedlinePinw/2+feedlineGapw+pinw/2)+bondPadPositions[j][0] - bondPadPositions[i][0]-couplingSpacing,endPtBondPad[1]+resonatorSpacing+2200-125)
                drawTightHanger(mask, start, startDirection, leftRight, pinw, gapw, couplingGap, couplingStraight, resonatorStraight1, meanderSize, turnEdge, resLength, hangerRadius = radius)
                
                #resonator 6 : middle right
                resLength = baseLength +5*step
                leftRight = -1
                meanderSize = 300.
                turnEdge = 200.
                nbends = 5
                resonatorStraight1 = resLength -couplingStraight - radius * numpy.pi/2 - (nbends-1)*turnEdge - nbends*meanderSize - nbends*2*radius * numpy.pi/2
                if step == 150:
                    start=  (endPtBondPad[0]-(feedlinePinw/2+feedlineGapw+pinw/2)+bondPadPositions[j][0] - bondPadPositions[i][0]-couplingSpacing,endPtBondPad[1]+resonatorSpacing+1900)
                else:
                    start=  (endPtBondPad[0]-(feedlinePinw/2+feedlineGapw+pinw/2)+bondPadPositions[j][0] - bondPadPositions[i][0]-couplingSpacing,endPtBondPad[1]+resonatorSpacing+1200-125)
                drawTightHanger(mask, start, startDirection, leftRight, pinw, gapw, couplingGap, couplingStraight, resonatorStraight1, meanderSize, turnEdge, resLength, hangerRadius = radius)
                
                #resonator 7 : lower right
                resLength = baseLength +6*step
                leftRight = -1
                meanderSize = 300.
                if step == 150:
                    turnEdge = 250.
                else:
                    turnEdge = 300.
                nbends = 5
                resonatorStraight1 = resLength -couplingStraight - radius * numpy.pi/2 - (nbends-1)*turnEdge - nbends*meanderSize - nbends*2*radius * numpy.pi/2
                if step == 150:
                    start=  (endPtBondPad[0]-(feedlinePinw/2+feedlineGapw+pinw/2)+bondPadPositions[j][0] - bondPadPositions[i][0]-couplingSpacing,endPtBondPad[1]+resonatorSpacing+900)
                else:
                    start=  (endPtBondPad[0]-(feedlinePinw/2+feedlineGapw+pinw/2)+bondPadPositions[j][0] - bondPadPositions[i][0]-couplingSpacing,endPtBondPad[1]+resonatorSpacing+125)
                drawTightHanger(mask, start, startDirection, leftRight, pinw, gapw, couplingGap, couplingStraight, resonatorStraight1, meanderSize, turnEdge, resLength, hangerRadius = radius)
                
                if step == 150:
                    resLength = baseLength +7*step
                    leftRight = -1
                    meanderSize = 300.
                    turnEdge = 75.
                    nbends = 5
                    resonatorStraight1 = resLength -couplingStraight - radius * numpy.pi/2 - (nbends-1)*turnEdge - nbends*meanderSize - nbends*2*radius * numpy.pi/2
                    start=  (endPtBondPad[0]-(feedlinePinw/2+feedlineGapw+pinw/2)+bondPadPositions[j][0] - bondPadPositions[i][0]-couplingSpacing,endPtBondPad[1]+resonatorSpacing-100)
                    drawTightHanger(mask, start, startDirection, leftRight, pinw, gapw, couplingGap, couplingStraight, resonatorStraight1, meanderSize, turnEdge, resLength, hangerRadius = radius)
                                    
                return
            
            def draw_BendRadiusTestChip(step = 150):
                #radius test chip
#                step = 150
                
                #draw the feedline
                startPtBondPad = bondPadPositions[stopBondPads[-1]]
                mask.last = current_pos
                endPtBondPad = mask.last
                i = stopBondPads[0]
                j = stopBondPads[1]
                CPWStraight(mask,4000,feedlinePinw,feedlineGapw)
                CPWBend(mask,-90,feedlinePinw,feedlineGapw,feedlineRadius)  
                CPWStraight(mask,bondPadPositions[j][0] - bondPadPositions[i][0]-2*feedlineRadius,feedlinePinw,feedlineGapw)
                CPWBend(mask,-90,feedlinePinw,feedlineGapw,feedlineRadius)  
                CPWStraight(mask,4000,feedlinePinw,feedlineGapw)
                mask.last = (mask.last[0], mask.last[1]-(endPtBondPad[1]-startPtBondPad[1]))
                
                
                #setup for the hangers
                leftRight = 1
                scalefactor= 1.
                pinw=scalefactor*pinw_Z0
                gapw=scalefactor*gapw_Z0
                startDirection = -90
                start=  (endPtBondPad[0]+(feedlinePinw/2+feedlineGapw+pinw/2+couplingSpacing),endPtBondPad[1]+resonatorSpacing+1100)
                
                
                baseLength = 7000/2.
                jog = 400
                
                #non-tight
                resLength = baseLength 
                leftRight = 1
                newRadius = 250.
                if step == 300:
                    resonatorStraight1 = 1100.
                    resonatorStraight2 = 400
                    resonatorStraight3 = jog-200
                else:
                    resonatorStraight1 = 1350.
                    resonatorStraight2 = 400
                    resonatorStraight3 = jog-200
                start=  (endPtBondPad[0]+(feedlinePinw/2+feedlineGapw+pinw/2+couplingSpacing),endPtBondPad[1]+resonatorSpacing+3200)
                drawMixedRadiusHanger(mask, start, startDirection, leftRight, pinw, gapw, couplingGap, couplingStraight, resonatorStraight1, resonatorStraight2, resonatorStraight3, resLength, newRadius)
                
                #non-tight
                resLength = baseLength + step
                leftRight = 1
                newRadius = 200.
                if step == 300:
                    resonatorStraight1 = 1400.
                    resonatorStraight2 = 600
                    resonatorStraight3 = jog-125
                else:
                    resonatorStraight1 = 1400.
                    resonatorStraight2 = 600
                    resonatorStraight3 = jog-125
                start=  (endPtBondPad[0]+(feedlinePinw/2+feedlineGapw+pinw/2+couplingSpacing),endPtBondPad[1]+resonatorSpacing+1750)
                drawMixedRadiusHanger(mask, start, startDirection, leftRight, pinw, gapw, couplingGap, couplingStraight, resonatorStraight1, resonatorStraight2, resonatorStraight3, resLength, newRadius)
                
                #non-tight
                resLength = baseLength + 2*step
                leftRight = 1
                newRadius = 150.
                if step == 300:
                    resonatorStraight1 = 1600.-75
                    resonatorStraight2 = 650
                    resonatorStraight3 = jog + 200
                else:
                    resonatorStraight1 = 1600.-175
                    resonatorStraight2 = 550
                    resonatorStraight3 = jog + 200
                start=  (endPtBondPad[0]+(feedlinePinw/2+feedlineGapw+pinw/2+couplingSpacing),endPtBondPad[1]+resonatorSpacing+250)
                drawMixedRadiusHanger(mask, start, startDirection, leftRight, pinw, gapw, couplingGap, couplingStraight, resonatorStraight1, resonatorStraight2, resonatorStraight3, resLength, newRadius)


                #switching to right
                
                #non-tight
                resLength = baseLength + 3*step
                leftRight = -1
                upDown = -1
                newRadius = 120.
                if step == 300:
                    resonatorStraight1 = 2000.
                    resonatorStraight2 = 1000
                    resonatorStraight3 = jog
                else:
                    resonatorStraight1 = 1700.
                    resonatorStraight2 = 700
                    resonatorStraight3 = jog
                start=  (endPtBondPad[0]-(feedlinePinw/2+feedlineGapw+pinw/2)+bondPadPositions[j][0] - bondPadPositions[i][0]-couplingSpacing,endPtBondPad[1]+resonatorSpacing+3200)
                drawMixedRadiusHanger(mask, start, startDirection, leftRight, pinw, gapw, couplingGap, couplingStraight, resonatorStraight1, resonatorStraight2, resonatorStraight3, resLength, newRadius, upDown)
                
                #non-tight
                resLength = baseLength + 4*step
                leftRight = -1
                upDown = -1
                newRadius = 100.
                if step == 300:
                    resonatorStraight1 = 2000.
                    resonatorStraight2 = 1200
                    resonatorStraight3 = jog
                else:
                    resonatorStraight1 = 1700.
                    resonatorStraight2 = 900
                    resonatorStraight3 = jog
                start=  (endPtBondPad[0]-(feedlinePinw/2+feedlineGapw+pinw/2)+bondPadPositions[j][0] - bondPadPositions[i][0]-couplingSpacing,endPtBondPad[1]+resonatorSpacing+1700)
                drawMixedRadiusHanger(mask, start, startDirection, leftRight, pinw, gapw, couplingGap, couplingStraight, resonatorStraight1, resonatorStraight2, resonatorStraight3, resLength, newRadius, upDown)
                
                #non-tight
                resLength = baseLength + 5*step
                leftRight = -1
                upDown = -1
#                radius = 80.
                newRadius = 80.
                if step == 300:
                    resonatorStraight1 = 2000.
                    resonatorStraight2 = 1200
                    resonatorStraight3 = jog+250
                else:
                    resonatorStraight1 = 1800.
                    resonatorStraight2 = 1000
                    resonatorStraight3 = jog+100
                start=  (endPtBondPad[0]-(feedlinePinw/2+feedlineGapw+pinw/2)+bondPadPositions[j][0] - bondPadPositions[i][0]-couplingSpacing,endPtBondPad[1]+resonatorSpacing+125)
                drawMixedRadiusHanger(mask, start, startDirection, leftRight, pinw, gapw, couplingGap, couplingStraight, resonatorStraight1, resonatorStraight2, resonatorStraight3, resLength, newRadius, upDown)

                #return to original default radius, incase it got messed up
                radius = defaultRadius
                return
            
            def draw_TtestChip(step = 450):
                #t test chip
#                step = 450
                
                #draw the feedline
                startPtBondPad = bondPadPositions[stopBondPads[-1]]
                mask.last = current_pos
                endPtBondPad = mask.last
                i = stopBondPads[0]
                j = stopBondPads[1]
                CPWStraight(mask,4000,feedlinePinw,feedlineGapw)
                CPWBend(mask,-90,feedlinePinw,feedlineGapw,feedlineRadius)  
                CPWStraight(mask,bondPadPositions[j][0] - bondPadPositions[i][0]-2*feedlineRadius,feedlinePinw,feedlineGapw)
                CPWBend(mask,-90,feedlinePinw,feedlineGapw,feedlineRadius)  
                CPWStraight(mask,4000,feedlinePinw,feedlineGapw)
                mask.last = (mask.last[0], mask.last[1]-(endPtBondPad[1]-startPtBondPad[1]))
                
                
                #setup for the hangers
                leftRight = 1
                scalefactor= 1.
                pinw=scalefactor*pinw_Z0
                gapw=scalefactor*gapw_Z0
                startDirection = -90
                start=  (endPtBondPad[0]+(feedlinePinw/2+feedlineGapw+pinw/2+couplingSpacing),endPtBondPad[1]+resonatorSpacing+1100)
                
                
                baseLength = 7000.
                jog = 400
                
                #pair1.1 : lower left
                resLength = baseLength
                leftRight = 1
                upDown = 1
                meanderSize = 500.
                finalStraight = 500.
                turnEdge = 0.
                nbends = 8
                resonatorStraight1 = resLength -couplingStraight - radius * numpy.pi/2 - (nbends-1)*turnEdge - nbends*meanderSize - nbends*2*radius * numpy.pi/2
                start=  (endPtBondPad[0]+(feedlinePinw/2+feedlineGapw+pinw/2+couplingSpacing),endPtBondPad[1]+resonatorSpacing+350)
                drawTTestTightHanger(mask, start, startDirection, leftRight, pinw, gapw, couplingGap, couplingStraight, resonatorStraight1, meanderSize, turnEdge, resLength, hangerRadius = radius, finalStraight = finalStraight, upDown = upDown)
                
                
                #pair 1.2: lower left
                resLength = baseLength 
                leftRight = 1
                upDown = -1
                meanderSize = 500.
                finalStraight = 500.
                turnEdge = 0.
                nbends = 8
                resonatorStraight1 = resLength -couplingStraight - radius * numpy.pi/2 - (nbends-1)*turnEdge - nbends*meanderSize - nbends*2*radius * numpy.pi/2
                start=  (endPtBondPad[0]+(feedlinePinw/2+feedlineGapw+pinw/2+couplingSpacing),endPtBondPad[1]+resonatorSpacing+550)
                drawTTestTightHanger(mask, start, startDirection, leftRight, pinw, gapw, couplingGap, couplingStraight, resonatorStraight1, meanderSize, turnEdge, resLength, hangerRadius = radius, finalStraight = finalStraight, upDown = upDown)
                
                
                #pair2.1 : upper left
                resLength = baseLength + step
                leftRight = 1
                upDown = 1
                meanderSize = 550.
                finalStraight = 500.
                turnEdge = 0.
                nbends = 8
                resonatorStraight1 = resLength -couplingStraight - radius * numpy.pi/2 - (nbends-1)*turnEdge - nbends*meanderSize - nbends*2*radius * numpy.pi/2
                start=  (endPtBondPad[0]+(feedlinePinw/2+feedlineGapw+pinw/2+couplingSpacing),endPtBondPad[1]+resonatorSpacing+2500)
                drawTTestTightHanger(mask, start, startDirection, leftRight, pinw, gapw, couplingGap, couplingStraight, resonatorStraight1, meanderSize, turnEdge, resLength, hangerRadius = radius, finalStraight = finalStraight, upDown = upDown)
                
                
                #pair 2.2: upper left
                resLength = baseLength +step
                leftRight = 1
                upDown = -1
                meanderSize = 550.
                finalStraight = 500.
                turnEdge = 0.
                nbends = 8
                resonatorStraight1 = resLength -couplingStraight - radius * numpy.pi/2 - (nbends-1)*turnEdge - nbends*meanderSize - nbends*2*radius * numpy.pi/2
                start=  (endPtBondPad[0]+(feedlinePinw/2+feedlineGapw+pinw/2+couplingSpacing),endPtBondPad[1]+resonatorSpacing+2700)
                drawTTestTightHanger(mask, start, startDirection, leftRight, pinw, gapw, couplingGap, couplingStraight, resonatorStraight1, meanderSize, turnEdge, resLength, hangerRadius = radius, finalStraight = finalStraight, upDown = upDown)


                #pair3.1 : right
                resLength = baseLength +2*step
                leftRight = -1
                upDown = 1
                meanderSize = 515.
                finalStraight = 500.
                turnEdge = 0.
                nbends = 9
                resonatorStraight1 = resLength -couplingStraight - radius * numpy.pi/2 - (nbends-1)*turnEdge - nbends*meanderSize - nbends*2*radius * numpy.pi/2
                start=  (endPtBondPad[0]-(feedlinePinw/2+feedlineGapw+pinw/2)+bondPadPositions[j][0] - bondPadPositions[i][0]-couplingSpacing,endPtBondPad[1]+resonatorSpacing+500)
                drawTTestTightHanger(mask, start, startDirection, leftRight, pinw, gapw, couplingGap, couplingStraight, resonatorStraight1, meanderSize, turnEdge, resLength, hangerRadius = radius, finalStraight = finalStraight, upDown = upDown)
                
                
                #pair 3.2: right
                resLength = baseLength +2*step
                leftRight = -1
                upDown = -1
                meanderSize = 515.
                finalStraight = 500.
                turnEdge = 0.
                nbends = 9
                resonatorStraight1 = resLength -couplingStraight - radius * numpy.pi/2 - (nbends-1)*turnEdge - nbends*meanderSize - nbends*2*radius * numpy.pi/2
                start=  (endPtBondPad[0]-(feedlinePinw/2+feedlineGapw+pinw/2)+bondPadPositions[j][0] - bondPadPositions[i][0]-couplingSpacing,endPtBondPad[1]+resonatorSpacing+2300)
                drawTTestTightHanger(mask, start, startDirection, leftRight, pinw, gapw, couplingGap, couplingStraight, resonatorStraight1, meanderSize, turnEdge, resLength, hangerRadius = radius, finalStraight = finalStraight, upDown = upDown)

                return
            
            
            
            

            if (idx == 0 or idx == 1):
                draw_BendTestChip(step = 150)

            if (idx == 2 or idx == 3):
#                draw_ParallelStraighTestChip_twoRadii(step = 150)
                draw_ParallelStraighTestChip(step = 150)
                
#            if (idx == 2):
#                draw_ParallelStraighTestChip_twoRadii(step = 150)
#                
#            if (idx == 3):
#                draw_ParallelStraighTestChip(step = 150)
                
#            if (idx == 4 or idx == 5):
#                draw_BendRadiusTestChip(step = 150)
                
            if (idx == 4):
                draw_BendRadiusTestChip(step = 150)
                
            if (idx == 5 or idx == 6 or idx == 7 or idx == 8):
                draw_TtestChip(step = 450)


if __name__=="__main__":
        chip = layout('A')
        mask=sdxf.Drawing()
        mask.blocks.append(chip)
        mask.append(sdxf.Insert(chip.name,point=(0,0)))
        filename = 'hangerResontorTest_tutorial.dxf'
        "Name the output file here"
        mask.saveas(filename)
        
        