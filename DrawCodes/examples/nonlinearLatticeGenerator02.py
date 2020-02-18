# -*- coding: utf-8 -*-
"""
Created on Tue., July 18, 2017

@author: Mattias Fitzpatrick

"""


import sdxf
import random

from MaskMakerPro import *
from math import sin,cos,pi,floor,asin,acos,tan,atan,sqrt

import re
from scipy import *
import pylab
import numpy
import time

from alphanum import alphanum_dict
from random import randrange

import pickle
import datetime
import os
import sys

from matplotlib.pyplot import *

# Code to define lattice info
from EuclideanLayoutGenerator2 import EuclideanLayout

from GeneralLayoutGenerator import *

from LayoutGenerator5 import PlanarLayout

# Helper functions for drawing the graph
from generalLayoutDrawHelpers import *


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
                    (self.topleft_corner[0]+outerEdgeBuffer+innerEdgeBuffer,self.bottomleft_corner[1]+outerEdgeBuffer+innerEdgeBuffer),
                    (self.bottomleft_corner[0]+outerEdgeBuffer+innerEdgeBuffer+chipSpacing,self.bottomleft_corner[1]+outerEdgeBuffer+innerEdgeBuffer),
                    (self.bottomleft_corner[0]+outerEdgeBuffer+innerEdgeBuffer+chipSpacing,self.bottomleft_corner[1]+outerEdgeBuffer+innerEdgeBuffer),
                    (self.bottomleft_corner[0]+outerEdgeBuffer+innerEdgeBuffer+2*chipSpacing,self.bottomleft_corner[1]+outerEdgeBuffer+innerEdgeBuffer),
                    (self.bottomleft_corner[0]+outerEdgeBuffer+innerEdgeBuffer+2*chipSpacing,self.bottomleft_corner[1]+outerEdgeBuffer+innerEdgeBuffer),
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
        scalefactor=7.5/10
        pinwLinear=scalefactor*10.
        gapwLinear=scalefactor*4.186
        x = 3000.
        r = 90.
        num_wiggles = 4
        up = 1
        reslength=7500

        "Cap Parameters"     
        cap_length = 80.0
        cap_gap = 2.0

        scalefactor=1.5
        pinw=scalefactor*10.
        gapw=scalefactor*4.186
        stop_gapw=10*gapw
        stop_pinw=8*pinw
        gapw_buffer = 0.0
        cap_gap_out = 100        
        cap_gap_ext = 0 #This values increases the gap for the outer capacitors
        
        "Bondpad Parameters"        
        bond_pad_length=350. #Length of Rectangular portion of bond pad
        launcher_scalefactor = 38
        launcher_pinw=launcher_scalefactor*10.  
        launcher_gapw=launcher_scalefactor*4.186      
        taper_length= 300.
        launcher_padding = launcher_gapw

        stub_length = 100.
        extension_length = 30.
        chipCornerToBondPad = [1200,1200]

        #####################################################################################################################
        "Find where the lattice coordinates are for Alicia's code"
        #####################################################################################################################
        
        mask=Structure(self,start=0)
        mask.defaults = {'pinw':pinw,'gapw':gapw,'bendradius':r} #Set these as defaults for simplified inherited input

        showLattice=0

        scalefactorSmall=2.0/10
        pinwSmall=scalefactorSmall*10.
        gapwSmall=scalefactorSmall*4.186

        scalefactorLarge=3.0/10
        pinwLarge=scalefactorLarge*10.
        gapwLarge=scalefactorLarge*4.186        

        for idx in range(12):
            if idx == 0:
                side = 1450.
                pinw = pinwSmall
                gapw = gapwSmall
                resonators = numpy.asarray([[0.0, 0.0,side,0.0]])
                layout = GeneralLayout(resonators)
                bondPadConnectionParameters = {'startVertices':[0,1],'resonatorIndices':[0,2],'stopBondPad':[0,1],'directLines':[1,1],'pinwLinear':pinwLinear,'gapwLinear':gapwLinear}
                additionalPositionOffset = [0.,-900.]
                theta = 0.0

            elif idx == 1:
                side = 1450.
                pinw = pinwLarge
                gapw = gapwLarge
                resonators = numpy.asarray([[0.0, 0.0,side,0.0]])
                layout = GeneralLayout(resonators)
                bondPadConnectionParameters = {'startVertices':[1,0],'resonatorIndices':[0,2],'stopBondPad':[2,3],'directLines':[1,1],'pinwLinear':pinwLinear,'gapwLinear':gapwLinear}
                additionalPositionOffset = [0.,1000.]
                theta = 0.0

            elif idx == 2:
                side = 1450.
                pinw = pinwSmall
                gapw = gapwSmall
                resonators = numpy.asarray([[0.0, 0.0,side,0.0],[side,0.0, 2*side,0.0],[2*side,0.0, 3*side,0.0]])
                layout = GeneralLayout(resonators)
                bondPadConnectionParameters = {'startVertices':[0,3],'resonatorIndices':[0,2],'stopBondPad':[0,1],'directLines':[1,1],'pinwLinear':pinwLinear,'gapwLinear':gapwLinear}
                additionalPositionOffset = [0.,-900.]
                theta = 0.0

            elif idx == 3:
                side = 1450.
                pinw = pinwLarge
                gapw = gapwLarge
                resonators = numpy.asarray([[0.0, 0.0,side,0.0],[side,0.0, 2*side,0.0],[2*side,0.0, 3*side,0.0]])
                layout = GeneralLayout(resonators)
                bondPadConnectionParameters = {'startVertices':[0,3],'resonatorIndices':[2,0],'stopBondPad':[3,2],'directLines':[1,1],'pinwLinear':pinwLinear,'gapwLinear':gapwLinear}
                additionalPositionOffset = [0.,1000.]
                theta = 0.0

            elif idx == 4:
                side = 1450.
                pinw = pinwSmall
                gapw = gapwSmall
                resonators = numpy.asarray([[0.0, 0.0,side,0.0],[side,0.0, 2*side,0.0],[2*side,0.0, 3*side,0.0]])
                layout = GeneralLayout(resonators)
                bondPadConnectionParameters = {'startVertices':[0,3],'resonatorIndices':[0,2],'stopBondPad':[0,1],'directLines':[1,1],'pinwLinear':pinwLinear,'gapwLinear':gapwLinear}
                additionalPositionOffset = [0.,-900.]
                theta = 0.0

            elif idx == 5:
                side = 1450.
                pinw = pinwLarge
                gapw = gapwLarge
                resonators = numpy.asarray([[0.0, 0.0,side,0.0],[side,0.0, 2*side,0.0],[2*side,0.0, 3*side,0.0]])
                layout = GeneralLayout(resonators)
                bondPadConnectionParameters = {'startVertices':[0,3],'resonatorIndices':[2,0],'stopBondPad':[3,2],'directLines':[1,1],'pinwLinear':pinwLinear,'gapwLinear':gapwLinear}
                additionalPositionOffset = [0.,1000.]
                theta = 0.0


            elif idx == 6:
                pinw = pinwSmall
                gapw = gapwSmall
                initialLayout = PlanarLayout(gon = 7, vertex = 3, side = 1500, radius_method = 'lin', modeType = 'FW')
                initialLayout.populate(2)
                resonators = initialLayout.get_all_resonators()
                resonators = resonators[0:initialLayout.gon,:]
                layout = GeneralLayout(resonators)
                bondPadConnectionParameters = {'startVertices':[0,1,6],'resonatorIndices':[2,2,2],'stopBondPad':[0,3,1],'directLines':[0,0,0],'pinwLinear':pinwLinear,'gapwLinear':gapwLinear}
                additionalPositionOffset = [0.,0.]
                theta = 0

            elif idx == 7:
                pinw = pinwSmall
                gapw = gapwSmall
                initialLayout = EuclideanLayout(xcells = 1, ycells = 1, lattice_type = 'PeterChain', side = x*sqrt(3), file_path = '', modeType = 'FW')
                layout = GeneralLayout(initialLayout.resonators)
                bondPadConnectionParameters = {'startVertices':[0,4,4],'resonatorIndices':[2,2,0],'stopBondPad':[0,1,2],'directLines':[0,1,0],'pinwLinear':pinwLinear,'gapwLinear':gapwLinear}
                additionalPositionOffset = [-180.,200.]
                theta = -35.

            elif idx == 8:
                pinw = pinwSmall
                gapw = gapwSmall
                initialLayout = PlanarLayout(gon = 7, vertex = 3, side = 1450, radius_method = 'lin', modeType = 'FW')
                initialLayout.populate(2)
                resonators = initialLayout.get_all_resonators()
                resonators = resonators[0:initialLayout.gon,:]
                layout = GeneralLayout(resonators)
                bondPadConnectionParameters = {'startVertices':[0,1,6],'resonatorIndices':[2,2,2],'stopBondPad':[0,3,1],'directLines':[0,0,0],'pinwLinear':pinwLinear,'gapwLinear':gapwLinear}
                additionalPositionOffset = [0.,0.]
                theta = 0

            elif idx == 9:
                pinw = pinwLarge
                gapw = gapwLarge
                initialLayout = PlanarLayout(gon = 7, vertex = 3, side = 1450, radius_method = 'lin', modeType = 'FW')
                initialLayout.populate(2)
                resonators = initialLayout.get_all_resonators()
                resonators = resonators[0:initialLayout.gon,:]
                layout = GeneralLayout(resonators)
                bondPadConnectionParameters = {'startVertices':[0,1,6],'resonatorIndices':[2,2,2],'stopBondPad':[0,3,1],'directLines':[0,0,0],'pinwLinear':pinwLinear,'gapwLinear':gapwLinear}
                additionalPositionOffset = [0.,0.]
                theta = 0

            elif idx == 10:
                pinw = pinwLarge
                gapw = gapwLarge
                initialLayout = EuclideanLayout(xcells = 1, ycells = 1, lattice_type = 'PeterChain', side = x*sqrt(3), file_path = '', modeType = 'FW')
                layout = GeneralLayout(initialLayout.resonators)
                bondPadConnectionParameters = {'startVertices':[0,4,4],'resonatorIndices':[2,2,0],'stopBondPad':[0,1,2],'directLines':[0,1,0],'pinwLinear':pinwLinear,'gapwLinear':gapwLinear}
                additionalPositionOffset = [-200.,200.]
                theta = -35.

            elif idx == 11:
                pinw = pinwLarge
                gapw = gapwLarge
                initialLayout = PlanarLayout(gon = 7, vertex = 3, side = 1450, radius_method = 'lin', modeType = 'FW')
                initialLayout.populate(2)
                resonators = initialLayout.get_all_resonators()
                resonators = resonators[0:initialLayout.gon,:]
                layout = GeneralLayout(resonators)
                bondPadConnectionParameters = {'startVertices':[0,1,6],'resonatorIndices':[2,2,2],'stopBondPad':[0,3,1],'directLines':[0,0,0],'pinwLinear':pinwLinear,'gapwLinear':gapwLinear}
                additionalPositionOffset = [0.,0.]
                theta = 0


            # construct a general layout
            layout.resonators = rotate_resonators(layout.resonators, theta)    
            
            # find the mean position of the resonators
            meanPos = [0.,0.]
            for res in layout.resonators:
                meanPos = [meanPos[0]+res[0]+res[2],meanPos[1]+res[1]+res[3]]
            meanPos = [meanPos[0]/(2*len(layout.resonators)),meanPos[1]/(2*len(layout.resonators))]

            # print('layout.resonators',layout.resonators)

            resonatorShift = chipStarts[idx]

            layout.resonators = shift_resonators(layout.resonators, resonatorShift[0]+smallChipSize/2-meanPos[0]+additionalPositionOffset[0], resonatorShift[1]+smallChipSize/2-meanPos[1]+additionalPositionOffset[1])

            layout = GeneralLayout(layout.resonators)

            # get the Hamiltonian
            layout.generate_Hamiltonian()
            H = layout.H

            # find the vertex dictionarys
            layout.generate_vertex_dict()
            vertexDict = layout.vertexDict
            vertexPositions = layout.coords


            if showLattice == 1:
                pylab.figure(4)
                pylab.clf()
                ax = pylab.subplot(1,1,1)
                layout.draw_resonator_lattice(ax)

                idx = 0
                for vtexPos in vertexPositions:
                    ax.text(vtexPos[0], vtexPos[1], str(idx), fontsize=6,color='red')
                    idx +=1
                ax.set_aspect('equal')
                # pylab.show()    

                pylab.figure(5)
                ax = pylab.subplot(1,1,1)
                ax.imshow(layout.H)
                pylab.show()



            # extract the resonators
            resonators = layout.resonators

            drawDictionaries = getDrawDictionaries(vertexDict,vertexPositions,resonators)
                
            startVertices = bondPadConnectionParameters['startVertices']
            stopBondPads = bondPadConnectionParameters['stopBondPad']
            directLines = bondPadConnectionParameters['directLines']
            startResonators = bondPadConnectionParameters['resonatorIndices']

            bufferDistance = 1;
            # print('vertexDict',vertexDict,'vertexPositions',vertexPositions,'resonators',resonators,
            #     'verboseVertexDict',drawDictionaries.verboseVertexDict,'vertexAngles',drawDictionaries.vertexAngles,
            #     'resonatorInputOuputAngles',drawDictionaries.resonatorInputOuputAngles,'capacitorAngles',drawDictionaries.capacitorAngles,
            #     'vertexDict',drawDictionaries.vertexDict)    

            drawThreewayCaps(mask,drawDictionaries,cap_length,cap_gap,gapw_buffer,stub_length,pinw,stop_pinw,gapw,stop_gapw,bondPadConnectionParameters,bufferDistance)

            drawResonators(mask,drawDictionaries,reslength,x,r,cap_length,pinw,gapw,num_wiggles,up)

            # define potential bond pad positions
            lowerLeftChipCorner = chipStarts[idx]
            bondPadPositions = [[lowerLeftChipCorner[0]+chipCornerToBondPad[0],lowerLeftChipCorner[1]+chipCornerToBondPad[1]],
                                [lowerLeftChipCorner[0]+smallChipSize-chipCornerToBondPad[0],lowerLeftChipCorner[1]+chipCornerToBondPad[1]],
                                [lowerLeftChipCorner[0]+smallChipSize-chipCornerToBondPad[0],lowerLeftChipCorner[1]+smallChipSize-chipCornerToBondPad[1]],
                                [lowerLeftChipCorner[0]+chipCornerToBondPad[0],lowerLeftChipCorner[1]+smallChipSize-chipCornerToBondPad[1]]]
            bondPadAngles = [90,90,270,270]

            for idx in range(len(bondPadPositions)):

                # find where that bond pad is
                bondPadPos = bondPadPositions[idx]
                bondPadAng = bondPadAngles[idx]
                bondPadEndPos = [bondPadPos[0]+(bond_pad_length+launcher_padding+taper_length)*numpy.cos(bondPadAng*numpy.pi/180.),bondPadPos[1]+(bond_pad_length+launcher_padding+taper_length)*numpy.sin(bondPadAng*numpy.pi/180.)]

                # figure out if that bond pad is connected
                if idx in stopBondPads:

                    
                    drawBondPad(mask,bondPadPos,bondPadAng,pinwLinear,gapwLinear,bond_pad_length,launcher_pinw,launcher_gapw,taper_length,launcher_padding)

                    index = stopBondPads.index(idx)
                    startVertex = startVertices[index]
                    startRes = startResonators[index]
                    directLine = directLines[index]
                    bondPadConnection(mask,drawDictionaries,extension_length,lowerLeftChipCorner,smallChipSize,edgeBufferDistance,directLine,r,bondPadEndPos,bondPadAng,startVertex,startRes,cap_length,pinwLinear,gapwLinear)
                    # print('idx',idx,'stopBondPads',stopBondPads[index],'startVertex',startVertex)


if __name__=="__main__":
        chip = layout('A')

        mask=sdxf.Drawing()
        mask.blocks.append(chip)
    
        mask.append(sdxf.Insert(chip.name,point=(0,0)))
   
        "Name the output file here"
        mask.saveas('smallNonlinearLattices02_test_tutorial.dxf')
        