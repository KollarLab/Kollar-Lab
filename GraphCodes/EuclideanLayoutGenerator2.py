#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 17 13:22:52 2018

@author: kollar2

modified from LayoutGnerator5 which makes hyperbolic lattices
Tried to keep as muchof the structure and syntax consistent.
Many many things are the same, but some things had to change because this one 
is build unit cell by unit sell and not itteration shell by itteration shell

v0 - basic Huse lattice only

v2 - adding kagome and composite Huse/kagome. UnitCell class will get a lot fo 
    new functions, including Bloch theory calculations
    
    Also added alternate versions on the Huse idea trying to insert other types of polygons.
    Tried 7,4 and 12,3 and 8,4.
    Composite cell function can handle original 7,5 Huse and the new variants.
    (Note: 7,4 is not 3-regular, but it is still triangle protected.)
    
    7-25-18 Added zorder as optional argument to all the plot functions
    
    8-14-18 adding funtions that allow new unit cells to be made either by subdividing an old cell
    or taking its line graph
    

UnitCell Class
    Object to conveniently hold and define a single unit cell. Will store the number
    of site, where they are, what the links are between them and neighboring unit cells,
    and which sites are needed to close an incomplete unit cell
    
    Supported Types:
        Huse (v0)
    
    Methods:
        ########
        #generating the cell
        ########
        _generate_kagome_cell
        _generate_Huse_cell
        _generate_PeterChain_cell
        _generate_PeterChain2_cell
        _generate_square_cell
        _generate_84Huse_cell
        _generate_74Huse_cell
        _generate_123Huse_cell
        _generate_Hk_composite_cell
        _generate_arbitrary_Cell
        
        ########
        #drawing the cell
        ########
        draw_resonators
        draw_resonator_end_points
        draw_sites
        draw_SDlinks
        _get_orientation_plot_points
        draw_site_orientations
        
        ########
        #auto construction functions for SD links
        ########
        _auto_generate_SDlinks
        _auto_generate_cell_SDlinks
        
        ########
        #Bloch theory function
        ########
        generate_Bloch_matrix
        compute_band_structure
        plot_band_cut
        plot_bloch_wave
        plot_bloch_wave_end_state
        
        ########
        #making new cells
        ########
        split_cell
        line_graph_cell #for now this only works for coordination numbers 3 or smaller
        #4 and up require more link matrix space to be allocated.
     
    Sample syntax:
        #####
        #creating unit cell
        #####
        from EuclideanLayoutGenerator import UnitCell
        #built-in cell
        testCell = UnitCell(lattice_type = 'Huse', side = 1)
        
        #custom cell
        testCell = UnitCell(lattice_type = 'name', side = 1, resonators = resonatorMat, a1 = vec1, a2 = vec2)
        


EuclideanLayout Class
    Chose your UnitCell type, wave type, and number of unit cells and make a lattice
     
     Methods:
        ###########
        #automated construction, saving, loading
        ##########
        populate (autoruns at creation time)
        save
        load
        
        ########
        #functions to generate the resonator lattice
        #######
        generateLattice
        _fix_edge_resonators (already stores some SD properties of fixed edge)
         
        #######
        #resonator lattice get /view functions
        #######
        get_xs
        get_ys
        draw_resonator_lattice
        draw_resonator_end_points
        get_all_resonators
        get_coords
        
        ########
        #functions to generate effective JC-Hubbard lattice (semiduals)
        ######## 
        generate_semiduals
        _fix_SDedge
        
        #######
        #get and view functions for the JC-Hubbard (semi-dual lattice)
        #######
        draw_SD_points
        draw_SDlinks
        get_semidual_points (semi-defunct)
        
        ######
        #Hamiltonian related methods
        ######
        generate_Hamiltonian
        get_eigs
        
        ##########
        #methods for calculating/looking at states and interactions
        #########
        get_SDindex (removed for now. Needs to be reimplemented in sensible fashion)
        build_local_state
        V_int
        V_int_map
        plot_layout_state
        plot_map_state
        get_end_state_plot_points
        plot_end_layout_state
        
        
    Sample syntax:
        #####
        #loading precalculated layout
        #####
        from EuclideanLayoutGenerator import EuclideanLayout
        testLattice = EuclideanLayout(file_path = 'Huse_4x4_FW.pkl')
        
        #####
        #making new layout
        #####
        from EuclideanLayoutGenerator import EuclideanLayout
        #from built-in cell
        testLattice = EuclideanLayout(xcells = 4, ycells = 4, lattice_type = 'Huse', side = 1, file_path = '', modeType = 'FW')
        
        #from custom cell
        testCell = UnitCell(lattice_type = 'name', side = 1, resonators = resonatorMat, a1 = vec1, a2 = vec2)
        testLattice = EuclideanLayout(xcells = 4, ycells = 4, modeType = 'FW', resonatorsOnly=False, initialCell = testCell)
        
        #####
        #saving computed layout
        #####
        testLattice.save( name = 'filename.pkl') #filename can be a full path, but must have .pkl extension
        



"""



import re
import scipy
import pylab
import numpy
import time

import pickle
import datetime
import os
import sys

import scipy.linalg


class UnitCell(object):
    def __init__(self, lattice_type, side = 1, resonators = '', a1 = [1,0], a2 = [0,1]):
        '''
        optional resonator and a1, a2 reciprocal lattice vector input arguments will only be used 
        if making a cell of non-built-in type using _generate_arbitrary_cell
        '''
        
        self.side = side*1.0
        
    #auto parse variants on Huse-type lattices
        match = re.search(r'(\d*)(Huse)(\d*)(\_*)(\d*)', lattice_type)
        if match:
            #Huse type lattice of some sort
            HuseType = match.groups()[0] + match.groups()[1]
            
            if match.groups()[3]  =='':
                #regular unit cell
                self.type = lattice_type
                generateMethod = getattr(self, '_generate_' + HuseType + '_cell')
                generateMethod(side= self.side)
            else:
                #composite unit cell
                defect_type = HuseType
                xtrans = int(match.groups()[2])
                ytrans = int(match.groups()[4])
                
                self.type = lattice_type
                self._generate_Hk_composite_cell(xtrans, ytrans, side = self.side, defect_type = defect_type)
        
        elif lattice_type == 'kagome':
            self.type = lattice_type
            self._generate_kagome_cell(self.side)
            
        elif lattice_type == 'PeterChain':
            self.type = lattice_type
            self._generate_PeterChain_cell(self.side)
            
        elif lattice_type == 'PeterChain_tail':
            self.type = lattice_type
            self._generate_PeterChain2_cell(self.side)
            
        elif lattice_type == 'square':
            self.type = lattice_type
            self._generate_square_cell(self.side)  
            
        else:
            #arbitrary lattice type
            self.type = lattice_type
            self._generate_arbitrary_cell(resonators, a1, a2)
            
            
    ########
    #generator functions for unit cells
    ########        
    def _generate_kagome_cell(self, side = 1):
        '''
        generate kagome-type unit cell
        '''
        #set up the sites
        self.numSites = 3
        xs = numpy.zeros(self.numSites)
        ys = numpy.zeros(self.numSites)
        
        #set up the lattice vectors
        self.a1 = numpy.asarray([self.side*numpy.sqrt(3)/2, self.side/2])
        self.a2 = numpy.asarray([0, self.side])
        dy = self.a1[1]/2
        dx = self.a1[0]/2
        xcorr = self.side/numpy.sqrt(3)/2/2
        
        #set up the positions of the sites of the effective lattice. ! look to newer functions for auto way to do these
        xs = numpy.asarray([-dx, -dx, 0])
        ys = numpy.asarray([dy, -dy, -2*dy])
        self.SDx = xs
        self.SDy = ys
        
        #set up the poisitions of all the resonators  and their end points
        self.resonators = numpy.zeros((self.numSites,4)) #pairs of resonator end points for each resonator
        self.coords = numpy.zeros((self.numSites,2)) #set of all resonator start points
        
        a = self.side/numpy.sqrt(3)
        b = self.a1[0]-a
        #xo,yo,x1,y1
        #define them so their orientation matches the chosen one. First entry is plus end, second is minus
        self.resonators[0,:] = [-a/2, 2*dy, -b-a/2,  0]
        self.resonators[1,:] = [-a/2-b, 0, -a/2,  -2*dy]
        self.resonators[2,:] = [a/2, -2*dy, -a/2,  -2*dy]
        
        self.coords = self.get_coords(self.resonators)
        
        
        #####manual population of the SD links
#        #matrix to hold all the bonds
#        #starting site, ending site, number units cells over in a1, number unit cells over in a2, initial end type, final end type
#        oldlinks = numpy.zeros((self.numSites*4, 4)) #without orientation
#        links = numpy.zeros((self.numSites*4, 6))   #with orientation
#        #orientation defines by +x or +y is the +end 
#        
#        #fill in the links
#        links[0,:] = [0,2,0,1,   1,0]
#        links[1,:] = [0,1,0,0,   0,1]
#        links[2,:] = [0,2,-1,1,  0,1]
#        links[3,:] = [0,1,0,1,   1,0]
#        
#        links[4,:] = [1,2,0,0,   0,0]
#        links[5,:] = [1,0,0,0,   1,0]
#        links[6,:] = [1,2,-1,1,  1,1]
#        links[7,:] = [1,0,0,-1,  0,1]
#        
#        links[8,:] = [2,1,0,0,   0,0]
#        links[9,:] = [2,0,0,-1,  0,1]
#        links[10,:] =[2,0,1,-1,  1,0]
#        links[11,:] =[2,1,1,-1,  1,1]
#        
#        oldlinks = links[:,0:4]
#        self.SDlinks = oldlinks
#        self.SDHWlinks = links
        
        #####auto population of the SD links
        self._auto_generate_SDlinks()
        
        
        
        #make note of which resonator you need in order to close the unit cell
        closure = {}
        
        #a1 direction (x)
        closure[(1,0)] =numpy.asarray([1])
        #-a1 direction (-x)
        closure[(-1,0)] =numpy.asarray([])
        
        #a2 direction (y)
        closure[(0,1)] =numpy.asarray([2])
        #-a2 direction (y)
        closure[(0,-1)] =numpy.asarray([])
        
         #a1,a2 direction (x,y)
        closure[(1,1)] =numpy.asarray([])
        #-a1,a2 direction (-x,y)
        closure[(-1,1)] =numpy.asarray([])
        #a1,-a2 direction (x,-y)
        closure[(1,-1)] =numpy.asarray([0])
        #-a1,-a2 direction (-x,-y)
        closure[(-1,-1)] =numpy.asarray([])
        self.closure = closure
        
        
        return
    
    def _generate_Huse_cell(self, side = 1):        
        '''
        generate standard Huse-type unit cell (7,5)
        '''
        
        #set up the sites
        self.numSites = 12
        xs = numpy.zeros(self.numSites)
        ys = numpy.zeros(self.numSites)
        
        #set up the lattice vectors
        self.a1 = numpy.asarray([self.side*numpy.sqrt(3)/2, self.side/2])
        self.a2 = numpy.asarray([0, 3*self.side])
        dy = self.a1[1]/2
        dx = self.a1[0]/2
        xcorr = self.side/numpy.sqrt(3)/2/2
        
        
        #set up the positions of the sites of the effective lattice
        xs = numpy.asarray([-dx, -dx,-dx,-dx,-dx,-dx, 0, -xcorr, +xcorr, 0, -xcorr, +xcorr])
        ys = numpy.asarray([5*dy, 3*dy, dy, -dy, -3*dy, -5*dy, -6*dy, 2*dy, 2*dy, 0, -2*dy, -2*dy])
        self.SDx = xs
        self.SDy = ys
        
        #set up the poisitions of all the resonators  and their end points
        self.resonators = numpy.zeros((self.numSites,4)) #pairs of resonator end points for each resonator
        self.coords = numpy.zeros((self.numSites,2)) #set of all resoantor start points
        
        a = self.side/numpy.sqrt(3)
        b = self.a1[0]-a
        #xo,yo,x1,y1
        #define them so their orientation matches the chosen one. First entry is plus end, second is minus
        self.resonators[0,:] = [-a/2, 6*dy, -b-a/2,  4*dy]
        self.resonators[1,:] = [-a/2-b, 4*dy, -a/2,  2*dy]
        self.resonators[2,:] = [-a/2, 2*dy, -b-a/2,  0]
        self.resonators[3,:] = [-a/2-b, 0, -a/2,  -2*dy]
        self.resonators[4,:] = [-a/2, -2*dy, -b-a/2,  -4*dy]
        self.resonators[5,:] = [-a/2-b, -4*dy, -a/2,  -6*dy]
        self.resonators[6,:] = [a/2, -6*dy, -a/2,  -6*dy]
        self.resonators[7,:] = [0, 2*dy, -a/2,  2*dy]
        self.resonators[8,:] = [a/2, 2*dy, 0,  2*dy]
        self.resonators[9,:] = [0, 2*dy, 0,  -2*dy]
        self.resonators[10,:] = [0, -2*dy, -a/2,  -2*dy]
        self.resonators[11,:] = [a/2, -2*dy, 0,  -2*dy]
        
        self.coords = self.get_coords(self.resonators)
        
        
        
        ######manual population of the SD links
#        #matrix to hold all the bonds
#        #starting site, ending site, number units cells over in a1, number unit cells over in a2, initial end type, final end type
#        oldlinks = numpy.zeros((self.numSites*4, 4)) #without orientation
#        links = numpy.zeros((self.numSites*4, 6))   #with orientation
#        #orientation defines by +x or +y is the +end 
#        
#        #fill in the links
#        links[0,:] = [0 , 6, 0,1,  1,0]
#        links[1,:] = [0 , 6, -1,1, 0,1]
#        links[2,:] = [0 , 5, 0,1,  1,0]
#        links[3,:] = [0 , 1, 0,0,  0,1]
#        
#        links[4,:] = [1,2, 0,0,    0,1]
#        links[5,:] = [1,0, 0,0,    1,0]
#        links[6,:] = [1,7, 0,0,    0,0]
#        links[7,:] = [1,6, -1,1,   1,1]
#        
#        links[8,:] = [2,1, 0,0,    1,0]
#        links[9,:] = [2,7, 0,0,    1,0]
#        links[10,:] = [2,3, 0,0,   0,1]
#        links[11,:] = [2,8, -1,0,  0,1]
#        
#        links[12,:] = [3,2, 0,0,   1,0]
#        links[13,:] = [3,4, 0,0,   0,1]
#        links[14,:] = [3,10, 0,0,  0,0]
#        links[15,:] = [3,8, -1,0,  1,1]
#        
#        links[16,:] = [4,3, 0,0,   1,0]
#        links[17,:] = [4,10, 0,0,  1,0]
#        links[18,:] = [4,5, 0,0,   0,1]
#        links[19,:] = [4,11,-1,0,  0,1]
#        
#        links[20,:] = [5,4, 0,0,   1,0]
#        links[21,:] = [5,6, 0,0,   0,0]
#        links[22,:] = [5,11,-1,0,  1,1]
#        links[23,:] = [5,0, 0,-1,  0,1]
#        
#        links[24,:] = [6,5, 0,0,   0,0]
#        links[25,:] = [6,0, 0,-1,  0,1]
#        links[26,:] = [6,1, 1,-1,  1,1]
#        links[27,:] = [6,0, 1,-1,  1,0]
#        
#        links[28,:] = [7,1, 0,0,   0,0]
#        links[29,:] = [7,2, 0,0,   0,1]
#        links[30,:] = [7,9, 0,0,   1,1]
#        links[31,:] = [7,8, 0,0,   1,0]
#        
#        links[32,:] = [8,7, 0,0,   0,1]
#        links[33,:] = [8,9, 0,0,   0,1]
#        links[34,:] = [8,2, 1,0,   1,0]
#        links[35,:] = [8,3, 1,0,   1,1]
#        
#        links[36,:] = [9,7, 0,0,   1,1]
#        links[37,:] = [9,8, 0,0,   1,0]
#        links[38,:] = [9,10, 0,0,  0,1]
#        links[39,:] = [9,11, 0,0,  0,0]
#        
#        links[40,:] = [10,9, 0,0,  1,0]
#        links[41,:] = [10,11,0,0,  1,0]
#        links[42,:] = [10,3, 0,0,  0,0]
#        links[43,:] = [10,4, 0,0,  0,1]
#        
#        links[44,:] = [11,9, 0,0,  0,0]
#        links[45,:] = [11,10,0,0,  0,1]
#        links[46,:] = [11,4, 1,0,  1,0]
#        links[47,:] = [11,5, 1,0,  1,1]
#        
#        oldlinks = links[:,0:4]
#        self.SDlinks = oldlinks
#        self.SDHWlinks = links
        
        #####auto population of the SD links
        self._auto_generate_SDlinks()
        
        
        #make note of which resonator you need in order to close the unit cell
        closure = {}
        
        #a1 direction (x)
        closure[(1,0)] =numpy.asarray([1,2,3,4,5])
        #-a1 direction (-x)
        closure[(-1,0)] =numpy.asarray([])
        
        #a2 direction (y)
        closure[(0,1)] =numpy.asarray([6])
        #-a2 direction (y)
        closure[(0,-1)] =numpy.asarray([])
        
         #a1,a2 direction (x,y)
        closure[(1,1)] =numpy.asarray([])
        #-a1,a2 direction (-x,y)
        closure[(-1,1)] =numpy.asarray([])
        #a1,-a2 direction (x,-y)
        closure[(1,-1)] =numpy.asarray([0])
        #-a1,-a2 direction (-x,-y)
        closure[(-1,-1)] =numpy.asarray([])
        self.closure = closure
        
        return
    
    def _generate_PeterChain_cell(self, side = 1):
        '''
        generate Pater-Sarnak chain unit cell
        '''
        #set up the sites
        self.numSites = 6
        xs = numpy.zeros(self.numSites)
        ys = numpy.zeros(self.numSites)
        
        #set up the lattice vectors
        self.a1 = numpy.asarray([self.side, 0])
        self.a2 = numpy.asarray([0, 2*self.side])
        dy = self.a1[1]/2
        dx = self.a1[0]/2
        xcorr = self.side/numpy.sqrt(3)/2/2
        
        
        #set up the poisitions of all the resonators  and their end points
        self.resonators = numpy.zeros((self.numSites,4)) #pairs of resonator end points for each resonator
        self.coords = numpy.zeros((self.numSites,2)) #set of all resonator start points
        
        a = self.side/(2*numpy.sqrt(2) + 1)
        b = numpy.sqrt(2)*a
        #xo,yo,x1,y1
        #define them so their orientation matches the chosen one. First entry is plus end, second is minus
        self.resonators[0,:] = [-a-b, 0, -b,  0]
        self.resonators[1,:] = [-b, 0, 0,  b]
        self.resonators[2,:] = [0, b, b,  0]
        self.resonators[3,:] = [-b, 0, 0,  -b]
        self.resonators[4,:] = [0, -b, b,  0]
        self.resonators[5,:] = [0, -b, 0,  b]
        
        self.coords = self.get_coords(self.resonators)
        
        #set up the positions of the sites of the effective lattice
        xs = numpy.zeros(self.numSites)
        ys = numpy.zeros(self.numSites)
        for rind in range(0, self.resonators.shape[0]):
            res = self.resonators[rind,:]
            xs[rind] = (res[0] + res[2])/2
            ys[rind] = (res[1] + res[3])/2
        self.SDx = xs
        self.SDy = ys
        
        
        #####auto population of the SD links
        self._auto_generate_SDlinks()
        
        
        
        #make note of which resonator you need in order to close the unit cell
        closure = {}
        
        #a1 direction (x)
        closure[(1,0)] =numpy.asarray([])
        #-a1 direction (-x)
        closure[(-1,0)] =numpy.asarray([])
        
        #a2 direction (y)
        closure[(0,1)] =numpy.asarray([])
        #-a2 direction (y)
        closure[(0,-1)] =numpy.asarray([])
        
         #a1,a2 direction (x,y)
        closure[(1,1)] =numpy.asarray([])
        #-a1,a2 direction (-x,y)
        closure[(-1,1)] =numpy.asarray([])
        #a1,-a2 direction (x,-y)
        closure[(1,-1)] =numpy.asarray([])
        #-a1,-a2 direction (-x,-y)
        closure[(-1,-1)] =numpy.asarray([])
        self.closure = closure
        
        
        return
    
    def _generate_PeterChain2_cell(self, side = 1):
        '''
        generate Pater-Sarnak chain unit cell
        '''
        #set up the sites
        self.numSites = 7
        xs = numpy.zeros(self.numSites)
        ys = numpy.zeros(self.numSites)
        
        #set up the lattice vectors
        self.a1 = numpy.asarray([self.side, 0])
        self.a2 = numpy.asarray([0, 2*self.side])
        dy = self.a1[1]/2
        dx = self.a1[0]/2
        xcorr = self.side/numpy.sqrt(3)/2/2
        
        
        #set up the poisitions of all the resonators  and their end points
        self.resonators = numpy.zeros((self.numSites,4)) #pairs of resonator end points for each resonator
        self.coords = numpy.zeros((self.numSites,2)) #set of all resonator start points
        
        a = self.side/(2*numpy.sqrt(2) + 2)
        b = numpy.sqrt(2)*a
        #xo,yo,x1,y1
        #define them so their orientation matches the chosen one. First entry is plus end, second is minus
        self.resonators[0,:] = [-a-b, 0, -b,  0]
        self.resonators[1,:] = [-b, 0, 0,  b]
        self.resonators[2,:] = [0, b, b,  0]
        self.resonators[3,:] = [-b, 0, 0,  -b]
        self.resonators[4,:] = [0, -b, b,  0]
        self.resonators[5,:] = [0, -b, 0,  b]
        self.resonators[6,:] = [b, 0, a+b, 0]
        
        self.coords = self.get_coords(self.resonators)
        
        #set up the positions of the sites of the effective lattice
        xs = numpy.zeros(self.numSites)
        ys = numpy.zeros(self.numSites)
        for rind in range(0, self.resonators.shape[0]):
            res = self.resonators[rind,:]
            xs[rind] = (res[0] + res[2])/2
            ys[rind] = (res[1] + res[3])/2
        self.SDx = xs
        self.SDy = ys
        
        
        #####auto population of the SD links
        self._auto_generate_SDlinks()
        
        
        
        #make note of which resonator you need in order to close the unit cell
        closure = {}
        
        #a1 direction (x)
        closure[(1,0)] =numpy.asarray([])
        #-a1 direction (-x)
        closure[(-1,0)] =numpy.asarray([])
        
        #a2 direction (y)
        closure[(0,1)] =numpy.asarray([])
        #-a2 direction (y)
        closure[(0,-1)] =numpy.asarray([])
        
         #a1,a2 direction (x,y)
        closure[(1,1)] =numpy.asarray([])
        #-a1,a2 direction (-x,y)
        closure[(-1,1)] =numpy.asarray([])
        #a1,-a2 direction (x,-y)
        closure[(1,-1)] =numpy.asarray([])
        #-a1,-a2 direction (-x,-y)
        closure[(-1,-1)] =numpy.asarray([])
        self.closure = closure
        
        
        return
    
    def _generate_square_cell(self, side = 1):
        '''
        generate sqare lattice unit cell
        '''
        #set up the sites
        self.numSites = 2
#        self.numSites = 4
        xs = numpy.zeros(self.numSites)
        ys = numpy.zeros(self.numSites)
        
        #set up the lattice vectors
        self.a1 = numpy.asarray([self.side, 0])
        self.a2 = numpy.asarray([0, self.side])
        dy = self.a1[1]/2
        dx = self.a1[0]/2
        xcorr = self.side/numpy.sqrt(3)/2/2
        
        
        #set up the poisitions of all the resonators  and their end points
        self.resonators = numpy.zeros((self.numSites,4)) #pairs of resonator end points for each resonator
        self.coords = numpy.zeros((self.numSites,2)) #set of all resonator start points
        
        a = self.side
        #xo,yo,x1,y1
        #define them so their orientation matches the chosen one. First entry is plus end, second is minus
        self.resonators[0,:] = [0, 0, 0, a]
        self.resonators[1,:] = [a, 0, 0, 0]
        
#        self.resonators[0,:] = [0, 0, 0, a/2.]
#        self.resonators[1,:] = [0, 0, a/2., 0]
#        self.resonators[2,:] = [-a/2., 0, 0, 0]
#        self.resonators[3,:] = [0, -a/2., 0, 0]
        
        self.coords = self.get_coords(self.resonators)
        
        #set up the positions of the sites of the effective lattice
        xs = numpy.zeros(self.numSites)
        ys = numpy.zeros(self.numSites)
        for rind in range(0, self.resonators.shape[0]):
            res = self.resonators[rind,:]
            xs[rind] = (res[0] + res[2])/2
            ys[rind] = (res[1] + res[3])/2
        self.SDx = xs
        self.SDy = ys
        
        
        #####auto population of the SD links
        self._auto_generate_SDlinks()
        
        
        
        #make note of which resonator you need in order to close the unit cell
        closure = {}
        
        #a1 direction (x)
        closure[(1,0)] =numpy.asarray([])
        #-a1 direction (-x)
        closure[(-1,0)] =numpy.asarray([])
        
        #a2 direction (y)
        closure[(0,1)] =numpy.asarray([])
        #-a2 direction (y)
        closure[(0,-1)] =numpy.asarray([])
        
         #a1,a2 direction (x,y)
        closure[(1,1)] =numpy.asarray([])
        #-a1,a2 direction (-x,y)
        closure[(-1,1)] =numpy.asarray([])
        #a1,-a2 direction (x,-y)
        closure[(1,-1)] =numpy.asarray([])
        #-a1,-a2 direction (-x,-y)
        closure[(-1,-1)] =numpy.asarray([])
        self.closure = closure
        
        
        return
    
    def _generate_84Huse_cell(self, side = 1):        
        '''
        generate 8,4 Huse variant unit cell
        '''
        #set up the sites
        self.numSites = 12
        xs = numpy.zeros(self.numSites)
        ys = numpy.zeros(self.numSites)
        
        #set up the lattice vectors
        self.a1 = numpy.asarray([self.side*numpy.sqrt(3)/2, self.side/2])
        self.a2 = numpy.asarray([0, 3*self.side])
        dy = self.a1[1]/2
        
        #set up the poisitions of all the resonators  and their end points
        self.resonators = numpy.zeros((self.numSites,4)) #pairs of resonator end points for each resonator
        self.coords = numpy.zeros((self.numSites,2)) #set of all resoantor start points
        
        a = self.side/numpy.sqrt(3)
        b = self.a1[0]-a
        #xo,yo,x1,y1
        #define them so their orientation matches the chosen one. First entry is plus end, second is minus
        self.resonators[0,:] = [-a/2, 6*dy, -b-a/2,  4*dy]
        self.resonators[1,:] = [-a/2-b, 4*dy, -a/2,  2*dy]
        self.resonators[2,:] = [-a/2, 2*dy, -b-a/2,  0]
        self.resonators[3,:] = [-a/2-b, 0, -a/2,  -2*dy]
        self.resonators[4,:] = [-a/2, -2*dy, -b-a/2,  -4*dy]
        self.resonators[5,:] = [-a/2-b, -4*dy, -a/2,  -6*dy]
        self.resonators[6,:] = [a/2, -6*dy, -a/2,  -6*dy]
        
        self.resonators[7,:] = [-a/2, 2*dy, -a/4,  0]
        self.resonators[8,:] = [a/4,0, -a/4,  0]
        self.resonators[9,:] = [ a/2,  2*dy,a/4, 0]
        self.resonators[10,:] = [-a/4, 0, -a/2, -2*dy]
        self.resonators[11,:] = [a/4, 0, a/2, -2*dy]
#        self.resonators[7,:] = [-a/2, 2*dy, -a/2,  0]
#        self.resonators[8,:] = [a/2,0, -a/2,  0]
#        self.resonators[9,:] = [ a/2,  2*dy,a/2, 0]
#        self.resonators[10,:] = [-a/2, 0, -a/2, -2*dy]
#        self.resonators[11,:] = [a/2, 0, a/2, -2*dy]
        
        self.coords = self.get_coords(self.resonators)
        
        #set up the positions of the sites of the effective lattice
#        xs = numpy.asarray([-dx, -dx,-dx,-dx,-dx,-dx, 0, -xcorr, +xcorr, 0, -xcorr, +xcorr])
#        ys = numpy.asarray([5*dy, 3*dy, dy, -dy, -3*dy, -5*dy, -6*dy, 2*dy, 2*dy, 0, -2*dy, -2*dy])
        xs = numpy.zeros(self.numSites)
        ys = numpy.zeros(self.numSites)
        for rind in range(0, self.resonators.shape[0]):
            res = self.resonators[rind,:]
            xs[rind] = (res[0] + res[2])/2
            ys[rind] = (res[1] + res[3])/2
        self.SDx = xs
        self.SDy = ys
        
        
        #####auto population of the SD links
        self._auto_generate_SDlinks()
        
        #remove bad links from 4-way coupler
        badLinks= []
        for lind in range(self.SDHWlinks.shape[0]):
            link = self.SDHWlinks[lind,:]
            site1 = link[0]
            site2 = link[1]
            if site1 ==7 and site2 ==10:
                badLinks.append(lind)
            if site1 ==10 and site2 ==7:
                badLinks.append(lind)
            if site1 ==9 and site2 ==8:
                badLinks.append(lind)
            if site1 ==8 and site2 ==9:
                badLinks.append(lind)
#        #mark the bad rows
#        self.SDHWlinks[badLinks,:] = -4
#        #excise the bad rows
#        self.SDHWlinks = self.SDHWlinks[~numpy.all(self.SDHWlinks == -4, axis=1)] 
#        #also store the old link format
#        oldlinks = self.SDHWlinks[:,0:4]
#        self.SDlinks = oldlinks
        
        
        #make note of which resonator you need in order to close the unit cell
        closure = {}
        
        #a1 direction (x)
        closure[(1,0)] =numpy.asarray([1,2,3,4,5])
        #-a1 direction (-x)
        closure[(-1,0)] =numpy.asarray([])
        
        #a2 direction (y)
        closure[(0,1)] =numpy.asarray([6])
        #-a2 direction (y)
        closure[(0,-1)] =numpy.asarray([])
        
         #a1,a2 direction (x,y)
        closure[(1,1)] =numpy.asarray([])
        #-a1,a2 direction (-x,y)
        closure[(-1,1)] =numpy.asarray([])
        #a1,-a2 direction (x,-y)
        closure[(1,-1)] =numpy.asarray([0])
        #-a1,-a2 direction (-x,-y)
        closure[(-1,-1)] =numpy.asarray([])
        self.closure = closure
        
        return
    
    def _generate_74Huse_cell(self, side = 1):        
        '''
        generate 7,4 Huse variant unit cell
        
        Note: this graph is not three regular, but it is triangle protected
        '''
        
        #set up the sites
        self.numSites = 11
        xs = numpy.zeros(self.numSites)
        ys = numpy.zeros(self.numSites)
        
        #set up the lattice vectors
        self.a1 = numpy.asarray([self.side*numpy.sqrt(3)/2, self.side/2])
        self.a2 = numpy.asarray([0, 3*self.side])
        dy = self.a1[1]/2
        dx = self.a1[0]/2
        xcorr = self.side/numpy.sqrt(3)/2/2
        
        #set up the poisitions of all the resonators  and their end points
        self.resonators = numpy.zeros((self.numSites,4)) #pairs of resonator end points for each resonator
        self.coords = numpy.zeros((self.numSites,2)) #set of all resoantor start points
        
        a = self.side/numpy.sqrt(3)
        b = self.a1[0]-a
        #xo,yo,x1,y1
        #define them so their orientation matches the chosen one. First entry is plus end, second is minus
        self.resonators[0,:] = [-a/2, 6*dy, -b-a/2,  4*dy]
        self.resonators[1,:] = [-a/2-b, 4*dy, -a/2,  2*dy]
        self.resonators[2,:] = [-a/2, 2*dy, -b-a/2,  0]
        self.resonators[3,:] = [-a/2-b, 0, -a/2,  -2*dy]
        self.resonators[4,:] = [-a/2, -2*dy, -b-a/2,  -4*dy]
        self.resonators[5,:] = [-a/2-b, -4*dy, -a/2,  -6*dy]
        self.resonators[6,:] = [a/2, -6*dy, -a/2,  -6*dy]
        
        self.resonators[7,:] = [-a/2, 2*dy, 0,  0]
        self.resonators[8,:] = [a/2, 2*dy, 0,  0]
        self.resonators[9,:] = [ 0,  0,-a/2, -2*dy]
        self.resonators[10,:] = [0, 0, a/2, -2*dy]
        
        self.coords = self.get_coords(self.resonators)
        
        #set up the positions of the sites of the effective lattice
#        xs = numpy.asarray([-dx, -dx,-dx,-dx,-dx,-dx, 0, -xcorr, +xcorr, 0, -xcorr, +xcorr])
#        ys = numpy.asarray([5*dy, 3*dy, dy, -dy, -3*dy, -5*dy, -6*dy, 2*dy, 2*dy, 0, -2*dy, -2*dy])
        xs = numpy.zeros(self.numSites)
        ys = numpy.zeros(self.numSites)
        for rind in range(0, self.resonators.shape[0]):
            res = self.resonators[rind,:]
            xs[rind] = (res[0] + res[2])/2
            ys[rind] = (res[1] + res[3])/2
        self.SDx = xs
        self.SDy = ys
        
        
        #####auto population of the SD links
        self._auto_generate_SDlinks()
        
        #remove bad links from 4-way coupler
        badLinks= []
        for lind in range(self.SDHWlinks.shape[0]):
            link = self.SDHWlinks[lind,:]
            site1 = link[0]
            site2 = link[1]
            if site1 ==7 and site2 ==10:
                badLinks.append(lind)
            if site1 ==10 and site2 ==7:
                badLinks.append(lind)
            if site1 ==9 and site2 ==8:
                badLinks.append(lind)
            if site1 ==8 and site2 ==9:
                badLinks.append(lind)
#        #mark the bad rows
#        self.SDHWlinks[badLinks,:] = -4
#        #excise the bad rows
#        self.SDHWlinks = self.SDHWlinks[~numpy.all(self.SDHWlinks == -4, axis=1)] 
#        #also store the old link format
#        oldlinks = self.SDHWlinks[:,0:4]
#        self.SDlinks = oldlinks
        
        
        #make note of which resonator you need in order to close the unit cell
        closure = {}
        
        #a1 direction (x)
        closure[(1,0)] =numpy.asarray([1,2,3,4,5])
        #-a1 direction (-x)
        closure[(-1,0)] =numpy.asarray([])
        
        #a2 direction (y)
        closure[(0,1)] =numpy.asarray([6])
        #-a2 direction (y)
        closure[(0,-1)] =numpy.asarray([])
        
         #a1,a2 direction (x,y)
        closure[(1,1)] =numpy.asarray([])
        #-a1,a2 direction (-x,y)
        closure[(-1,1)] =numpy.asarray([])
        #a1,-a2 direction (x,-y)
        closure[(1,-1)] =numpy.asarray([0])
        #-a1,-a2 direction (-x,-y)
        closure[(-1,-1)] =numpy.asarray([])
        self.closure = closure
        
        return
    
    def _generate_123Huse_cell(self, side = 1): 
        '''
        generate 12,3 Huse variant unit cell
        '''
        
        #set up the sites
        self.numSites = 9
        xs = numpy.zeros(self.numSites)
        ys = numpy.zeros(self.numSites)
        
        #set up the lattice vectors
        self.a1 = numpy.asarray([self.side*numpy.sqrt(3)/2, self.side/2])
        self.a2 = numpy.asarray([0, 3*self.side])
        dy = self.a1[1]/2
        
        #set up the poisitions of all the resonators  and their end points
        self.resonators = numpy.zeros((self.numSites,4)) #pairs of resonator end points for each resonator
        self.coords = numpy.zeros((self.numSites,2)) #set of all resoantor start points
        
        a = self.side/numpy.sqrt(3)
        b = self.a1[0]-a
        #xo,yo,x1,y1
        #define them so their orientation matches the chosen one. First entry is plus end, second is minus
        self.resonators[0,:] = [-a/2, 6*dy, -b-a/2,  4*dy]
        self.resonators[1,:] = [-a/2-b, 4*dy, -a/2,  2*dy]
        self.resonators[2,:] = [-a/2, 2*dy, -b-a/2,  0]
        self.resonators[3,:] = [-a/2-b, 0, -a/2,  -2*dy]
        self.resonators[4,:] = [-a/2, -2*dy, -b-a/2,  -4*dy]
        self.resonators[5,:] = [-a/2-b, -4*dy, -a/2,  -6*dy]
        self.resonators[6,:] = [a/2, -6*dy, -a/2,  -6*dy]
        
        self.resonators[7,:] = [-a/2, 2*dy, -a/2,  -2*dy]
        self.resonators[8,:] = [a/2, 2*dy, a/2,  -2*dy]
        
        self.coords = self.get_coords(self.resonators)
        
        #set up the positions of the sites of the effective lattice
#        xs = numpy.asarray([-dx, -dx,-dx,-dx,-dx,-dx, 0, -xcorr, +xcorr, 0, -xcorr, +xcorr])
#        ys = numpy.asarray([5*dy, 3*dy, dy, -dy, -3*dy, -5*dy, -6*dy, 2*dy, 2*dy, 0, -2*dy, -2*dy])
        xs = numpy.zeros(self.numSites)
        ys = numpy.zeros(self.numSites)
        for rind in range(0, self.resonators.shape[0]):
            res = self.resonators[rind,:]
            xs[rind] = (res[0] + res[2])/2
            ys[rind] = (res[1] + res[3])/2
        self.SDx = xs
        self.SDy = ys
        
#        #alternate drawing
#        self.SDx[7] = -a/15
#        self.SDx[8] = a/15
        
        
        #####auto population of the SD links
        self._auto_generate_SDlinks()
         
        #make note of which resonator you need in order to close the unit cell
        closure = {}
        
        #a1 direction (x)
        closure[(1,0)] =numpy.asarray([1,2,3,4,5])
        #-a1 direction (-x)
        closure[(-1,0)] =numpy.asarray([])
        
        #a2 direction (y)
        closure[(0,1)] =numpy.asarray([6])
        #-a2 direction (y)
        closure[(0,-1)] =numpy.asarray([])
        
         #a1,a2 direction (x,y)
        closure[(1,1)] =numpy.asarray([])
        #-a1,a2 direction (-x,y)
        closure[(-1,1)] =numpy.asarray([])
        #a1,-a2 direction (x,-y)
        closure[(1,-1)] =numpy.asarray([0])
        #-a1,-a2 direction (-x,-y)
        closure[(-1,-1)] =numpy.asarray([])
        self.closure = closure
        
        return
    
    def _generate_Hk_composite_cell(self, xtrans, ytrans, side = 1, defect_type = 'Huse'):
        '''
        make a composite unit cell.
        
        if xtrans =1, huse cells will touch in x direction
        if x trans = 2, Huse cells will have a collumn of hexagona inbetween them
        
        same for ytrans
        
        haven't yet taken care of the closure properly
        
        defect_type tells it what kind of variation to use. Currently accepts all Huse-type cells
        Hopefully will auto accept other tings as they are added, but some funny business may crop up
        if certtain special symetries or things aren't respected
        '''
        cellH = UnitCell(defect_type, side = side)
        cellk = UnitCell('kagome', side = side)
        self.cellH = cellH
        self.cellk = cellk
        
        #set up the sites
        self.numSites = cellH.numSites + (xtrans-1)*cellk.numSites*3 +  (xtrans)*(ytrans-1)*cellk.numSites*3
        
        #set up the lattice vectors
        self.a1 = cellH.a1*xtrans
        self.a2 = cellH.a2*ytrans
        
        #allocate for the resonators
        self.resonators = numpy.zeros((self.numSites,4)) #pairs of resonator end points for each resonator
        self.coords = numpy.zeros((self.numSites,2)) #set of all resoantor start points
        
        #maks for shifting resonators
        xmask = numpy.zeros((cellk.numSites,4))
        ymask = numpy.zeros((cellk.numSites,4))
        xmask[:,0] = 1
        xmask[:,2] = 1
        ymask[:,1] = 1
        ymask[:,3] = 1
        
        #compile all the resonators
        rind = 0
        #all the others
        for indx in range(0, xtrans):
            for indy in range(0, ytrans):
                if (indx==0) and (indy ==0):
                    #the Huse cell
                    self.resonators[rind:rind+cellH.resonators.shape[0],:] = cellH.resonators
                    rind =  rind+ cellH.resonators.shape[0]
                else:  
                    for subind in range(0,3):
                        xOffset = indx*cellH.a1[0] + indy*cellH.a2[0] + (subind-1)*cellk.a2[0]
                        yOffset = indx*cellH.a1[1] + indy*cellH.a2[1] + (subind-1)*cellk.a2[1]
                        
                        ress = cellk.resonators + xOffset*xmask + yOffset*ymask
                        self.resonators[rind:rind+ress.shape[0],:] = ress
                        rind =  rind+ ress.shape[0]
                        
        #####auto population of the SD links
        self._auto_generate_SDlinks()
        
        #set up the positions of the sites of the effective lattice
        x0 = self.resonators[:,0]
        y0 = self.resonators[:,1]
        x1 = self.resonators[:,2]
        y1 = self.resonators[:,3]
        self.SDx = (x0+x1)/2.
        self.SDy = (y0+y1)/2.
        
        self.coords = self.get_coords(self.resonators)
        
        #make note of which resonator you need in order to close the unit cell
        ####!!!!!!!!!! incomplete 
        closure = {}
        
        #a1 direction (x)
        closure[(1,0)] =numpy.asarray([])
        #-a1 direction (-x)
        closure[(-1,0)] =numpy.asarray([])
        
        #a2 direction (y)
        closure[(0,1)] =numpy.asarray([])
        #-a2 direction (y)
        closure[(0,-1)] =numpy.asarray([])
        
         #a1,a2 direction (x,y)
        closure[(1,1)] =numpy.asarray([])
        #-a1,a2 direction (-x,y)
        closure[(-1,1)] =numpy.asarray([])
        #a1,-a2 direction (x,-y)
        closure[(1,-1)] =numpy.asarray([])
        #-a1,-a2 direction (-x,-y)
        closure[(-1,-1)] =numpy.asarray([])
        self.closure = closure
            
        
        return
    
    def _generate_arbitrary_cell(self, resonators, a1 = [1,0], a2 = [0,1]):
        '''
        generate arbitrary unit cell
        
        it needs to take in a set of resonators
        and possibly reciprocal lattice vectors
        
        it will multiply everything by self.side, so make sure resonators agrees with a1, and a2
        '''
        if resonators == '':
            raise ValueError, 'not a built-in unit cell type and no resonators given'
        else:
#            print resonators.shape
            if resonators.shape[1] != 4:
                raise ValueError, 'provided resonators are not the right shape'
            
        if a1.shape != (2,):
            raise ValueError, 'first reciprocal lattice vector has invalid shape'
            
        if a2.shape != (2,):
            raise ValueError, 'first reciprocal lattice vector has invalid shape'
        
        #set up the sites
        self.numSites = resonators.shape[0]
        
        #set up the lattice vectors
        self.a1 = numpy.asarray(a1)
        self.a2 = numpy.asarray(a2)
        
        
        #set up the poisitions of all the resonators  and their end points
        self.resonators = numpy.zeros((self.numSites,4)) #pairs of resonator end points for each resonator
        self.resonators=resonators*self.side
        
        self.coords = self.get_coords(self.resonators)
        
        #set up the positions of the sites of the effective lattice
        x0 = self.resonators[:,0]
        y0 = self.resonators[:,1]
        x1 = self.resonators[:,2]
        y1 = self.resonators[:,3]
        self.SDx = (x0+x1)/2.
        self.SDy = (y0+y1)/2.
        
        #####auto population of the SD links
        self._auto_generate_SDlinks()
        
        
        
        #make note of which resonator you need in order to close the unit cell
        #this is not handled well with an arbitrary cell
        closure = {}
        
        #a1 direction (x)
        closure[(1,0)] =numpy.asarray([])
        #-a1 direction (-x)
        closure[(-1,0)] =numpy.asarray([])
        
        #a2 direction (y)
        closure[(0,1)] =numpy.asarray([])
        #-a2 direction (y)
        closure[(0,-1)] =numpy.asarray([])
        
         #a1,a2 direction (x,y)
        closure[(1,1)] =numpy.asarray([])
        #-a1,a2 direction (-x,y)
        closure[(-1,1)] =numpy.asarray([])
        #a1,-a2 direction (x,-y)
        closure[(1,-1)] =numpy.asarray([])
        #-a1,-a2 direction (-x,-y)
        closure[(-1,-1)] =numpy.asarray([])
        self.closure = closure
        
        
        return
    
    def get_coords(self, resonators, roundDepth = 3):
        '''
        take in a set of resonators and calculate the set of end points.
        
        Will round all coordinates the the specified number of decimals.
        
        Should remove all redundancies.
        '''
        
        coords_overcomplete = numpy.zeros((resonators.shape[0]*2, 1)).astype('complex')
        coords_overcomplete =  numpy.concatenate((resonators[:,0], resonators[:,2])) + 1j * numpy.concatenate((resonators[:,1], resonators[:,3]))
        
        coords_complex = numpy.unique(numpy.round(coords_overcomplete, roundDepth))
    
        coords = numpy.zeros((coords_complex.shape[0],2))
        coords[:,0] = numpy.real(coords_complex)
        coords[:,1] = numpy.imag(coords_complex)
        
        return coords
    

    #######
    #draw functions
    #######  
    def draw_resonators(self, ax, color = 'g', alpha = 1 , linewidth = 0.5, zorder = 1):
        '''
        draw each resonator as a line
        '''
        for res in range(0,self.resonators.shape[0] ):
            [x0, y0, x1, y1]  = self.resonators[res,:]
            ax.plot([x0, x1],[y0, y1] , color = color, alpha = alpha, linewidth = linewidth, zorder = zorder)
        return
    
    def draw_resonator_end_points(self, ax, color = 'g', edgecolor = 'k',  marker = 'o' , size = 10, zorder = 1):
        '''will double draw some points'''
        x0s = self.resonators[:,0]
        y0s = self.resonators[:,1]
        
        x1s = self.resonators[:,2]
        y1s = self.resonators[:,3]
        
        pylab.sca(ax)
        pylab.scatter(x0s, y0s ,c =  color, s = size, marker = marker, edgecolors = edgecolor, zorder = zorder)
        pylab.scatter(x1s, y1s ,c =  color, s = size, marker = marker, edgecolors = edgecolor, zorder = zorder)
        return
      
    def draw_sites(self, ax, color = 'g', edgecolor = 'k',  marker = 'o' , size = 10, zorder=1):
        '''
        draw sites of the semidual (effective lattice)
        '''
        xs = self.SDx
        ys = self.SDy
        pylab.sca(ax)
        pylab.scatter(xs, ys ,c =  color, s = size, marker = marker, edgecolors = edgecolor, zorder = zorder)
        ax.set_aspect('equal')
        return
    
    def draw_SDlinks(self, ax, color = 'firebrick', linewidth = 0.5, HW = False, minus_color = 'goldenrod', zorder = 1, alpha = 1):
        '''
        draw all the links of the semidual lattice
        
        if extra is True it will draw only the edge sites required to fix the edge of the tiling
        
        set HW to true if you want the links color coded by sign
        minus_color sets the sign of the negative links
        '''

        links = self.SDHWlinks[:]
        
        for link in range(0, links.shape[0]):
            [startSite, endSite, deltaA1, deltaA2]  = links[link,0:4]
            startSite = int(startSite)
            endSite = int(endSite)
            
            [x0,y0] = [self.SDx[startSite], self.SDy[startSite]]
            [x1,y1] = numpy.asarray([self.SDx[endSite], self.SDy[endSite]]) + deltaA1*self.a1 + deltaA2*self.a2
            
            if HW:
                ends = links[link,4:6]
                if ends[0]==ends[1]:
                    #++ or --, use normal t
                    ax.plot([x0, x1],[y0, y1] , color = color, linewidth = linewidth, zorder = zorder, alpha = alpha)
                else:
                    #+- or -+, use inverted t
                    ax.plot([x0, x1],[y0, y1] , color = minus_color, linewidth = linewidth, zorder = zorder, alpha = alpha)
            else :
                ax.plot([x0, x1],[y0, y1] , color = color, linewidth = linewidth, zorder = zorder, alpha = alpha)
                
        return
    
    def _get_orientation_plot_points(self,scaleFactor = 0.5):
        '''
        find end coordinate locations part way along each resonator so that
        they can be used to plot the field at both ends of the resonator.
        
        Scale factor says how far appart the two points will be: +- sclaeFactor.2 of the total length
        
        returns the polt points as collumn matrix
        '''
        if scaleFactor> 1:
            raise ValueError, 'scale factor too big'
            
            
        size = len(self.SDx)
        plot_points = numpy.zeros((size*2, 2))
        
        resonators = self.resonators
        for ind in range(0, size):
            [x0, y0, x1, y1]  = resonators[ind, :]
            xmean = (x0+x1)/2
            ymean = (y0+y1)/2
            
            xdiff = x1-x0
            ydiff = y1-y0
            
            px0 = xmean - xdiff*scaleFactor/2
            py0 = ymean - ydiff*scaleFactor/2
            
            px1 = xmean + xdiff*scaleFactor/2
            py1 = ymean + ydiff*scaleFactor/2
            
            
            plot_points[2*ind,:] = [px0,py0]
            plot_points[2*ind+1,:] = [px1,py1]
            ind = ind+1
            
        return plot_points
    
    def draw_site_orientations(self,ax, title = 'state weight', colorbar = False, plot_links = False, cmap = 'jet_r', scaleFactor = 0.5, mSizes = 60, zorder = 1):
        Amps = numpy.ones(len(self.SDx))
        Probs = numpy.abs(Amps)**2
        mSizes = Probs * len(Probs)*30
        mColors = Amps
       
        mSizes = 60
        
        #build full state with value on both ends of the resonators 
        mColors_end = numpy.zeros(len(Amps)*2)
        mColors_end[0::2] = mColors

        #put opposite sign on other side
        mColors_end[1::2] = -mColors
#        mColors_end[1::2] = 5
        
        cm = pylab.cm.get_cmap(cmap)
        
        #get coordinates for the two ends of the resonator
        plotPoints = self._get_orientation_plot_points(scaleFactor = scaleFactor)
        xs = plotPoints[:,0]
        ys = plotPoints[:,1]
        
        pylab.sca(ax)
#        pylab.scatter(xs, ys,c =  mColors_end, s = mSizes, marker = 'o', edgecolors = 'k', cmap = cm, vmin = -1, vmax = 1, zorder = zorder)
        pylab.scatter(xs, ys,c =  mColors_end, s = mSizes, marker = 'o', edgecolors = 'k', cmap = cm, vmin = -1.5, vmax = 2.0, zorder = zorder)
        if colorbar:
            cbar = pylab.colorbar(fraction=0.046, pad=0.04)
            cbar.set_label('phase (pi radians)', rotation=270)
              
        if plot_links:
            self.draw_SDlinks(ax,linewidth = 0.5, color = 'firebrick', zorder = zorder)
        
        pylab.title(title, fontsize=8)
        ax.xaxis.set_visible(False)
        ax.yaxis.set_visible(False)
        ax.set_aspect('equal')
        
        return mColors_end
    
    #####
    #auto construction functions for SD links
    ######
    def _auto_generate_SDlinks(self):
        '''
        start from all the resonators of a unit cell auto generate the full link matrix,
        including neighboring cells
        '''
        xmask = numpy.zeros((self.numSites,4))
        ymask = numpy.zeros((self.numSites,4))
        
        xmask[:,0] = 1
        xmask[:,2] = 1
        
        ymask[:,1] = 1
        ymask[:,3] = 1
        
        if self.type[0:2] == '74':
            self.SDHWlinks = numpy.zeros((self.numSites*4+4,6))
        elif self.type == 'square':
            self.SDHWlinks = numpy.zeros((self.numSites*6,6))
        else:
#            self.SDHWlinks = numpy.zeros((self.numSites*4,6))
            self.SDHWlinks = numpy.zeros((self.numSites*8,6)) #temporary hack to allow some line graph games
        
        lind = 0
        for da1 in range(-1,2):
            for da2 in range(-1,2):
                links = self._auto_generate_cell_SDlinks(da1, da2)
                newLinks = links.shape[0]
                self.SDHWlinks[lind:lind+newLinks,:] = links
                lind = lind + newLinks
        
        #remove blank links (needed for some types of arbitrary cells)
        self.SDHWlinks = self.SDHWlinks[~numpy.all(self.SDHWlinks == 0, axis=1)] 
        
        #also store the old link format
        oldlinks = self.SDHWlinks[:,0:4]
        self.SDlinks = oldlinks 
        
        return
            
    def _auto_generate_cell_SDlinks(self, deltaA1, deltaA2):
        '''
        function to autogenerate the links between two sets of resonators
        deltaA1 and deltaA2 specify how many lattice vectors the two cells are seperated by
        in the first (~x) and second  (~y) lattice directions
        
        could be twice the same set, or it could be two different unit cells.
        
        will return a matrix of all the links [start, target, deltaA1, deltaA2, start_polarity, end_polarity]
        
        '''
        ress1 = self.resonators
        len1 = ress1.shape[0]
        
        #find the new unit cell
        xmask = numpy.zeros((self.numSites,4))
        ymask = numpy.zeros((self.numSites,4))
        xmask[:,0] = 1
        xmask[:,2] = 1
        ymask[:,1] = 1
        ymask[:,3] = 1
        xOffset = deltaA1*self.a1[0] + deltaA2*self.a2[0]
        yOffset = deltaA1*self.a1[1] + deltaA2*self.a2[1]
        ress2 = ress1 + xOffset*xmask + yOffset*ymask

        #place to store the links
        linkMat = numpy.zeros((len1*4+len1*4,6))
        
        #find the links
        
        #round the coordinates to prevent stupid mistakes in finding the connections
        plusEnds = numpy.round(ress2[:,0:2],3)
        minusEnds = numpy.round(ress2[:,2:4],3)
        
        extraLinkInd = 0
        for resInd in range(0,ress1.shape[0]):
            res = numpy.round(ress1[resInd,:],3)
            x1 = res[0]
            y1 = res[1]
            x0 = res[2]
            y0 = res[3]

            plusPlus = numpy.where((plusEnds == (x1, y1)).all(axis=1))[0]
            minusMinus = numpy.where((minusEnds == (x0, y0)).all(axis=1))[0]
            
            plusMinus = numpy.where((minusEnds == (x1, y1)).all(axis=1))[0] #plus end of new res, minus end of old
            minusPlus = numpy.where((plusEnds == (x0, y0)).all(axis=1))[0]
            
            for ind in plusPlus:
                if ind == resInd:
                    #self link
                    pass
                else:
                    linkMat[extraLinkInd,:] = [resInd, ind, deltaA1, deltaA2, 1,1]
                    extraLinkInd = extraLinkInd+1
                    
            for ind in minusMinus:
                if ind == resInd:
                    #self link
                    pass
                else:
                    linkMat[extraLinkInd,:] = [resInd, ind, deltaA1, deltaA2,  0,0]
                    extraLinkInd = extraLinkInd+1
                    
            for ind in plusMinus:
                linkMat[extraLinkInd,:] = [resInd, ind, deltaA1, deltaA2,  1,0]
                extraLinkInd = extraLinkInd+1
                
            for ind in minusPlus:
                linkMat[extraLinkInd,:] = [ resInd, ind, deltaA1, deltaA2,  0,1]
                extraLinkInd = extraLinkInd+1
        
        #clean the skipped links away 
        linkMat = linkMat[~numpy.all(linkMat == 0, axis=1)]  
        
        return linkMat
    
    ######
    #Bloch theory calculation functions
    ######
    def generate_Bloch_matrix(self, kx, ky, modeType = 'FW', t = 1, phase = 0):
        BlochMat = numpy.zeros((self.numSites, self.numSites))*(0 + 0j)
        
        for lind in range(0, self.SDHWlinks.shape[0]):
            link = self.SDHWlinks[lind,:]
            startInd = int(link[0]) #within the unit cell
            targetInd = int(link[1])
            deltaA1 = int(link[2])
            deltaA2   = int(link[3])
            startPol = int(link[4])
            targetPol = int(link[5])
            
            polarity = startPol^targetPol #xor of the two ends. Will be one when the two ends are different
            if phase == 0: #all the standard FW HW cases
                if modeType == 'HW':
                    signum =(-1.)**(polarity)
                elif modeType == 'FW':
                    signum = 1.
                else:
                    raise ValueError, 'Incorrect mode type. Must be FW or HW.'
            else: #artificially break TR symmetry
                if modeType == 'HW':
                    signum =(-1.)**(polarity)
                    if signum < 0:
                        if startInd > targetInd:
                            phaseFactor = numpy.exp(1j *phase) #e^i phi in one corner
                        elif startInd < targetInd:
                            phaseFactor = numpy.exp(-1j *phase) #e^-i phi in one corner, so it's Hermitian
                        else:
                            phaseFactor = 1
                            
                        signum = signum*phaseFactor
                        
                elif modeType == 'FW':
                    signum = 1.
                else:
                    raise ValueError, 'Incorrect mode type. Must be FW or HW.'
            
            #corrdiates of origin site
            x0 = self.SDx[startInd]
            y0 = self.SDy[startInd]
            
            #coordinates of target site
            x1 = self.SDx[targetInd] + deltaA1*self.a1[0] + deltaA2*self.a2[0]
            y1 = self.SDy[targetInd] + deltaA1*self.a1[1] + deltaA2*self.a2[1]
            
            deltaX = x1-x0
            deltaY = y1-y0
            
            phaseFactor = numpy.exp(1j*kx*deltaX)*numpy.exp(1j*ky*deltaY)
            BlochMat[startInd, targetInd] = BlochMat[startInd, targetInd]+ t*phaseFactor*signum
        return BlochMat
    
    def compute_band_structure(self, kx_0, ky_0, kx_1, ky_1, numsteps = 100, modeType = 'FW', returnStates = False, phase  = 0):
        '''
        from scipy.linalg.eigh:
        The normalized selected eigenvector corresponding to the eigenvalue w[i] is the column v[:,i].
        
        This returns same format with two additional kx, ky indices
        '''
        
        kxs = numpy.linspace(kx_0, kx_1,numsteps)
        kys = numpy.linspace(ky_0, ky_1,numsteps)
        
        bandCut = numpy.zeros((self.numSites, numsteps))
        
        stateCut = numpy.zeros((self.numSites, self.numSites, numsteps)).astype('complex')
        
        for ind in range(0, numsteps):
            kvec = [kxs[ind],kys[ind]]
            
            H = self.generate_Bloch_matrix(kvec[0], kvec[1], modeType = modeType, phase  = phase)
        
            #Psis = numpy.zeros((self.numSites, self.numSites)).astype('complex')
            Es, Psis = scipy.linalg.eigh(H)
            
            bandCut[:,ind] = Es
            stateCut[:,:,ind] = Psis
        if returnStates:
            return kxs, kys, bandCut, stateCut
        else:
            return kxs, kys, bandCut
    
    def plot_band_cut(self, ax, bandCut, colorlist = '', zorder = 1, dots = False, linewidth = 2.5):
        if colorlist == '':
            colorlist = ['firebrick', 'dodgerblue', 'blueviolet', 'mediumblue', 'goldenrod', 'cornflowerblue']
        
        pylab.sca(ax)
        
        for ind in range(0,self.numSites):
            colorInd = numpy.mod(ind, len(colorlist))
            if dots:
                pylab.plot(bandCut[ind,:], color = colorlist[colorInd] , marker = '.', markersize = '5', linestyle = '', zorder = zorder)
            else:
                pylab.plot(bandCut[ind,:], color = colorlist[colorInd] , linewidth = linewidth, linestyle = '-', zorder = zorder)
#            pylab.plot(bandCut[ind,:], '.')
        pylab.title('some momentum cut')
        pylab.ylabel('Energy')
        pylab.xlabel('k_something')
    
    def plot_bloch_wave(self, state_vect, ax, title = 'state weight', colorbar = False, plot_links = False, cmap = 'Wistia', zorder = 1):
        '''
        plot a state (wavefunction) on the graph of semidual points
        
        Only really works for full-wave solutions
        '''
        Amps = state_vect
        Probs = numpy.abs(Amps)**2
        mSizes = Probs * len(Probs)*30
        mColors = numpy.angle(Amps)/numpy.pi
        
        #move the branch cut to -0.5
        outOfRange = numpy.where(mColors< -0.5)[0]
        mColors[outOfRange] = mColors[outOfRange] + 2
        
        
        cm = pylab.cm.get_cmap(cmap)
        
        pylab.sca(ax)
        pylab.scatter(self.SDx, self.SDy,c =  mColors, s = mSizes, marker = 'o', edgecolors = 'k', cmap = cm, vmin = -0.5, vmax = 1.5, zorder = zorder)
        if colorbar:
            print 'making colorbar'
            cbar = pylab.colorbar(fraction=0.046, pad=0.04)
            cbar.set_label('phase (pi radians)', rotation=270)
              
        if plot_links:
            self.draw_SDlinks(ax, linewidth = 0.5, color = 'firebrick', zorder = zorder)
        
        pylab.title(title, fontsize=8)
        ax.xaxis.set_visible(False)
        ax.yaxis.set_visible(False)
        ax.set_aspect('equal')
        return
    
    def plot_bloch_wave_end_state(self, state_vect, ax, modeType, title = 'state weight', colorbar = False, plot_links = False, cmap = 'Wistia', scaleFactor = 0.5, zorder = 1):
        '''
        plot a state (wavefunction) on the graph of semidual points, but with a 
        value plotted for each end of the resonator
        
        If you just want a single value for the resonator use plot_layout_state
        
        Takes states defined on only one end of each resonator. Will autogenerate 
        the value on other end based on mode type.
        
        
        SOMETHING may be hinky with the range and flipping the sign
        
        '''
        Amps = state_vect
        Probs = numpy.abs(Amps)**2
        mSizes = Probs * len(Probs)*30
        mColors = numpy.angle(Amps)/numpy.pi
        
        #build full state with value on both ends of the resonators
        mSizes_end = numpy.zeros(len(Amps)*2)
        mSizes_end[0::2] = mSizes
        mSizes_end[1::2] = mSizes
        
        mColors_end = numpy.zeros(len(Amps)*2)
        mColors_end[0::2] = mColors
        if modeType == 'FW':
            mColors_end[1::2] = mColors
        elif modeType == 'HW':
            #put opposite phase on other side
            oppositeCols = mColors + 1
            #rectify the phases back to between -0.5 and 1.5 pi radians
            overflow = numpy.where(oppositeCols > 1.5)[0]
            newCols = oppositeCols
            newCols[overflow] = oppositeCols[overflow] - 2
            
            mColors_end[1::2] = newCols
        else:
            raise ValueError, 'You screwed around with the mode type. It must be FW or HW.'
        
        cm = pylab.cm.get_cmap(cmap)
        
        #get coordinates for the two ends of the resonator
        plotPoints = self._get_orientation_plot_points(scaleFactor = scaleFactor)
        xs = plotPoints[:,0]
        ys = plotPoints[:,1]
        
        pylab.sca(ax)
        pylab.scatter(xs, ys,c =  mColors_end, s = mSizes_end, marker = 'o', edgecolors = 'k', cmap = cm, vmin = -0.5, vmax = 1.5, zorder = zorder)
        if colorbar:
            print 'making colorbar'
            cbar = pylab.colorbar(fraction=0.046, pad=0.04)
            cbar.set_label('phase (pi radians)', rotation=270)
              
        if plot_links:
            self.draw_SDlinks(ax,linewidth = 0.5, color = 'firebrick', zorder = zorder)
        
        pylab.title(title, fontsize=8)
        ax.xaxis.set_visible(False)
        ax.yaxis.set_visible(False)
        ax.set_aspect('equal')
        
#        return plotPoints
        return mColors
    
    def split_cell(self, splitIn = 2, name = 'TBD'):
        resMat = self.resonators
        
        oldNum = resMat.shape[0]
    
        if type(splitIn) != int:
            raise ValueError, 'need an integer split'
        newNum = oldNum*splitIn
        
        newResonators = numpy.zeros((newNum,4))
        
        for rind in range(0, oldNum):
            oldRes = resMat[rind,:]
            xstart = oldRes[0]
            ystart = oldRes[1]
            xend = oldRes[2]
            yend = oldRes[3]
            
            xs = numpy.linspace(xstart, xend, splitIn+1)
            ys = numpy.linspace(ystart, yend, splitIn+1)
            for sind in range(0, splitIn):
                newResonators[splitIn*rind + sind,:] = [xs[sind], ys[sind], xs[sind+1], ys[sind+1]]
                
        newCell = UnitCell(name, resonators = newResonators, a1 = self.a1, a2 = self.a2)
        return newCell
    
    def line_graph_cell(self, name = 'TBD', resonatorsOnly = False):
        newResonators = numpy.zeros((self.SDHWlinks.shape[0], 4))
        
        for lind in range(0, self.SDHWlinks.shape[0]):
            link = self.SDHWlinks[lind,:]
            startInd = int(link[0]) #within the unit cell
            targetInd = int(link[1])
            
            deltaA1 = int(link[2])
            deltaA2   = int(link[3])
            
            startPol = int(link[4])
            targetPol = int(link[5])
            
            if (deltaA1,deltaA2) == (-1,1):
#                print 'skipping -1,1'
                pass
            elif (deltaA1,deltaA2) == (-1,0):
#                print 'skipping -1,0'
                pass
            elif (deltaA1,deltaA2) == (-1,-1):
#                print 'skipping -1,-1'
                pass
            elif (deltaA1,deltaA2) == (0,-1):
#                print 'skipping 0,-1'
                pass
            else:
                if (deltaA1,deltaA2) == (0,0) and  startInd > targetInd:
                    pass
                    #don't want to double count going the other way within the cell
                    #links to neighboring cells won't get double counted in this same way
                else:
                    #corrdiates of origin site
                    x0 = self.SDx[startInd]
                    y0 = self.SDy[startInd]
                    
                    #coordinates of target site
                    x1 = self.SDx[targetInd] + deltaA1*self.a1[0] + deltaA2*self.a2[0]
                    y1 = self.SDy[targetInd] + deltaA1*self.a1[1] + deltaA2*self.a2[1]
                    
                    res = numpy.asarray([x0, y0, x1, y1])
                    newResonators[lind, :] = res
                    
        
        #clean out balnk rows that were for redundant resonators
        newResonators = newResonators[~numpy.all(newResonators == 0, axis=1)]  

        if resonatorsOnly:
            return newResonators
        else:
            newCell = UnitCell(name, resonators = newResonators, a1 = self.a1, a2 = self.a2)
            return newCell
        
        
        

class EuclideanLayout(object):
    def __init__(self, xcells = 4, ycells = 4, lattice_type = 'Huse', side = 1, file_path = '', modeType = 'FW', resonatorsOnly=False, initialCell = ''):
        '''
        
        '''
        
        if file_path != '':
            self.load(file_path)
        else:
            #create plank planar layout object with the bare bones that you can build on
            self.xcells = xcells
            self.ycells = ycells
            self.side = side*1.0

            self.lattice_type = lattice_type
            
            if type(initialCell) == UnitCell:
                #use the unit cell object provided
                self.unitcell = initialCell
            else:
                #use a built in unit cell specified by keyword
                #starting unit cell
                self.unitcell = UnitCell(self.lattice_type, self.side)

            
            if not ((modeType == 'FW') or (modeType  == 'HW')):
                raise ValueError, 'Invalid mode type. Must be FW or HW'
            self.modeType = modeType
            
            self.populate(resonatorsOnly)
            
            
    ###########
    #automated construction, saving, loading
    ##########
    def populate(self, resonatorsOnly=False, Hamiltonian = True, save = False, save_name = ''):
        '''
        fully populate the structure up to itteration = MaxItter
        
        if Hamiltonian = False will not generate H
        save is obvious
        '''
         
        #make the resonator lattice
        self.generate_lattice(self.xcells, self.ycells)
        
        if not resonatorsOnly:
            #make the JC-Hubbard lattice
            self.generate_semiduals()
            
            if Hamiltonian:
                self.generate_Hamiltonian()
            
            if save:
                self.save(save_name)
            
        return
    
    def save(self, name = ''):
        '''
        save structure to a pickle file
        
        if name is blank, will use dafualt name
        '''
        if self.modeType == 'HW':
            waveStr = 'HW'
        else:
            waveStr = ''
            
        if name == '':
            name = str(self.lattice_type) + '_' + str(self.xcells) + 'x ' + str(self.ycells) + '_' + waveStr + '.pkl'
        
        savedict = self.__dict__
        pickle.dump(savedict, open(name, 'wb'))
        return
    
    def load(self, file_path):
        '''
        laod structure from pickle file
        '''
        pickledict = pickle.load(open(file_path, "rb" ) )
        
        for key in pickledict.keys():
            setattr(self, key, pickledict[key])
           
        #handle the case of old picle files that do not have a mode type property  
        #they are all calculated for the full wave
        if not 'modeType' in self.__dict__.keys():
            print 'Old pickle file. Pre FW-HW.'
            self.modeType = 'FW'
        return    
            
    #######
    #resonator lattice get /view functions
    #######       
    def get_xs(self):
        '''
        return x coordinates of all the resonator end points
        '''
        return self.coords[:,0]
    
    def get_ys(self):
        '''
        return y coordinates of all the resonator end points
        '''
        return self.coords[:,1]    
    
    def draw_resonator_lattice(self, ax, color = 'g', alpha = 1 , linewidth = 0.5, extras = False, zorder = 1):
        if extras == True:
            resonators = self.extraResonators
        else:
            resonators = self.resonators
            
        for res in range(0,resonators.shape[0] ):
            [x0, y0, x1, y1]  = resonators[res,:]
            ax.plot([x0, x1],[y0, y1] , color = color, alpha = alpha, linewidth = linewidth, zorder = zorder)
        return
            
    def draw_resonator_end_points(self, ax, color = 'g', edgecolor = 'k',  marker = 'o' , size = 10, zorder = 1):
        '''will double draw some points'''
        x0s = self.resonators[:,0]
        y0s = self.resonators[:,1]
        
        x1s = self.resonators[:,2]
        y1s = self.resonators[:,3]
        
        pylab.sca(ax)
        pylab.scatter(x0s, y0s ,c =  color, s = size, marker = marker, edgecolors = edgecolor, zorder = zorder)
        pylab.scatter(x1s, y1s ,c =  color, s = size, marker = marker, edgecolors = edgecolor, zorder = zorder)
        return   

    def get_all_resonators(self, maxItter = -1):
        '''
        function to get all resonators as a pair of end points
        
        each resontator returned as a row with four entries.
        (orientation is important to TB calculations)
        x0,y0,x1,y1
        
        '''
        return self.resonators
    
    def get_coords(self, resonators, roundDepth = 3):
        '''
        take in a set of resonators and calculate the set of end points.
        
        Will round all coordinates the the specified number of decimals.
        
        Should remove all redundancies.
        '''
        
        coords_overcomplete = numpy.zeros((resonators.shape[0]*2, 1)).astype('complex')
        coords_overcomplete =  numpy.concatenate((resonators[:,0], resonators[:,2])) + 1j * numpy.concatenate((resonators[:,1], resonators[:,3]))
        
        coords_complex = numpy.unique(numpy.round(coords_overcomplete, roundDepth))
    
        coords = numpy.zeros((coords_complex.shape[0],2))
        coords[:,0] = numpy.real(coords_complex)
        coords[:,1] = numpy.imag(coords_complex)
        
        return coords
    
    
    ########
    #functions to generate the resonator lattice
    #######
    def generate_lattice(self, xsize = -1, ysize = -1):
        '''
        Hopefully will become a general function to fill out the lattice. Has some issues
        with the edges right now
        
        it is important not to reverse the order of the endpoints of the lattice. These indicate
        the orinetation of the site, and will be needed to fill in the extra links to fix the edge of the lattice
        in HW mode
        '''
        if xsize <0:
            xsize = self.xcells
        if ysize <0:
            ysize = self.ycells
            
        #make sure that the object has the right size recorded
        self.xcells = xsize
        self.ycells = ysize
        
        self.resonators = numpy.zeros((xsize*ysize*self.unitcell.numSites, 4))
        
        #need a place to store extra resonators that live on the edge of the lattice
        self.extraResonators = numpy.zeros((xsize*ysize*self.unitcell.numSites, 4)) 
        self.extraSDx = numpy.zeros(xsize*ysize*self.unitcell.numSites) #these are easier to calculate here rather than later even though they don't fit thematically
        self.extraSDy = numpy.zeros(xsize*ysize*self.unitcell.numSites)
        
        xmask = numpy.zeros((self.unitcell.numSites,4))
        ymask = numpy.zeros((self.unitcell.numSites,4))
        
        xmask[:,0] = 1
        xmask[:,2] = 1
        
        ymask[:,1] = 1
        ymask[:,3] = 1
        
        ind = 0
        extraInd = 0
        for indx in range(0,xsize):
            for indy in range(0,ysize):
                xOffset = indx*self.unitcell.a1[0] + indy*self.unitcell.a2[0]
                yOffset = indx*self.unitcell.a1[1] + indy*self.unitcell.a2[1]
                self.resonators[ind:ind+self.unitcell.numSites, :] = self.unitcell.resonators + xOffset*xmask + yOffset*ymask
                
                if indy ==0:
                    #bottom row of lattice sites
                    xx = 0
                    yy = -1
                    indOut = self._fix_edge_resonators(extraInd,indx, indy, xx, yy)
                    extraInd = indOut
                    
                    xx = 1
                    yy = -1
                    indOut = self._fix_edge_resonators(extraInd,indx, indy, xx, yy)
                    extraInd = indOut
                
                if indx ==(xsize-1):
                    #right hand edge
                    xx = 1
                    yy = 0
                    indOut = self._fix_edge_resonators(extraInd,indx, indy, xx, yy)
                    extraInd = indOut
                    
                    xx = 1
                    yy = 1
                    indOut = self._fix_edge_resonators(extraInd,indx, indy, xx, yy)
                    extraInd = indOut
                    
                    if indy != 0:
                        xx = 1
                        yy = -1
                        indOut = self._fix_edge_resonators(extraInd,indx, indy, xx, yy)
                        extraInd = indOut
                
                if indy ==(ysize-1):
                    #top row of lattice sites
                    xx = 0
                    yy = 1
                    indOut = self._fix_edge_resonators(extraInd,indx, indy, xx, yy)
                    extraInd = indOut
                    
                    xx = -1
                    yy = 1
                    indOut = self._fix_edge_resonators(extraInd,indx, indy, xx, yy)
                    extraInd = indOut
                    
                if indy ==0:
                    #left0-hand edge
                    xx = -1
                    yy = 0
                    indOut = self._fix_edge_resonators(extraInd,indx, indy, xx, yy)
                    extraInd = indOut
                    
                    xx = -1
                    yy = -1
                    indOut = self._fix_edge_resonators(extraInd,indx, indy, xx, yy)
                    extraInd = indOut
                
                ind = ind + self.unitcell.numSites
        
        #clean the blank resonators away
        self.extraResonators = self.extraResonators[~numpy.all(self.extraResonators == 0, axis=1)]
        self.extraSDx = self.extraSDx[0:extraInd] 
        self.extraSDy = self.extraSDy[0:extraInd]
        
        #combine regular unit cell resonators and the extra edge ones
        self.resonators = numpy.concatenate((self.resonators, self.extraResonators))
        
        self.coords = self.get_coords(self.resonators)
        
        return
    
    def _fix_edge_resonators(self, extraInd,indx, indy, xx, yy):
        ''' very dirty function for fixing and adding missing resonators at the edge of the lattice
        mess with it at your own peril.
        
        Should only be called internally from generat_lattice
        
        also generates the extra center points of the resonators, because it's easier to do it all together.
        '''
        xmask = numpy.zeros((self.unitcell.numSites,4))
        ymask = numpy.zeros((self.unitcell.numSites,4))
        
        xmask[:,0] = 1
        xmask[:,2] = 1
        
        ymask[:,1] = 1
        ymask[:,3] = 1
        
        tempSites = self.unitcell.closure[(xx,yy)]
        xOffset = (indx+xx)*self.unitcell.a1[0] + (indy+yy)*self.unitcell.a2[0]
        yOffset = (indx+xx)*self.unitcell.a1[1] + (indy+yy)*self.unitcell.a2[1]
        for site in tempSites:
            self.extraResonators[extraInd, :] = self.unitcell.resonators[site,:] + xOffset*xmask[0,:] + yOffset*ymask[0,:]
            
            self.extraSDx[extraInd] = self.unitcell.SDx[site] + xOffset
            self.extraSDy[extraInd] = self.unitcell.SDy[site] + yOffset
            extraInd = extraInd + 1
            
        return extraInd
        
        

    ########
    #functions to generate effective JC-Hubbard lattice
    ########
    def generate_semiduals(self):
        '''
        main workhorse function to generate the JC-Hubbard lattice.
        This is the one you shold call. All the others are workhorses that it uses.
        
        Will loop through the existing and create attributes for the 
        JC-Hubbard lattice (here jokingly called semi-dual) and fill them
        '''
        xsize = self.xcells
        ysize = self.ycells
        
        self.SDx = numpy.zeros(xsize*ysize*self.unitcell.numSites)
        self.SDy = numpy.zeros(xsize*ysize*self.unitcell.numSites)
        
        #self.SDlinks = numpy.zeros((xsize*ysize*self.unitcell.numSites*4, 2))
        
        if self.lattice_type == 'square':
            self.SDHWlinks = numpy.zeros((xsize*ysize*self.unitcell.numSites*6, 4))
        
            self.extraSDHWlinks = numpy.zeros((xsize*ysize*self.unitcell.numSites*6, 4))
        else:
#            self.SDHWlinks = numpy.zeros((xsize*ysize*self.unitcell.numSites*4, 4))
#            
#            self.extraSDHWlinks = numpy.zeros((xsize*ysize*self.unitcell.numSites*4, 4))
            
            #temporary hack to allow playing with larger coordination numbers. Otherwise
            #there was not enough space allocated
            #what is really needed is something where the max coordination number of the unit
            #cell is used to do this properly.
            self.SDHWlinks = numpy.zeros((xsize*ysize*self.unitcell.numSites*8, 4))
            
            self.extraSDHWlinks = numpy.zeros((xsize*ysize*self.unitcell.numSites*8, 4))
            
        
        
        #set up for getting the positions of the semidual points
        xmask = numpy.zeros((self.unitcell.numSites,4))
        ymask = numpy.zeros((self.unitcell.numSites,4))
        
        xmask[:,0] = 1
        xmask[:,2] = 1
        
        ymask[:,1] = 1
        ymask[:,3] = 1

        #links will be done by site index, which will include the unit cell number
        latticelinkInd = 0
        ind = 0
        for indx in range(0,xsize):
            for indy in range(0,ysize):
                currCell = [indx, indy]
                
                xOffset = indx*self.unitcell.a1[0] + indy*self.unitcell.a2[0]
                yOffset = indx*self.unitcell.a1[1] + indy*self.unitcell.a2[1]
                self.SDx[ind:ind+self.unitcell.numSites] = self.unitcell.SDx + xOffset
                self.SDy[ind:ind+self.unitcell.numSites] = self.unitcell.SDy + yOffset
                
                ind = ind + self.unitcell.numSites
                
                for link in range(0, self.unitcell.SDlinks.shape[0]):
                    [startSite, targetSite, deltaA1, deltaA2, startEnd, targetEnd]  = self.unitcell.SDHWlinks[link,:]
                    targetCell = [indx + deltaA1, indy + deltaA2]
#                    print [startSite, targetSite, deltaA1, deltaA2, startEnd, targetEnd]
#                    print currCell
#                    print targetCell
                    if (targetCell[0]<0) or (targetCell[1]<0) or (targetCell[0]>xsize-1) or (targetCell[1]>ysize-1):
                        #this cell is outside of the simulation. Leave it
#                        print 'passing by'
                        pass
                    else:
                        startInd = startSite + currCell[0]*self.unitcell.numSites*ysize + currCell[1]*self.unitcell.numSites
                        targetInd = targetSite + targetCell[0]*self.unitcell.numSites*ysize + targetCell[1]*self.unitcell.numSites
                        self.SDHWlinks[latticelinkInd,:] = [startInd, targetInd, startEnd, targetEnd]
#                        print [startInd, targetInd, startEnd, targetEnd]
                        latticelinkInd = latticelinkInd +1  
#                    print '   '
        
        #fix the edge
        self._fix_SDedge()
        
        #clean the skipped links away 
        self.SDHWlinks = self.SDHWlinks[~numpy.all(self.SDHWlinks == 0, axis=1)]  
        self.extraSDHWlinks = self.extraSDHWlinks[~numpy.all(self.extraSDHWlinks == 0, axis=1)] 
        
        #add the extra links to the lattice
        self.SDHWlinks = numpy.concatenate((self.SDHWlinks , self.extraSDHWlinks))
        
        #make the truncated SD links
        self.SDlinks = self.SDHWlinks[:,0:2]
        
        #add the extra sites to the lattice
        self.SDx = numpy.concatenate((self.SDx, self.extraSDx))
        self.SDy = numpy.concatenate((self.SDy, self.extraSDy))
        
        return
    
    def _fix_SDedge(self):
        '''function to loop over the extra edge resonators and add their links in '''
        originalLatticeSize = self.xcells*self.ycells*self.unitcell.numSites
        
        #round the coordinates to prevent stupid mistakes in finding the connections
        plusEnds = numpy.round(self.resonators[:,0:2],3)
        minusEnds = numpy.round(self.resonators[:,2:4],3)
        
        extraLinkInd = 0
        for resInd in range(0,self.extraResonators.shape[0]):
            res = numpy.round(self.extraResonators[resInd,:],3)
            x1 = res[0]
            y1 = res[1]
            x0 = res[2]
            y0 = res[3]

            plusPlus = numpy.where((plusEnds == (x1, y1)).all(axis=1))[0]
            minusMinus = numpy.where((minusEnds == (x0, y0)).all(axis=1))[0]
            
            plusMinus = numpy.where((minusEnds == (x1, y1)).all(axis=1))[0] #plus end of new res, minus end of old
            minusPlus = numpy.where((plusEnds == (x0, y0)).all(axis=1))[0]
            
            for ind in plusPlus:
                if ind == originalLatticeSize+ resInd:
                    #self link
                    pass
                else:
                    self.extraSDHWlinks[extraLinkInd,:] = [originalLatticeSize+ resInd, ind,  1,1]
                    extraLinkInd = extraLinkInd+1
                    
                    #reverse link
                    self.extraSDHWlinks[extraLinkInd,:] = [ind, originalLatticeSize+ resInd,  1,1]
                    extraLinkInd = extraLinkInd+1
            for ind in minusMinus:
                if ind == originalLatticeSize+ resInd:
                    #self link
                    pass
                else:
                    self.extraSDHWlinks[extraLinkInd,:] = [originalLatticeSize+ resInd, ind,  0,0]
                    extraLinkInd = extraLinkInd+1
                    
                    #reverse link
                    self.extraSDHWlinks[extraLinkInd,:] = [ind, originalLatticeSize+ resInd,  0,0]
                    extraLinkInd = extraLinkInd+1
                
            for ind in plusMinus:
                self.extraSDHWlinks[extraLinkInd,:] = [originalLatticeSize+ resInd, ind,  1,0]
                extraLinkInd = extraLinkInd+1
                
                #reverse link
                self.extraSDHWlinks[extraLinkInd,:] = [ind, originalLatticeSize+ resInd,  0,1]
                extraLinkInd = extraLinkInd+1
                
            for ind in minusPlus:
                self.extraSDHWlinks[extraLinkInd,:] = [originalLatticeSize+ resInd, ind,  0,1]
                extraLinkInd = extraLinkInd+1
                
                #reverse link
                self.extraSDHWlinks[extraLinkInd,:] = [ind, originalLatticeSize+ resInd,  1,0]
                extraLinkInd = extraLinkInd+1
            
        return
    
    def draw_SD_points(self, ax, color = 'g', edgecolor = 'k',  marker = 'o' , size = 10,  extra = False, zorder = 1):
        '''
        draw the locations of all the semidual sites
        
        if extra is True it will draw only the edge sites required to fix the edge of the tiling
        '''
        if extra == True:
            xs = self.extraSDx
            ys = self.extraSDy
        else:
            xs = self.SDx
            ys = self.SDy
        
        pylab.sca(ax)
        pylab.scatter(xs, ys ,c =  color, s = size, marker = marker, edgecolors = edgecolor, zorder = zorder)
        
        return

    def draw_SDlinks(self, ax, color = 'firebrick', linewidth = 0.5, extra = False, minus_links = False, minus_color = 'goldenrod', zorder = 1, alpha = 1):
        '''
        draw all the links of the semidual lattice
        
        if extra is True it will draw only the edge sites required to fix the edge of the tiling
        
        set minus_links to true if you want the links color coded by sign
        minus_color sets the sign of the negative links
        '''
        if extra == True:
            xs = self.SDx
            ys = self.SDy
            links = self.extraSDHWlinks[:]
        else:
            xs = self.SDx
            ys = self.SDy
            links = self.SDHWlinks[:]
        
        for link in range(0, links.shape[0]):
            [startInd, endInd]  = links[link,0:2]
            startInd = int(startInd)
            endInd = int(endInd)
            
            [x0,y0] = [xs[startInd], ys[startInd]]
            [x1,y1] = [xs[endInd], ys[endInd]]
            
            if  minus_links == True and self.modeType == 'HW':
                ends = links[link,2:4]
                if ends[0]==ends[1]:
                    #++ or --, use normal t
                    ax.plot([x0, x1],[y0, y1] , color = color, linewidth = linewidth, zorder = zorder, alpha = alpha)
                else:
                    #+- or -+, use inverted t
                    ax.plot([x0, x1],[y0, y1] , color = minus_color, linewidth = linewidth, zorder = zorder, alpha = alpha)
            else :
                ax.plot([x0, x1],[y0, y1] , color = color, linewidth = linewidth, zorder = zorder, alpha = alpha)
                
        return

    def get_semidual_points(self):
        '''
        get all the semidual points in a given itteration.
        
        Mostly vestigial for compatibility
        '''
        return[self.SDx, self.SDy]

    
    ######
    #Hamiltonian related methods
    ######
    def generate_Hamiltonian(self, t = 1, internalBond = 1000):
        '''
        create the effective tight-binding Hamiltonian
        
        Also calculated as stores eigenvectors and eigenvalues for that H
        
        
        Will use FW or HW TB coefficients depending on self.modeType
        '''
        

        self.t = t
        self.internalBond = 1000*self.t
        
        totalSize = len(self.SDx)
            
        self.H = numpy.zeros((totalSize, totalSize))
        self.H_HW = numpy.zeros((totalSize*2, totalSize*2)) #vestigial
        
        #loop over the links and fill the Hamiltonian
        for link in range(0, self.SDlinks.shape[0]):
            [sourceInd, targetInd] = self.SDlinks[link, :]
            [sourceEnd, targetEnd] = self.SDHWlinks[link, 2:]
            source = int(sourceInd)
            target = int(targetInd)
            sourceEnd = int(sourceEnd)
            targetEnd = int(targetEnd)
            
            
            if self.modeType == 'FW':
                self.H[source, target] = self.t
            elif self.modeType == 'HW':
                polarity = sourceEnd^targetEnd #xor of the two ends. Will be one when the two ends are different
                signum =(-1.)**(polarity) #will be zero when  two ends are same, and minus 1 otherwise
                self.H[source, target] = self.t * signum
            else:
                raise ValueError, 'You screwed around with the mode type. It must be FW or HW.'
            self.H_HW[2*source + sourceEnd, 2*target+targetEnd] = 2*self.t
                
        #fix the bonds between the two ends of the same site
        for site in range(0, totalSize):
            self.H_HW[2*site, 2*site+1] = self.internalBond
            self.H_HW[2*site+1, 2*site] = self.internalBond
                
        self.Es, self.Psis = scipy.linalg.eigh(self.H)
        self.Eorder = numpy.argsort(self.Es)
        
        return
    
    def get_eigs(self):
        '''
        returns eigenvectors and eigenvalues
        '''
        return [self.Es, self.Psis, self.Eorder]
    
#    def get_SDindex(self,num, itt, az = True):
#        '''
#        get the index location of a semidual point. 
#        
#        Point spcified by
#        something TBD
#        
#        (useful for making localized states at specific sites)
#        '''
#        
#        return currInd

    def build_local_state(self, site):
        '''
        build a single site state at any location on the lattice.
        
        site is the absolute index coordinate of the lattice site
        (use get_SDindex to obtain this in a halfway sensible fashion)
        '''
        if site >= len(self.SDx):
            raise ValueError, 'lattice doesnt have this many sites'
            
        state = numpy.zeros(len(self.SDx))*(0+0j)
        
        state[site] = 1.
        
        return state
    
    def V_int(self, ind1, ind2, states):
        '''
        Calculate total interaction enegery of two particles at lattice sites
        indexed by index 1 and index 2
        
        states is the set of eigenvectors that you want to include e.g. [0,1,2,3]
        '''
        psis_1 = self.Psis[ind1,states]
        psis_2 = self.Psis[ind2,states]
        
        return numpy.dot(numpy.conj(psis_2), psis_1)

    def V_int_map(self, source_ind, states = []):
        '''
        calculate a map of the interaction energy for a given location of the first
        qubit.
        Lattice sites specified by index in semidual points array
        
        must also specify which igenstates to include. default = all
        '''
        if states == []:
            states = scipy.arange(0, len(self.Es),1)
        
        int_vals = numpy.zeros(len(self.SDx))
        for ind2 in range(0, len(self.SDx)):
            int_vals[ind2] = self.V_int(source_ind, ind2,states)
        
        return int_vals
    
    def plot_layout_state(self, state_vect, ax, title = 'state weight', colorbar = False, plot_links = False, cmap = 'Wistia', zorder = 1):
        '''
        plot a state (wavefunction) on the graph of semidual points
        '''
        Amps = state_vect
        Probs = numpy.abs(Amps)**2
        mSizes = Probs * len(Probs)*30
        mColors = numpy.angle(Amps)/numpy.pi
        
        cm = pylab.cm.get_cmap(cmap)
        
        pylab.sca(ax)
        pylab.scatter(self.SDx, self.SDy,c =  mColors, s = mSizes, marker = 'o', edgecolors = 'k', cmap = cm, vmin = -0.5, vmax = 1.5, zorder = zorder)
        if colorbar:
            cbar = pylab.colorbar(fraction=0.046, pad=0.04)
            cbar.set_label('phase (pi radians)', rotation=270)
              
        if plot_links:
            self.draw_SDlinks(ax, linewidth = 0.5, color = 'firebrick', zorder = zorder)
        
        pylab.title(title, fontsize=8)
        ax.xaxis.set_visible(False)
        ax.yaxis.set_visible(False)
        ax.set_aspect('equal')
        return

    def plot_map_state(self, map_vect, ax, title = 'ineraction weight', colorbar = False, plot_links = False, cmap = 'winter', autoscale = False, scaleFactor = 0.5, zorder = 1):
        '''plot an interaction map on the graph
        '''
        Amps = map_vect
        
        mSizes = 100
        mColors = Amps
        
    #    cm = pylab.cm.get_cmap('seismic')
        cm = pylab.cm.get_cmap(cmap)
    #    cm = pylab.cm.get_cmap('RdBu')
        
        
        vals = numpy.sort(mColors)
        peak = vals[-1]
        second_biggest = vals[-2]
        
        if autoscale:
            vmax = peak
            if self.modeType == 'HW':
                vmin = -vmax
            else:
                vmin = vals[0]
        else:
            vmax = second_biggest
            vmin = -second_biggest
        
        if self.modeType == 'FW':
            pylab.sca(ax)
            pylab.scatter(self.SDx, self.SDy,c =  mColors, s = mSizes, marker = 'o', edgecolors = 'k', cmap = cm, vmax = vmax, vmin = vmin, zorder = zorder)
        elif self.modeType == 'HW':
            #build full state with value on both ends of the resonators
            mColors_end = numpy.zeros(len(Amps)*2)
            mColors_end[0::2] = mColors
            mColors_end[1::2] = -mColors
            
            #get coordinates for the two ends of the resonator
            plotPoints = self.get_end_state_plot_points(scaleFactor = scaleFactor)
            xs = plotPoints[:,0]
            ys = plotPoints[:,1]
            
#            mColors_end = scipy.arange(1.,len(Amps)*2+1,1)/300
#            print mColors_end.shape
#            print Amps.shape
            
            #plot
            pylab.sca(ax)
            pylab.scatter(xs, ys,c =  mColors_end, s = mSizes/1.4, marker = 'o', edgecolors = 'k', cmap = cm, vmax = vmax, vmin = vmin, zorder = zorder)
        else:
            raise ValueError, 'You screwed around with the mode type. It must be FW or HW.'
            
        
        if colorbar:
            cbar = pylab.colorbar(fraction=0.046, pad=0.04)
            cbar.set_label('interaction energy (AU)', rotation=270)
              
        if plot_links:
            self.draw_SDlinks(ax, linewidth = 0.5, color = 'firebrick', zorder = zorder)
        
        pylab.title(title, fontsize=8)
        ax.xaxis.set_visible(False)
        ax.yaxis.set_visible(False)
        ax.set_aspect('equal')
        return
    
    def get_end_state_plot_points(self,scaleFactor = 0.5):
        '''
        find end coordinate locations part way along each resonator so that
        they can be used to plot the field at both ends of the resonator.
        (Will retun all values up to specified itteration. Default is the whole thing)
        
        Scale factor says how far appart the two points will be: +- sclaeFactor.2 of the total length
        
        returns the polt points as collumn matrix
        '''
        if scaleFactor> 1:
            raise ValueError, 'scale factor too big'
            
            
        size = len(self.SDx)
        plot_points = numpy.zeros((size*2, 2))
        
        resonators = self.get_all_resonators()
        for ind in range(0, size):
            [x0, y0, x1, y1]  = resonators[ind, :]
            xmean = (x0+x1)/2
            ymean = (y0+y1)/2
            
            xdiff = x1-x0
            ydiff = y1-y0
            
            px0 = xmean - xdiff*scaleFactor/2
            py0 = ymean - ydiff*scaleFactor/2
            
            px1 = xmean + xdiff*scaleFactor/2
            py1 = ymean + ydiff*scaleFactor/2
            
            
            plot_points[2*ind,:] = [px0,py0]
            plot_points[2*ind+1,:] = [px1,py1]
            ind = ind+1
            
        return plot_points
    
    def plot_end_layout_state(self, state_vect, ax, title = 'state weight', colorbar = False, plot_links = False, cmap = 'Wistia', scaleFactor = 0.5, zorder = 1):
        '''
        plot a state (wavefunction) on the graph of semidual points, but with a 
        value plotted for each end of the resonator
        
        If you just want a single value for the resonator use plot_layout_state
        
        Takes states defined on only one end of each resonator. Will autogenerate 
        the value on other end based on mode type.
        
        '''
        Amps = state_vect
        Probs = numpy.abs(Amps)**2
        mSizes = Probs * len(Probs)*30
        mColors = numpy.angle(Amps)/numpy.pi
        
        #build full state with value on both ends of the resonators
        mSizes_end = numpy.zeros(len(Amps)*2)
        mSizes_end[0::2] = mSizes
        mSizes_end[1::2] = mSizes
        
        mColors_end = numpy.zeros(len(Amps)*2)
        mColors_end[0::2] = mColors
        if self.modeType == 'FW':
            mColors_end[1::2] = mColors
        elif self.modeType == 'HW':
            #put opposite phase on other side
            oppositeCols = mColors + 1
            #rectify the phases back to between -0.5 and 1.5 pi radians
            overflow = numpy.where(oppositeCols > 1.5)[0]
            newCols = oppositeCols
            newCols[overflow] = oppositeCols[overflow] - 2
            
            mColors_end[1::2] = newCols
        else:
            raise ValueError, 'You screwed around with the mode type. It must be FW or HW.'
        
        cm = pylab.cm.get_cmap(cmap)
        
        #get coordinates for the two ends of the resonator
        plotPoints = self.get_end_state_plot_points(scaleFactor = scaleFactor)
        xs = plotPoints[:,0]
        ys = plotPoints[:,1]
        
        pylab.sca(ax)
        pylab.scatter(xs, ys,c =  mColors_end, s = mSizes_end, marker = 'o', edgecolors = 'k', cmap = cm, vmin = -0.5, vmax = 1.5, zorder = zorder)
        if colorbar:
            cbar = pylab.colorbar(fraction=0.046, pad=0.04)
            cbar.set_label('phase (pi radians)', rotation=270)
              
        if plot_links:
            self.draw_SDlinks(ax,linewidth = 0.5, color = 'firebrick', zorder = zorder)
        
        pylab.title(title, fontsize=8)
        ax.xaxis.set_visible(False)
        ax.yaxis.set_visible(False)
        ax.set_aspect('equal')
        
#        return plotPoints
        return mColors
    
    
    
    
    
    
    
    
    
    
    

if __name__=="__main__":  
#    Cell = True
    Cell = False

    Lattice = True
#    Lattice = False    


    ####Cell mode sub options
    K0States = False  #display or not 


    ####Lattice mode sub options
    LatticeHamiltonian = False #display or not 
    LatticeInteractionMap = False #display or not 
    
    
    
    
    ##################################
    ##################################
    ##################################
    ##################################
    ###########
    #lattice testing  and examples
    ##################################
    ##################################
    ##################################
    ##################################
    if Cell:
        testCell = UnitCell('Huse')
        
        #pylab.rcParams.update({'font.size': 14})
        #pylab.rcParams.update({'font.size': 8})
        
        modeType = 'FW'
        #modeType = 'HW'
        
#        testCell = UnitCell('Huse')
#        testCell = UnitCell('74Huse')
#        testCell = UnitCell('84Huse')
#        testCell = UnitCell('123Huse')
#        testCell = UnitCell('kagome')
        
#        testCell = UnitCell('Huse2_1')
        #testCell = UnitCell('Huse2_2')
        #testCell = UnitCell('Huse3_1')
        #testCell = UnitCell('Huse3_3')
        
        #testCell = UnitCell('84Huse2_1')
        
        #testCell = UnitCell('PeterChain')
        #testCell = UnitCell('PaterChain_tail')
        
        ######
        #test the unit cell
        #######
        pylab.figure(1)
        pylab.clf()
        ax = pylab.subplot(1,2,1)
        testCell.draw_sites(ax)
        pylab.title('Sites of Huse Cell')
        
        ax = pylab.subplot(1,2,2)
        testCell.draw_sites(ax,color = 'goldenrod', edgecolor = 'k',  marker = 'o' , size = 20)
        testCell.draw_resonators(ax, color = 'cornflowerblue', linewidth = 1)
        testCell.draw_SDlinks(ax, color = 'firebrick', linewidth = 1)
        testCell.draw_resonator_end_points(ax, color = 'deepskyblue', edgecolor = 'k',  marker = 'o' , size = 20)
        pylab.title('Links of Unit Cell')
        pylab.show()
        
        
        ######
        #show the orientations
        ######
        #alternate version
        fig = pylab.figure(2)
        pylab.clf()
        ax = pylab.subplot(1,1,1)
        testCell.draw_resonators(ax, color = 'cornflowerblue', linewidth = 1)
        testCell.draw_resonator_end_points(ax, color = 'indigo', edgecolor = 'indigo',  marker = '+' , size = 20)
        testCell.draw_site_orientations(ax, title = 'unit cell convention', colorbar = False, plot_links = False, cmap = 'jet', scaleFactor = 0.5)
        testCell.draw_SDlinks(ax, linewidth = 1.5, HW = True , minus_color = 'goldenrod')
        pylab.title('site orientations : ' + testCell.type)
        #ax.set_aspect('auto')
        ax.set_aspect('equal')
        #    fig.savefig('HW.png', dpi = 200)
        
        pylab.show()
        

        
        #####
        #testing bloch theory
        ####
        
        Hmat = testCell.generate_Bloch_matrix(0,0,  modeType = modeType)
        pylab.figure(3)
        pylab.clf()
        ax = pylab.subplot(1,2,1)
        pylab.imshow(numpy.abs(Hmat))
        pylab.title('|H|')
        
        ax = pylab.subplot(1,2,2)
        pylab.imshow(numpy.real(Hmat - numpy.transpose(numpy.conj(Hmat))))
        pylab.title('H - Hdagger')
        
        pylab.show()
        
        
        
        #kx_x, ky_y, cutx = testCell.compute_band_structure(-2*numpy.pi, 0, 2*numpy.pi, 0, numsteps = 100, modeType = modeType)
        #kx_y, ky_y, cuty = testCell.compute_band_structure(0, -8./3*numpy.pi, 0, 8./3*numpy.pi, numsteps = 100, modeType = modeType)
        kx_x, ky_y, cutx = testCell.compute_band_structure(-2*numpy.pi, 0, 2*numpy.pi, 0, numsteps = 100, modeType = modeType)
        kx_y, ky_y, cuty = testCell.compute_band_structure(0, -2.5*numpy.pi, 0, 2.5*numpy.pi, numsteps = 100, modeType = modeType)
        
        fig2 = pylab.figure(4)
        pylab.clf()
        ax = pylab.subplot(1,2,1)
        testCell.plot_band_cut(ax, cutx)
        pylab.title('xcut')
        
        ax = pylab.subplot(1,2,2)
        testCell.plot_band_cut(ax, cuty)
        pylab.title('ycut')
        
        titleStr = testCell.type + ', modeType: ' + modeType + ' (Made with UnitCell class)' 
        pylab.suptitle(titleStr)
        
        pylab.show()
        

        
        #####
        #look at above gap state at k= 0
        #####
        if K0States:
            Es, Psis = scipy.linalg.eigh(Hmat)
            
            stateInd = 0
            aboveGap = Psis[:,stateInd]
            print Es[stateInd]
            print aboveGap
            
            pylab.figure(5)
            pylab.clf()
            
            ax = pylab.subplot(1,1,1)
            #testCell.draw_sites(ax,color = 'goldenrod', edgecolor = 'k',  marker = 'o' , size = 20)
            testCell.draw_SDlinks(ax, color = 'firebrick', linewidth = 1)
            testCell.draw_resonators(ax, color = 'cornflowerblue', linewidth = 1)
            testCell.draw_resonator_end_points(ax, color = 'deepskyblue', edgecolor = 'k',  marker = 'o' , size = 20)
            #testCell.plot_bloch_wave(aboveGap*2, ax, title = 'state weight', colorbar = False, plot_links = False, cmap = 'Wistia')
            temp = testCell.plot_bloch_wave_end_state(aboveGap*2, ax,modeType = modeType,  title = modeType + '_' + str(stateInd), colorbar = False, plot_links = False, cmap = 'Wistia')
            ax.set_aspect('equal')
            pylab.show()
            
            
            ####try to plot all the unit cell wave functions. Doesn't work very well. You can't see anything
            #pylab.figure(6)
            #pylab.clf()
            #for ind in range(0, testCell.numSites):
            #    ax = pylab.subplot(1,testCell.numSites,ind+1)
            #    testCell.draw_SDlinks(ax, color = 'firebrick', linewidth = 1)
            #    testCell.draw_resonators(ax, color = 'cornflowerblue', linewidth = 1)
            #    testCell.draw_resonator_end_points(ax, color = 'deepskyblue', edgecolor = 'k',  marker = 'o' , size = 20)
            ##    testCell.plot_bloch_wave(Psis[:,ind], ax, title = 'state weight', colorbar = False, plot_links = False, cmap = 'Wistia')
            #    testCell.plot_bloch_wave_end_state(Psis[:,ind], ax,modeType = modeType,  title = str(ind), colorbar = False, plot_links = False, cmap = 'Wistia')
            #    ax.set_aspect('equal')
            #pylab.show()
        else:
            pylab.figure(5)
            pylab.clf()
            
            pylab.figure(6)
            pylab.clf()
        
    
    
    
    
    
    ##################################
    ##################################
    ##################################
    ##################################
    ###########
    #lattice testing  and examples
    ##################################
    ##################################
    ##################################
    ##################################
    if Lattice:
    #    testLattice = EuclideanLayout(4,3,lattice_type = 'Huse', modeType = 'FW')
    #    testLattice = EuclideanLayout(2,1,lattice_type = 'Huse', modeType = 'FW')
        
    #    testLattice = EuclideanLayout(4,3,lattice_type = 'Huse', modeType = 'HW')
    #    testLattice = EuclideanLayout(4,2,lattice_type = 'Huse', modeType = 'HW')
    #    testLattice = EuclideanLayout(2,2,lattice_type = 'Huse', modeType = 'HW')
    #    testLattice = EuclideanLayout(1,1,lattice_type = 'Huse', modeType = 'HW')
    
    
    #    testLattice = EuclideanLayout(4,3,lattice_type = 'Huse', modeType = 'FW', side = 500)
        
        testLattice = EuclideanLayout(4,4,lattice_type = 'kagome', modeType = 'FW')
        
#        testLattice = EuclideanLayout(2,3,lattice_type = 'Huse2_1', modeType = 'FW')
    
#        testLattice = EuclideanLayout(1,1,lattice_type = '84Huse2_1', modeType = 'FW')
        
    #    testLattice = EuclideanLayout(2,1,lattice_type = '84Huse', modeType = 'FW')
#        testLattice = EuclideanLayout(4,3,lattice_type = '74Huse', modeType = 'FW')
#        testLattice = EuclideanLayout(4,3,lattice_type = '123Huse', modeType = 'FW')

#        testLattice = EuclideanLayout(3,3,lattice_type = 'square', modeType = 'FW')
    
        ######
        #test the unit cell
        #######
        pylab.figure(1)
        pylab.clf()
    
        ######
        #test the generate functions
        #######
    #    testLattice.generate_lattice()
    #    testLattice.generate_semiduals()
    #    testLattice.generate_Hamiltonian()
    
        debugMode = False
    #    debugMode = True
        
        ######
        #test the lattice and SD lattice constructions
        #######
        pylab.figure(2)
        pylab.clf()
        ax = pylab.subplot(1,2,1)
        testLattice.draw_resonator_lattice(ax, color = 'cornflowerblue', linewidth = 1)
        testLattice.draw_resonator_end_points(ax, color = 'deepskyblue', edgecolor = 'k',  marker = 'o' , size = 20)
        
        if debugMode:
            testLattice.draw_resonator_lattice(ax, color = 'indigo', linewidth = 1, extras = True)
            [x0, y0, x1, y1]  = testLattice.extraResonators[0,:]
    #        ax.plot([x0, x1],[y0, y1] , color = 'firebrick', alpha = 1, linewidth = 1)
            [x0, y0, x1, y1]  = testLattice.resonators[6,:]
    #        ax.plot([x0, x1],[y0, y1] , color = 'indigo', alpha = 1, linewidth = 1)
        
        pylab.title('Resonators of Huse Lattice')
        
        ax = pylab.subplot(1,2,2)
        testLattice.draw_SD_points(ax, color = 'deepskyblue', edgecolor = 'k',  marker = 'o' , size = 20)
        testLattice.draw_SDlinks(ax, color = 'firebrick', linewidth = 1)
        
        if debugMode:
            testLattice.draw_SD_points(ax, color = 'indigo', edgecolor = 'k',  marker = 'o' , size = 20, extra = True)
            testLattice.draw_SDlinks(ax, color = 'cornflowerblue', linewidth = 1, extra = True)
        #    pylab.scatter(testLattice.extraSDx,testLattice.extraSDy ,c =  'indigo', s = 25, marker ='o', edgecolors = 'k')
        pylab.title('Links of the Huse Lattice')
        pylab.show()
        
        
        ######
        #test the Hamiltonian
        #######
        eigNum = 168
        eigNum = 167
        eigNum = 0
        if LatticeHamiltonian:
            pylab.figure(3)
            pylab.clf()
            ax = pylab.subplot(1,2,1)
            pylab.imshow(testLattice.H,cmap = 'winter')
            pylab.title('Hamiltonian')
            
            ax = pylab.subplot(1,2,2)
            pylab.imshow(testLattice.H - numpy.transpose(testLattice.H),cmap = 'winter')
            pylab.title('H - Htranspose')
            pylab.show()
            
    
            
            xs = scipy.arange(0,len(testLattice.Es),1)
            eigAmps = testLattice.Psis[:,testLattice.Eorder[eigNum]]
            
            pylab.figure(4)
            pylab.clf()
            ax1 = pylab.subplot(1,2,1)
            pylab.plot(testLattice.Es, 'b.')
            pylab.plot(xs[eigNum],testLattice.Es[testLattice.Eorder[eigNum]], color = 'firebrick' , marker = '.', markersize = '10' )
            pylab.title('eigen spectrum')
            pylab.ylabel('Energy (t)')
            pylab.xlabel('eigenvalue number')
            
            ax2 = pylab.subplot(1,2,2)
            titleStr = 'eigenvector weight : ' + str(eigNum)
            testLattice.plot_layout_state(eigAmps, ax2, title = titleStr, colorbar = True, plot_links = True, cmap = 'Wistia')
            
            pylab.show()
        else:
            pylab.figure(3)
            pylab.clf()
            
            pylab.figure(4)
            pylab.clf()
        
        
        ######
        #test the layout plotters (center dot)
        #######
        
        pylab.figure(5)
        pylab.clf()
        stateInd = eigNum
        state1 = testLattice.Psis[:,stateInd]
        if testLattice.xcells < 4 and testLattice.ycells <3:
            state2 = testLattice.build_local_state(7)
        else:
#            state2 = testLattice.build_local_state(47)
            state2 = testLattice.build_local_state(4)
        
        
        ax = pylab.subplot(1,2,1)
        testLattice.plot_layout_state(state1, ax, title = 'eigenstate', colorbar = False, plot_links = True, cmap = 'Wistia')
        
        ax = pylab.subplot(1,2,2)
        testLattice.plot_layout_state(state2/10, ax, title = 'local state', colorbar = False, plot_links = True, cmap = 'Wistia')
        
        pylab.show()
        
        
        ######
        #test the interaction funtions
        #######
        if LatticeInteractionMap:
            #    interactionStates = scipy.arange(0,len(testLattice.Es),1)
            if testLattice.xcells < 4 and testLattice.ycells <3:
                interactionStates = scipy.arange(0,4,1)
                site1 = 1
                site2 = 5
            else:
                interactionStates = scipy.arange(0,47,1)
                site1 = 10
                site2 = 54
            
            
            
            V0 = testLattice.V_int(site1, site1, interactionStates)
            VV = testLattice.V_int(site1, site2, interactionStates)
            print V0
            print VV
            
            Vmap0 = testLattice.V_int_map(site2, interactionStates)
            Vmap1 = testLattice.V_int_map(site2, interactionStates[0:4])
            
            pylab.figure(6)
            pylab.clf()
            ax = pylab.subplot(1,2,1)
            testLattice.plot_map_state(Vmap0, ax, title = 'ineraction weight: all FB states, hopefully', colorbar = True, plot_links = True, cmap = 'winter', autoscale = False)
            pylab.scatter([testLattice.SDx[site2]], [testLattice.SDy[site2]], c =  'gold', s = 150, edgecolors = 'k')
            
            ax = pylab.subplot(1,2,2)
            testLattice.plot_map_state(Vmap1, ax, title = 'ineraction weight: first 4', colorbar = True, plot_links = True, cmap = 'winter', autoscale = False)
            pylab.scatter([testLattice.SDx[site2]], [testLattice.SDy[site2]], c =  'gold', s = 150, edgecolors = 'k')
            
            pylab.show()
        else:
            pylab.figure(6)
            pylab.clf()
        
        
        ######
        #test visualization functions for shwing both ends of the resonators
        #######
        state_uniform = numpy.ones(len(testLattice.SDx))/numpy.sqrt(len(testLattice.SDx))
        
        pylab.figure(7)
        pylab.clf()
        ax = pylab.subplot(1,2,1)
    #    testLattice.plot_layout_state(state1, ax, title = 'eigenstate', colorbar = False, plot_links = True, cmap = 'Wistia')
        testLattice.plot_layout_state(state_uniform, ax, title = 'eigenstate', colorbar = False, plot_links = True, cmap = 'Wistia')
        
        ax = pylab.subplot(1,2,2)
        endplot_points = testLattice.get_end_state_plot_points()
    #    testLattice.plot_end_layout_state(state1, ax, title = 'end weights', colorbar = False, plot_links = True, cmap = 'Wistia', scaleFactor = 0.5)
        testLattice.plot_end_layout_state(state_uniform, ax, title = 'end weights', colorbar = False, plot_links = True, cmap = 'Wistia', scaleFactor = 0.5)
        
        pylab.show()
        
        
        
    #    #####
    #    #checking conventions
    #    #####
    #    
    #    pylab.figure(17)
    #    pylab.clf()
    #    ax = pylab.subplot(1,2,1)
    #    testLattice.draw_resonator_lattice(ax, color = 'cornflowerblue', linewidth = 1)
    #    testLattice.draw_resonator_end_points(ax, color = 'indigo', edgecolor = 'indigo',  marker = '+' , size = 20)
    ##    testLattice.plot_end_layout_state(state_uniform, ax, title = 'end weights', colorbar = False, plot_links = False, cmap = 'Wistia', scaleFactor = 0.5)
    #    testLattice.plot_end_layout_state(state_uniform*1.4, ax, title = 'unit cell convention', colorbar = False, plot_links = False, cmap = 'jet', scaleFactor = 0.5)
    #    testLattice.draw_SDlinks(ax, linewidth = 1, extra = False, minus_links = True, minus_color = 'goldenrod')
    #    pylab.title('site orientations')
    #    
    #    ax = pylab.subplot(1,2,2)
    #    pylab.imshow(testLattice.H,cmap = 'winter')
    #    pylab.title('Hamiltonian')
    #    pylab.show()
    #    
    #    pylab.figure(19)
    #    pylab.clf()
    #    ax = pylab.subplot(1,1,1)
    #    testLattice.draw_resonator_lattice(ax, color = 'cornflowerblue', linewidth = 1)
    #    testLattice.draw_resonator_end_points(ax, color = 'indigo', edgecolor = 'indigo',  marker = '+' , size = 20)
    ##    testLattice.plot_end_layout_state(state_uniform, ax, title = 'end weights', colorbar = False, plot_links = False, cmap = 'Wistia', scaleFactor = 0.5)
    #    testLattice.plot_end_layout_state(state_uniform*1.4, ax, title = 'unit cell convention', colorbar = False, plot_links = False, cmap = 'jet', scaleFactor = 0.5)
    #    testLattice.draw_SDlinks(ax, linewidth = 1, extra = False, minus_links = True, minus_color = 'goldenrod')
    #    pylab.title('site orientations')
    #    ax.set_aspect('auto')
    #    pylab.show()
        
        #alternate version
        fig = pylab.figure(19)
        pylab.clf()
        ax = pylab.subplot(1,1,1)
        testLattice.draw_resonator_lattice(ax, color = 'cornflowerblue', linewidth = 1)
        testLattice.draw_resonator_end_points(ax, color = 'indigo', edgecolor = 'indigo',  marker = '+' , size = 20)
        testLattice.plot_end_layout_state(state_uniform, ax, title = 'unit cell convention', colorbar = False, plot_links = False, cmap = 'jet', scaleFactor = 0.5)
        testLattice.draw_SDlinks(ax, linewidth = 1.5, extra = False, minus_links = True, minus_color = 'goldenrod')
        pylab.title('site orientations')
#        ax.set_aspect('auto')
        ax.set_aspect('equal')
    #    fig.savefig('HW.png', dpi = 200)
        pylab.show()
    
        #show lattice and medial
        fig = pylab.figure(20)
        pylab.clf()
        ax = pylab.subplot(1,1,1)
    #    testLattice.draw_resonator_lattice(ax, color = 'cornflowerblue', linewidth = 2)
        testLattice.draw_resonator_lattice(ax, color = 'firebrick', linewidth = 2)
        testLattice.draw_SDlinks(ax, linewidth = 2, extra = False, minus_links = False, color = 'goldenrod')
        pylab.title('site orientations')
        ax.set_aspect('auto')
#        ax.set_aspect('equal')
        ax.axis('off')
    #    fig.savefig('HL.png', dpi = 200)
        pylab.show()
    
    #    #show just the medial
    #    fig = pylab.figure(21)
    #    pylab.clf()
    #    ax = pylab.subplot(1,1,1)
    #    testLattice.draw_SDlinks(ax, linewidth = 1.5, extra = False, minus_links = False, color = 'mediumblue')
    ##    ax.set_aspect('auto')
    #    ax.set_aspect('equal')
    #    ax.axis('off')
    ##    fig.savefig('Kagome.png', dpi = 200)
    #    pylab.show()
        
        
        
    #        #show lattice and medial
    #    fig = pylab.figure(21)
    #    pylab.clf()
    #    ax = pylab.subplot(1,2,1)
    #    testLattice.draw_resonator_lattice(ax, color = 'firebrick', linewidth = 2)
    #    testLattice.draw_SDlinks(ax, linewidth = 2, extra = False, minus_links = False, color = 'goldenrod')
    #    pylab.title('original resonators')
    ##    ax.set_aspect('auto')
    #    ax.set_aspect('equal')
    #    ax.axis('off')
    #    
    #    ax = pylab.subplot(1,2,2)
    #    testLattice.draw_SD_points(ax, color = 'dodgerblue', edgecolor = 'k',  marker = 'o' , size = 10)
    #    pylab.title('SD sites')
    #    ax.set_aspect('equal')
    #    ax.axis('off')
    #    
    ##    fig.savefig('HL.png', dpi = 200)
        
        
    #    
    













