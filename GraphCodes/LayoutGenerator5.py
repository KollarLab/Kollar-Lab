#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 28 12:51:45 2017

@author: kollar2

v3 - adding square and possibly hexagonal planar lattices togain intuition
     Also fix a bug with the radial linkes at self.vertex>3
     
     switch to scipy.linalg.eigh for diagonalization. It's taylored to hermitian and symmetric 
     and gives much better result. Numpy.linalg.eig just plain messes up and gives
     non-othogonal eigenvectors
     
     adding a lot fo plot functions so I don't have to keep pasting the code

v4 - generalizing the tiling and link building rules. No longer have seperate functions for different tilings.
     Should now be able to handle every tiling (except for spherical ones)
     
     Also added more automatic plot functions.
     
     And a populate function the automatically fills everything up to a specified itteration
     
     Also come functions for calculating interactions
     
v5 - Added capability to deal with HW resonators properly

    Also modfied how SD points are generated. Previously zeorth azimuthal points were visibly
    offset from how the resoantors are drawn. Moved the SD points to lie on resonator chords
    rather than on the circle. Have not changed higher itterations because the difference
    is not visible there. 2-20-18
    
    7-25-18  Fixed bug in resonatorsOnly option. resonators only will now generate the semiduals
    They are needed for determining the sizes of things using get_all_resonators.
         Also added zorder optional argument to all the draw functions
     
     Methods:
        ###########
        #automated construction, saving, loading
        ##########
        populate
        save
        load
        
        ########
        #functions to generate the resonator lattice
        #######
        itter_generate
        _itter_generate_full_general
         
        #######
        #resonator lattice get /view functions
        #######
        get_xs
        get_ys
        get_radials
        draw_radials
        draw_azimuthals
        draw_all_azimuthals
        draw_all_radials
        draw_resonator_lattice
        draw_resonator_end_points
        get_all_resonators
        
        ########
        #functions to generate effective JC-Hubbard lattice (semiduals)
        ######## 
        generate_semiduals
        _azimuthal_links_full_general
        _radial_links_full_general
        _radial_links_non_triangle (defunct)
        
        #######
        #get and view functions for the JC-Hubbard (semi-dual lattice)
        #######
        draw_SDlinks
        get_semidual_points (semi-defunct)
        get_all_semidual_points
        
        ######
        #Hamiltonian related methods
        ######
        generate_Hamiltonian
        get_sub_Hamiltonian
        get_eigs
        
        ##########
        #methods for calculating/looking at states and interactions
        #########
        sh
        build_local_state_az
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
        from LayoutGenerator4 import PlanarLayout
        test = PlanarLayout(file_path = '7gon_3vertex_ 2.pkl')
        
        #####
        #making new layout
        #####
        from LayoutGenerator5 import PlanarLayout
        test = PlanarLayout(gon = 7, vertex = 3,side = 1, radius_method = 'lin', modeType = 'FW')
        test.populate(maxItter = 1)
        
        #####
        #saving computed layout
        #####
        test.save( name = 'filename.pkl') #filename can be a full path, but must have .pkl extension
        

   
     
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




class PlanarLayout(object):
    def __init__(self, gon = 4, vertex = 4,side = 1, radius_method = 'lin', file_path = '', modeType = 'FW'):
        if file_path != '':
            self.load(file_path)
        else:
            #create plank planar layout object with the bare bones that you can build on
            self.gon = gon
            self.vertex = vertex
            self.side = side*1.0
            self.alpha = 2.*numpy.pi/gon #inside angle
            self.theta = numpy.pi - self.alpha #turn angle at vertex
            self.radius = self.side/numpy.cos(self.theta/2)/2
            
            #starting plaquette
    #        self.x0 = numpy.asarray([numpy.sqrt(side),numpy.sqrt(side),-numpy.sqrt(side),-numpy.sqrt(side) ])
    #        self.y0 = numpy.asarray([numpy.sqrt(side),-numpy.sqrt(side),numpy.sqrt(side),-numpy.sqrt(side) ])
    #        self.edge_angles = scipy.arange(0, 2*numpy.pi, self.alpha) 
            self.edge_angles = scipy.arange(0, 2*numpy.pi, self.alpha) +self.alpha/2
    #        self.diag_angles = self.edge_angles + self.alpha/2
            self.diag_angles = self.edge_angles - self.alpha/2
            self.x0 = numpy.cos(self.diag_angles)*self.radius
            self.y0 = numpy.sin(self.diag_angles)*self.radius
            
            #lattice
            self.coords = numpy.zeros((len(self.x0),2)) #only the azimuthals stored here
            self.coords[:,0] = self.x0
            self.coords[:,1] = self.y0
            
            #dictionary to hold the points and the radials at each itteration of the generation
            self.radii = {}
            self.radials = {} #angular location of the radials
            self.points = {} #angular location of the points
            self.radials[0] = numpy.asarray([[], []]).reshape(0,2) #format = starting angle, ending angle
            self.points[0] = self.diag_angles
            self.radii[0] = self.radius
            self.itter = 0
            
            self.radius_method = radius_method
            
            if not ((modeType == 'FW') or (modeType  == 'HW')):
                raise ValueError('Invalid mode type. Must be FW or HW')
            self.modeType = modeType 
    
    ###########
    #automated construction, saving, loading
    ##########
    def populate(self, maxItter, resonatorsOnly = False, Hamiltonian = True, save = False, save_name = ''):
        '''
        fully populate the structure up to itteration = MaxItter
        
        if Hamiltonian = False will not generate H
        save is obvious
        '''
        if self.itter != 0:
            raise ValueError('This is not an empty planar layout')
         
        #make the resonator lattice
        for itt in range(1, maxItter+1):
            self.itter_generate()
        
        #make the JC-Hubbard lattice
        self.generate_semiduals() #need this so that various sizes are defined. Does not do big link computations
        if not resonatorsOnly:
#            #make the JC-Hubbard lattice
#            self.generate_semiduals()
            
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
            name = str(self.gon) + 'gon_' + str(self.vertex) + 'vertex_ ' + str(self.itter) + waveStr + '.pkl'
        
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
            print('Old pickle file. Pre FW-HW.')
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
    
    def get_radials(self, itteration):
        '''
        get the matrix of all radial connections at a given itteration
        
        radials are stored as a collumn matrix. First collumn in the angular coordinate
        of the starting point of the radial. Second collumn is the angular coordinate 
        of the end point
        '''
        radialAngles = self.radials[itteration]
        outerRadius = self.radii[itteration]
        if len(radialAngles) != 0:
            innerRadius = self.radii[itteration-1]
            
        rxs_in = innerRadius*numpy.cos(radialAngles[:,0])
        rys_in = innerRadius*numpy.sin(radialAngles[:,0])
        
        rxs_out = outerRadius*numpy.cos(radialAngles[:,1])
        rys_out = outerRadius*numpy.sin(radialAngles[:,1])
        return [rxs_in, rys_in, rxs_out, rys_out]
    
    def draw_radials(self, ax, itteration, color = 'g', alpha = 1, linewidth = 0.5, zorder = 1):
        '''
        draw all radial connections
        for a specific itteration
        '''
        if itteration == 0:
            return
        else:
            [rxs_in, rys_in, rxs_out, rys_out] = self.get_radials(itteration)
            for ind in range(0, len(rxs_out)):
                ax.plot([rxs_in[ind], rxs_out[ind]],[rys_in[ind], rys_out[ind]] , color = color, alpha = alpha, linewidth = linewidth, zorder = zorder)
                
    def draw_azimuthals(self, ax, itteration, color = 'g', alpha = 1, linewidth = 0.5, zorder = 1):
        '''
        for a specific itteration
        '''
        points = self.points[itteration]
        xs = self.radii[itteration]*numpy.cos(points)
        ys = self.radii[itteration]*numpy.sin(points)
        for ind in range(0, len(points)):
            x1 = xs[ind]
            y1 = ys[ind]
            x2 = xs[numpy.mod(ind+1, len(points))]
            y2 = ys[numpy.mod(ind+1, len(points))]
            ax.plot([x1, x2],[y1, y2] , color = color, alpha = alpha, linewidth = linewidth, zorder = zorder)
        return
    
    def draw_all_azimuthals(self, ax, maxItter = -1, color = 'g', alpha = 1 , linewidth = 0.5, zorder = 1):
        if maxItter > self.itter:
            raise ValueError('dont have this many itterations')
        elif maxItter <0:
            maxItter = self.itter
            
        for itteration in range(0, maxItter+1):
            self.draw_azimuthals(ax, itteration, color = color, alpha = alpha, linewidth = linewidth, zorder = zorder)
        return
                
    def draw_all_radials(self, ax, maxItter = -1, color = 'g', alpha = 1 , linewidth = 0.5, zorder = 1):
        if maxItter > self.itter:
            raise ValueError('dont have this many itterations')
        elif maxItter <0:
            maxItter = self.itter
            
        for itteration in range(1, maxItter + 1):
            self.draw_radials(ax, itteration, color = color, alpha = alpha, linewidth = linewidth, zorder = zorder)
        return

    def draw_resonator_lattice(self, ax, mode = 'line', maxItter = -1, color = 'g', alpha = 1 , linewidth = 0.5, zorder = 1):
        if maxItter > self.itter:
            raise ValueError('dont have this many itterations')
        elif maxItter <0:
            maxItter = self.itter
            
        if mode == 'line':
            self.draw_all_radials(ax, maxItter = maxItter, color = color, alpha = alpha, linewidth = linewidth, zorder = zorder)
            self.draw_all_azimuthals(ax, maxItter = maxItter, color = color, alpha = alpha, linewidth = linewidth, zorder = zorder)
        return
    
#    def draw_resonator_end_points(self, ax, maxItter = -1, color = 'g', alpha = 1 , linewidth = 0.5):
#        if maxItter > self.itter:
#            raise ValueError, 'dont have this many itterations'
#        elif maxItter <0:
#            maxItter = self.itter
#            
#        #plot the azimuthal points
#        xs = self.coords[:,0]
#        ys = self.coords[:,1]
#        pylab.sca(ax)
#        pylab.scatter(xs, ys ,c =  color, s = 10, marker = 'o', edgecolors = 'k')
#
#        return
    def draw_resonator_end_points(self, ax, maxItter = -1, color = 'g', edgecolor = 'k',  marker = 'o' , size = 10, zorder = 1):
        if maxItter > self.itter:
            raise ValueError('dont have this many itterations')
        elif maxItter <0:
            maxItter = self.itter
            
        #plot the azimuthal points
        xs = self.coords[:,0]
        ys = self.coords[:,1]
        pylab.sca(ax)
        pylab.scatter(xs, ys ,c =  color, s = size, marker = marker, edgecolors = edgecolor, zorder = zorder)

        return
    
    def get_all_resonators(self, maxItter = -1):
        '''
        function to get all resonators as a pair of end points
        
        each resontator returned as a row with four entries.
        (orientation is important to TB calculations)
        x0,y0,x1,y1
        
        '''
        if maxItter > self.itter:
            raise ValueError('dont have this many itterations')
        elif maxItter <0:
            maxItter = self.itter
        
        if maxItter < self.itter:
            size = self.get_SDindex(0,maxItter+1,az = True)
        else:
            size = len(self.SDx)
        resonators = numpy.zeros((size, 4))
        ind = 0
        for itteration in range(0, maxItter + 1):
            #get the azimuthal resonators
            points = self.points[itteration]
            xs = self.radii[itteration]*numpy.cos(points)
            ys = self.radii[itteration]*numpy.sin(points)
#            print itteration
#            print xs
#            print ys
            for az in range(0, len(points)):
                x1 = xs[az]
                y1 = ys[az]
                x2 = xs[numpy.mod(az+1, len(points))]
                y2 = ys[numpy.mod(az+1, len(points))]
                resonators[ind,:] = [x1,y1,x2,y2]
                ind = ind +1
            if itteration>=1:
                #get the radial resonators
#                [rxs_in, rys_in, rxs_out, rys_out] = self.get_radials(itteration)
                temp = numpy.asarray(self.get_radials(itteration))
                [cols, rows] = temp.shape
                resonators[ind: ind+rows, :] = numpy.transpose(temp)
                ind = ind + rows
                      
        return resonators
    
    ########
    #functions to generate the resonator lattice
    #######      
    def itter_generate(self):
        '''
        generate the next itteration of the resonator lattice
        '''
        self._itter_generate_full_general()

    
    def _itter_generate_full_general(self):
        '''
        full general (I hope) itteration function to generate the lattice
        by propagation of the tiling rule.
        
        Should be abe to handle and flat or hyperbolic tiling.
        Spherical tilings tend to crash eventudally because the tiling rule ends
        '''
        
        print('finished itteration equals = ' + str(self.itter))
        print('attempting itteration ' + str(self.itter+1))
        numVertices = len(self.points[self.itter])
        numFaces = len(self.points[self.itter])
        numRadials = self.radials[self.itter].shape[0]
        newItter = self.itter+1
        
        #find the new radial lines
        outgoing_radials_at_each_vertex = numpy.zeros(numVertices) #array to tell how manyrials emerge from each point
        for ind in range(0, numVertices):
            #see if this point is the end of a previous itterations radial
            currentAngle = numpy.mod(self.points[self.itter][ind], 2*numpy.pi)
            incoming = len(numpy.where(numpy.mod(self.radials[self.itter][:,1],2*numpy.pi)==currentAngle)[0])
            outgoing_radials_at_each_vertex[ind] = self.vertex-2-incoming
                
        print(outgoing_radials_at_each_vertex[0:7])
        print('\n')
        
        
        #find the number of new points
        currentIndexP = 0
        currentIndexR = 0
        for ind in range(0, numVertices):
            for rad in range(0, int(outgoing_radials_at_each_vertex[ind])):
                #do each gon that only has a corner touching the previous itteration (mostly)
                #move to the next site
                currentIndexP = currentIndexP + (self.gon-2)
            #back up one site, because the next gon touches the previous itteration at a whole face, so we need a smaller step
            currentIndexP = currentIndexP - 1
            
        numNewRadials = int(numpy.sum(outgoing_radials_at_each_vertex))
        numNewPoints = int(currentIndexP)
#        print numNewRadials
#        print numNewPoints
        
        
        #find the l ocations of the new azimuthal vertices
#        self.points[newItter ] = 0 + scipy.arange(0,2*numpy.pi, 2*numpy.pi/numNewPoints) #don't want 2pi included
#        self.points[newItter ] = -self.alpha/2 + scipy.arange(0,2*numpy.pi, 2*numpy.pi/numNewPoints)
        if self.gon == 7:
#            self.points[newItter ] = scipy.arange(0,2*numpy.pi, 2*numpy.pi/numNewPoints) - 0*self.alpha/2 + 0.25*(self.itter)**0.5*self.alpha/2.
            self.points[newItter ] = scipy.arange(0,numNewPoints,1)*2*numpy.pi/numNewPoints - 0*self.alpha/2 + 0.25*(self.itter)**0.5*self.alpha/2.
        elif self.gon == 5:
            if self.itter == 0:
#                self.points[newItter ] = scipy.arange(0,2*numpy.pi, 2*numpy.pi/numNewPoints) + 0.
                self.points[newItter ] = scipy.arange(0,numNewPoints,1)*2*numpy.pi/numNewPoints + 0.
            else:
#                self.points[newItter ] = scipy.arange(0,2*numpy.pi, 2*numpy.pi/numNewPoints) + 1.0*self.alpha/2.
#                self.points[newItter ] = scipy.arange(0,2*numpy.pi, 2*numpy.pi/numNewPoints) - 0*self.alpha/2 + 0.25*(self.itter)**0.5*self.alpha/2.
                self.points[newItter ] = scipy.arange(0,numNewPoints,1)*2*numpy.pi/numNewPoints - 0*self.alpha/2 + 0.25*(self.itter)**0.5*self.alpha/2.
        elif self.gon == 3:  
#            self.points[newItter ] =  scipy.arange(0,2*numpy.pi, 2*numpy.pi/numNewPoints) - self.alpha/2/1.5 #old as of 5-23-18
            self.points[newItter ] =  scipy.arange(0,numNewPoints,1)*2*numpy.pi/numNewPoints - self.alpha/2/1.5 #old as of 5-23-18
        else:
#            self.points[newItter ] =  scipy.arange(0,2*numpy.pi, 2*numpy.pi/numNewPoints) - self.alpha/2/1.5 #old as of 5-23-18
#            self.points[newItter ] = scipy.arange(0,2*numpy.pi, 2*numpy.pi/numNewPoints) - 0*self.alpha/2 + 0.25*(self.itter)**0.5*self.alpha/2.
            self.points[newItter ] = scipy.arange(0,numNewPoints,1)*2*numpy.pi/numNewPoints - 0*self.alpha/2 + 0.25*(self.itter)**0.5*self.alpha/2.
        if self.radius_method == 'lin':
            self.radii[newItter ] = self.radii[self.itter] + self.radius
        if self.radius_method == 'exp':
            self.radii[newItter ] = self.radius*numpy.exp(self.itter+1)
        if self.radius_method == 'Mattias':
#            MattiasRadii = [1.0, 3, 4.5, 7.,8.]
#            MattiasRadii = [1.0, 3.1, 4.8, 7.,8.]
            MattiasRadii = [1.0, 3.2, 5.1, 7.,8.]
            self.radii[newItter ] = self.radius*MattiasRadii[newItter]
        
        #store the cartesian locations of the new vertices
        xs = self.radii[newItter]*numpy.cos(self.points[newItter])
        ys = self.radii[newItter]*numpy.sin(self.points[newItter])
        newCoords = numpy.zeros((len(xs),2))
        newCoords[:,0] = xs
        newCoords[:,1] = ys
        self.coords = numpy.concatenate((self.coords, newCoords))
        
        
        self.radials[newItter] = numpy.zeros((numNewRadials,2))
        currentIndexP = 0
        currentIndexR = 0
        for ind in range(0, numVertices):
            for rad in range(0, int(outgoing_radials_at_each_vertex[ind])):
#                print ind
#                print currentIndexR
##                print currentIndexP
#                print '\n'
                #do each gon that only has a corner touching the previous itteration (mostly)
                self.radials[newItter][currentIndexR, :] = numpy.asarray([self.points[self.itter][ind], self.points[newItter][currentIndexP]])
                #ove to the next site
                currentIndexP = numpy.mod(currentIndexP + (self.gon-2), numNewPoints)
#                currentIndexP = currentIndexP + (self.gon-2)
                currentIndexR = currentIndexR + 1
            
            #back up one site, because the next gon touches the previous itteration at a whole face, so we need a smaller step
#            currentIndexP = numpy.mod(currentIndexP - 1, numNewPoints)
            if ind == 0 and int(outgoing_radials_at_each_vertex[ind]) ==0:
                #haven't actually stepped, so I don't want to remove anything
                pass
            else:
                currentIndexP = numpy.mod(currentIndexP - 1, numNewPoints)
                

        
        self.itter = self.itter+1
        return
    
    ########
    #functions to generate effective JC-Hubbard lattice
    ########    
    def generate_semiduals(self):
        '''
        main workhorse function to generate the JC-Hubbard lattice.
        This is the one you shold call. All the others are workhorses that it uses.
        
        Will loop through the existing itterations and create attributes for the 
        JC-Hubbard lattice (here jokingly called semi-dual) and fill them
        '''
        maxItter = self.itter
        
        self.SDradials = {}
        self.SDpoints = {}
        self.SDradii = {}  #angular location of the midpoints on the radials. I don't think so
        
        self.SDcartesian = {} #cartesian location of all the semidual nodes. To be filled later
        self.SDlinks = {} #4 entries: itteration of first point, index of first point, itteration of second, index of second
        self.SDHWlinks = {} #same as above, but 6 entries to keep track of each end of resonator for half wave resonators.
                            #5th entry is end of sourse, 6th entry is end of target
        
        #midline points on the circles
        for i in range(0, maxItter+1):
            self.SDpoints[i] = self.points[i] + (self.points[i][1] - self.points[i][0])/2.
            
        
            
        #radii of the midline points on the radials
        self.SDradii[0] = numpy.asarray([])
        for i in range(1, maxItter+1):
            self.SDradii[i] = (self.radii[i] + self.radii[i-1])/2.
            
        #angles of the midline points on the radials
        self.SDradials[0] = numpy.asarray([])
        for i in range(1, maxItter+1):
            x0s = self.radii[i-1]*numpy.cos(self.radials[i][:,0])
            y0s = self.radii[i-1]*numpy.sin(self.radials[i][:,0])
            
            x1s = self.radii[i]*numpy.cos(self.radials[i][:,1])
            y1s = self.radii[i]*numpy.sin(self.radials[i][:,1])
            
            xs = (x0s + x1s)/2.
            ys = (y0s + y1s)/2.
            
            thetas = numpy.arctan2(ys,xs)
            self.SDradials[i] = thetas 
            
        #get cartesian coordinates for everything
        for i in range(0, maxItter+1):
            [xs, ys] = self.get_semidual_points(i) #retreive an array with them all smooshed together
            self.SDcartesian[i] = numpy.zeros((len(xs), 2))
            self.SDcartesian[i][:,0] = xs
            self.SDcartesian[i][:,1] = ys
            
        #generate the links         
        self._azimuthal_links_full_general(maxItter)
        
        #the trickier radial links 
        for i in range(1,maxItter+1):
            self._radial_links_full_general(i)
#            if self.gon > 3:
##                self._radial_links_non_triangle(i)
#                self._radial_links_triangle(i)
#            else:
#                self._radial_links_triangle(i)
                
        self.get_all_semidual_points()
        return
    
    def _azimuthal_links_full_general(self, maxItter):
        '''
        calculate and store all links between points on Azimuthal rings
        '''
        self.numAzLinks = {} #I need to store the number of points skipped

        #loop over each itteration
        for itter in range(0, maxItter+1):
            numAz = len(self.points[itter]) #number azimuthal points at this itteration
            numRad = self.radials[itter].shape[0]
                    
            #find the new radial lines
            outgoing_radials_at_each_vertex = numpy.zeros(numAz) #array to tell how manyrials emerge from each point
            for ind in range(0, numAz):
                #see if this point is the end of a previous itterations radial
                currentAngle = numpy.mod(self.points[itter][ind], 2*numpy.pi)
                incoming = len(numpy.where(numpy.mod(self.radials[itter][:,1],2*numpy.pi)==currentAngle)[0])
                outgoing_radials_at_each_vertex[ind] = self.vertex-2-incoming
            incoming_radials_at_each_vertex  = self.vertex - outgoing_radials_at_each_vertex - 2
#            print outgoing_radials_at_each_vertex
#            print incoming_radials_at_each_vertex                  

            #will skip any point where there is an incoming radial (these are cyclic coouplers. No diagonal terms for now)        
#            skipped_points = len(numpy.where(incoming_radials_at_each_vertex != 0)[0]) 
            skipped_maybe1 = numpy.where(incoming_radials_at_each_vertex != 0)[0] 
            skipped_maybe2 = numpy.where(outgoing_radials_at_each_vertex != 0)[0]
            skipped_sites= numpy.intersect1d(skipped_maybe1, skipped_maybe2)
            skipped_points = len(skipped_sites)
            self.numAzLinks[itter] = numAz - skipped_points
            
            #allocate space
#            self.SDlinks[itter] = numpy.zeros(((numAz-skipped_points)*2+ numRad*8, 4))
            self.SDlinks[itter] = numpy.zeros(((numAz-skipped_points)*2+ numRad*10, 4)) #temporary
            self.SDHWlinks[itter] = numpy.zeros(((numAz-skipped_points)*2+ numRad*10, 6))
            
            #fill in azimuthal links
            currInd = 0
            for vertex in range(0, numAz):
                
                next_vertex = numpy.mod(vertex+1, numAz)
#                if incoming_radials_at_each_vertex[vertex] !=0 and incoming_radials_at_each_vertex[vertex] !=0:
                if incoming_radials_at_each_vertex[vertex] !=0 and outgoing_radials_at_each_vertex[vertex] !=0:
                    pass
                else:
                    #-1 takes care of how the links are indexed compared to the original vertices
                    self.SDlinks[itter][currInd,:] = [itter, numpy.mod(vertex-1, numAz), itter, numpy.mod(next_vertex-1, numAz)]
                    self.SDlinks[itter][currInd+1,:] = [itter, numpy.mod(next_vertex-1, numAz), itter, numpy.mod(vertex-1, numAz)]
                    self.SDHWlinks[itter][currInd,:] = [itter, numpy.mod(vertex-1, numAz), itter, numpy.mod(next_vertex-1, numAz), 1,0] #0 or 1 is for each end of the resonator
                    self.SDHWlinks[itter][currInd+1,:] = [itter, numpy.mod(next_vertex-1, numAz), itter, numpy.mod(vertex-1, numAz),0,1]
                    currInd  = currInd + 2
            
        return
    
    def _radial_links_full_general(self, itteration):
        '''
        calculate and store all links involving points radials
        '''
        print('radial link itteration ' + str(itteration))
        i = itteration
        
        #now add the radials
        numAz = len(self.points[i])
        numAz_minus1 = len(self.points[i-1])
        numRad = self.radials[i].shape[0]
        numRad_minus1 = self.radials[i-1].shape[0]
#        blankRads  = numpy.zeros((numRad*8, 4))
        blankRads  = numpy.zeros((numRad*10, 4))   # temporary
        blankRadsHW  = numpy.zeros((numRad*10, 6))
        
        #find the outgoing radial lines from the previous itteration
        outgoing_radials_at_each_vertex = numpy.zeros(numAz_minus1) #array to tell how manyrials emerge from each point
        for ind in range(0, numAz_minus1):
            #see if this point is the end of a previous itterations radial
            currentAngle = numpy.mod(self.points[i-1][ind], 2*numpy.pi)
            incoming = len(numpy.where(numpy.mod(self.radials[i-1][:,1],2*numpy.pi)==currentAngle)[0])
            outgoing_radials_at_each_vertex[ind] = self.vertex-2-incoming
            
        #find a good starting point in the radials, because the first one might be in the middle of a cluster
        currentIndexR = 0
        for spoke in range(0, numRad_minus1):
            startInd = numpy.where(self.radials[i][numpy.mod(currentIndexR, numRad),0] == self.points[i-1])[0][0]
            startInd_previous = numpy.where(self.radials[i][numpy.mod(currentIndexR-1, numRad),0] == self.points[i-1])[0][0]
            if startInd == startInd_previous:
                #current radial is in the middle of a converging cluster.
                #step back one radial to try to get closer to the start of the cluster
                currentIndexR = numpy.mod(currentIndexR-1, numRad)
            else:
                #we've reached the end of the cluster
                break
        
        #do the links to the previous azimuthal ring
        currInd = 0 #index for the rows of the link matrix
        print('starting radial  = ' + str(currentIndexR)  + ' of ' + str(numRad))
        print('outgoing radials: '+ str(outgoing_radials_at_each_vertex[0:7]))
        for vertex in range(0,numAz_minus1):
            if outgoing_radials_at_each_vertex[vertex] == 0:
                pass
                #no links to be made
            else:
                #there are outgoing radial links, so I have to do something
                
                #loop over the outgoing radials
                for rad in range(0, int(outgoing_radials_at_each_vertex[vertex])):
#                    print currentIndexR
#                    print 'rad = ' + str(rad)
                    startInd = numpy.where(self.radials[i][currentIndexR,0] == self.points[i-1])[0][0]
                    stopInd = numpy.where(self.radials[i][currentIndexR,1] == self.points[i])[0][0]
                    if rad == 0:
#                        print 'first radial'
                        #first outgoing radial. 
                        #connect to the previous azimuthal on theprevious ring
                        blankRads[currInd,:] = numpy.asarray([i-1, numpy.mod(startInd-1, numAz_minus1), i, numAz + currentIndexR])
                        blankRadsHW[currInd,:] = numpy.asarray([i-1, numpy.mod(startInd-1, numAz_minus1), i, numAz + currentIndexR, 1, 0])
                        currInd = currInd + 1
                        
                        #### reverse links
                        blankRads[currInd,:] = numpy.asarray([i, numAz + currentIndexR, i-1, numpy.mod(startInd-1, numAz_minus1)])
                        blankRadsHW[currInd,:] = numpy.asarray([i, numAz + currentIndexR, i-1, numpy.mod(startInd-1, numAz_minus1), 0,1])
                        currInd = currInd + 1
                        
                    else:
#                        print 'not first radial'
                        #not a first radial
                        #connect to the previous radial
                        blankRads[currInd,:] = numpy.asarray([i, numAz + numpy.mod(currentIndexR-1, numRad), i,numAz + currentIndexR])
                        blankRadsHW[currInd,:] = numpy.asarray([i, numAz + numpy.mod(currentIndexR-1, numRad), i,numAz + currentIndexR,0,0])
                        currInd = currInd + 1
                        #### reverse link
                        blankRads[currInd,:] = numpy.asarray([i,numAz + currentIndexR, i, numAz + numpy.mod(currentIndexR-1, numRad)])
                        blankRadsHW[currInd,:] = numpy.asarray([i,numAz + currentIndexR, i, numAz + numpy.mod(currentIndexR-1, numRad),0,0])
                        currInd = currInd + 1
                     
                        
                    if rad +1 == (int(outgoing_radials_at_each_vertex[vertex])):
#                        print 'last radial'
                        #last radial
                        #connect to the next azimuthal on the previousring (not radials effectively start from half-integer sites )
                        blankRads[currInd,:] = numpy.asarray([i-1, numpy.mod(startInd, numAz_minus1), i, numAz + currentIndexR])
                        blankRadsHW[currInd,:] = numpy.asarray([i-1, numpy.mod(startInd, numAz_minus1), i, numAz + currentIndexR,0,0])
                        currInd = currInd + 1
                        
                        #### reverse links
                        blankRads[currInd,:] = numpy.asarray([i, numAz + currentIndexR, i-1, numpy.mod(startInd, numAz_minus1)])
                        blankRadsHW[currInd,:] = numpy.asarray([i, numAz + currentIndexR, i-1, numpy.mod(startInd, numAz_minus1),0,0])
                        currInd = currInd + 1

                    currentIndexR = numpy.mod(currentIndexR +1, numRad)
                    
                    
        #find the incoming radial lines from the current itteration
        incoming_radials_at_each_vertex = numpy.zeros(numAz)
        for ind in range(0, numAz):
            #see if this point is the end of a previous itterations radial
            currentAngle = numpy.mod(self.points[i][ind], 2*numpy.pi)
            incoming = len(numpy.where(numpy.mod(self.radials[i][:,1],2*numpy.pi)==currentAngle)[0])
            incoming_radials_at_each_vertex[ind] = incoming
            
        #find a good starting point in the radials, because the first one might be in the middle of a cluster
        currentIndexR = 0
        for spoke in range(0, numRad):
            stopInd = numpy.where(self.radials[i][numpy.mod(currentIndexR, numRad),1] == self.points[i])[0][0]
            stopInd_previous = numpy.where(self.radials[i][numpy.mod(currentIndexR-1, numRad),1] == self.points[i])[0][0]
            if stopInd == stopInd_previous:
                #current radial is in the middle of a converging cluster.
                #step back one radial to try to get closer to the start of the cluster
                currentIndexR = numpy.mod(currentIndexR-1, numRad)
            else:
                #we've reached the end of the cluster
                break
        
        #do the forward links to the current azimuthal ring
        print('starting radial  = ' + str(currentIndexR) + ' of ' + str(numRad))
        print('incoming radials: '+ str(incoming_radials_at_each_vertex[0:7]))
        for vertex in range(0,numAz):
            if incoming_radials_at_each_vertex[vertex] == 0:
                pass
                #no links to be made
            else:
                #there are incoming radial links, so I have to do something
                
                #loop over the incoming radials
#                print int(incoming_radials_at_each_vertex[vertex])
                for rad in range(0, int(incoming_radials_at_each_vertex[vertex])):
#                    print currentIndexR
#                    print 'rad = ' + str(rad)
                    startInd = numpy.where(self.radials[i][currentIndexR,0] == self.points[i-1])[0][0]
                    stopInd = numpy.where(self.radials[i][currentIndexR,1] == self.points[i])[0][0]
                    if rad == 0:
                        #print 'first incoming radial'
                        #first outgoing radial. 
                        #connect to the previous azimuthal on the next ring
                        blankRads[currInd,:] = numpy.asarray([i, numpy.mod(stopInd-1, numAz), i, numAz + currentIndexR])
                        blankRadsHW[currInd,:] = numpy.asarray([i, numpy.mod(stopInd-1, numAz), i, numAz + currentIndexR,1,1])
                        currInd = currInd + 1
                        
                        #### reverse links
                        blankRads[currInd,:] = numpy.asarray([i, numAz + currentIndexR, i, numpy.mod(stopInd-1, numAz)])
                        blankRadsHW[currInd,:] = numpy.asarray([i, numAz + currentIndexR, i, numpy.mod(stopInd-1, numAz),1,1])
                        currInd = currInd + 1
                        
                    else:
                        #print 'not first incoming radial'
                        #not a first radial
                        #connect to the previous radial
                        blankRads[currInd,:] = numpy.asarray([i, numAz + numpy.mod(currentIndexR-1, numRad), i,numAz + currentIndexR])
                        blankRadsHW[currInd,:] = numpy.asarray([i, numAz + numpy.mod(currentIndexR-1, numRad), i,numAz + currentIndexR,1,1])
                        currInd = currInd + 1
                        #### reverse link
                        blankRads[currInd,:] = numpy.asarray([i,numAz + currentIndexR, i, numAz + numpy.mod(currentIndexR-1, numRad)])
                        blankRadsHW[currInd,:] = numpy.asarray([i,numAz + currentIndexR, i, numAz + numpy.mod(currentIndexR-1, numRad),1,1])
                        currInd = currInd + 1
                     
                        
                    if rad +1 == (int(incoming_radials_at_each_vertex[vertex])):
                        #print 'last incoming radial'
                        #last radial
                        #connect to the next azimuthal on the current ring (not radials effectively start from half-integer sites )
                        blankRads[currInd,:] = numpy.asarray([i, numpy.mod(stopInd, numAz), i, numAz + currentIndexR])
                        blankRadsHW[currInd,:] = numpy.asarray([i, numpy.mod(stopInd, numAz), i, numAz + currentIndexR,0,1])
                        currInd = currInd + 1
                        
                        #### reverse links
                        blankRads[currInd,:] = numpy.asarray([i, numAz + currentIndexR, i, numpy.mod(stopInd, numAz)])
                        blankRadsHW[currInd,:] = numpy.asarray([i, numAz + currentIndexR, i, numpy.mod(stopInd, numAz), 1,0])
                        currInd = currInd + 1

                    #step to the next radial
                    currentIndexR = numpy.mod(currentIndexR +1, numRad)

        
        #put it all together                
        self.SDlinks[i][(self.numAzLinks[i]*2):, :] = blankRads
        self.SDHWlinks[i][(self.numAzLinks[i]*2):, :] = blankRadsHW
        
        #cut out the extra zero rows
        self.SDlinks[i] = self.SDlinks[i][~numpy.all(self.SDlinks[i] == 0, axis=1)]    
        self.SDHWlinks[i] = self.SDHWlinks[i][~numpy.all(self.SDHWlinks[i] == 0, axis=1)]  

        return
    
    def _radial_links_non_triangle(self, itteration):
        '''
        vestigial function that can do radial links for all tilings except triangular.
        Triangular has the difficult property that ridials can converge at their ends.
        This property is not handled in the function.
        
        As of 7-20-17 this function is no longer in use.
        '''
        print(itteration)
        i = itteration
        
        #now add the radials
        numAz = len(self.points[i])
        numAz_minus1 = len(self.points[i-1])
        numRad = self.radials[i].shape[0]
#        blankRads  = numpy.zeros((numRad*8, 4))
        blankRads  = numpy.zeros((numRad*10, 4))   # temporary
        
        #find the outgoing radial lines from the previous itteration
        outgoing_radials_at_each_vertex = numpy.zeros(numAz_minus1) #array to tell how manyrials emerge from each point
        for ind in range(0, numAz_minus1):
            #see if this point is the end of a previous itterations radial
            currentAngle = numpy.mod(self.points[i-1][ind], 2*numpy.pi)
            incoming = len(numpy.where(numpy.mod(self.radials[i-1][:,1],2*numpy.pi)==currentAngle)[0])
            outgoing_radials_at_each_vertex[ind] = self.vertex-2-incoming
        
        currInd = 0
        currentIndexR = 0
        print(outgoing_radials_at_each_vertex[0:7])
        for vertex in range(0,numAz_minus1):
            if outgoing_radials_at_each_vertex[vertex] == 0:
                pass
                #no links to be made
            else:
                #there are outgoing radial links, so I have to do something
                
                #loop over the incoming radials
#                print int(outgoing_radials_at_each_vertex[vertex])
                for rad in range(0, int(outgoing_radials_at_each_vertex[vertex])):
#                    print currentIndexR
#                    print 'rad = ' + str(rad)
                    startInd = numpy.where(self.radials[i][currentIndexR,0] == self.points[i-1])[0][0]
                    stopInd = numpy.where(self.radials[i][currentIndexR,1] == self.points[i])[0][0]
                    if rad == 0:
#                        print 'first radial'
                        #first outgoing radial. 
                        #connect to the previous azimuthal on the previous ring
                        blankRads[currInd,:] = numpy.asarray([i-1, numpy.mod(startInd-1, numAz_minus1), i, numAz + currentIndexR])
                        currInd = currInd + 1
                        #connect to the previous azimuthal on current ring
                        blankRads[currInd,:] = numpy.asarray([i, numpy.mod(stopInd-1, numAz), i, numAz + currentIndexR])
                        currInd = currInd + 1
                        
                        #### reverse links
                        blankRads[currInd,:] = numpy.asarray([i, numAz + currentIndexR, i-1, numpy.mod(startInd-1, numAz_minus1)])
                        currInd = currInd + 1
                        blankRads[currInd,:] = numpy.asarray([ i, numAz + currentIndexR, i, numpy.mod(stopInd-1, numAz)])
                        currInd = currInd + 1
                        
                    else:
#                        print 'not first radial'
                        #not a first radial
                        #connect to the previous radial
                        blankRads[currInd,:] = numpy.asarray([i, numAz + numpy.mod(currentIndexR-1, numRad), i,numAz + currentIndexR])
                        currInd = currInd + 1
                        #### reverse link
                        blankRads[currInd,:] = numpy.asarray([i,numAz + currentIndexR, i, numAz + numpy.mod(currentIndexR-1, numRad)])
                        currInd = currInd + 1
                     
                    ######
                    #connections that should always me done
                    #####
                    #connect to the next azimuthal on current ring
                    blankRads[currInd,:] = numpy.asarray([i, numpy.mod(stopInd, numAz), i, numAz + currentIndexR])
                    currInd = currInd + 1
                    #connect to the previous azimuthal on current ring
                    blankRads[currInd,:] = numpy.asarray([i, numpy.mod(stopInd-1, numAz), i, numAz + currentIndexR])
                    currInd = currInd + 1
                    
                    #### reverse
                    blankRads[currInd,:] = numpy.asarray([i, numAz + currentIndexR, i, numpy.mod(stopInd, numAz)])
                    currInd = currInd + 1
                    #connect to the previous azimuthal on current ring
                    blankRads[currInd,:] = numpy.asarray([i, numAz + currentIndexR, i, numpy.mod(stopInd-1, numAz)])
                    currInd = currInd + 1
                        
                    if rad +1 == (int(outgoing_radials_at_each_vertex[vertex])):
#                        print 'last radial'
                        #last radial
                        #connect to the next azimuthal on the previousring (not radials effectively start from half-integer sites )
                        blankRads[currInd,:] = numpy.asarray([i-1, numpy.mod(startInd, numAz_minus1), i, numAz + currentIndexR])
                        currInd = currInd + 1
                        
                        #### reverse links
                        blankRads[currInd,:] = numpy.asarray([i, numAz + currentIndexR, i-1, numpy.mod(startInd, numAz_minus1)])
                        currInd = currInd + 1

                    currentIndexR = currentIndexR +1
                        
        self.SDlinks[i][(self.numAzLinks[i]*2):, :] = blankRads
        
        #cut out the extra zero rows
        self.SDlinks[i] = self.SDlinks[i][~numpy.all(self.SDlinks[i] == 0, axis=1)]    

        return
    
    #######
    #get and view functions for the JC-Hubbard (semi-dual lattice)
    #######
    def draw_SDlinks(self, ax, itteration = 2, color = 'firebrick', linewidth = 0.5, minus_links = False, minus_color = 'goldenrod', alpha = 1, zorder = 1):
        '''
        draw all the links of the semidual lattice
        '''
        #'lightskyblue'
        #'firebrick'
        if itteration == -1:
            return
        else:
            if self.itter == 1:
                itteration = 1
            for itt in range(0, itteration+1):
                for link in range(0, self.SDlinks[itt].shape[0]):
                    [startItt, startInd, stopItt, stopInd]  = self.SDlinks[itt][link,:]
                    startInd = int(startInd)
                    stopInd = int(stopInd)
                    startItt = int(startItt)
                    stopItt = int(stopItt)
                    
                    [x0, y0] = self.SDcartesian[startItt][startInd,:] 
                    [x1, y1] = self.SDcartesian[stopItt][stopInd,:] 
                
                    if  minus_links == True and self.modeType == 'HW':
                        ends = self.SDHWlinks[itt][link,4:]
                        if ends[0]==ends[1]:
                            #++ or --, use normal t
                            ax.plot([x0, x1],[y0, y1] , color = color, linewidth = linewidth, alpha = alpha, zorder = zorder)
                        else:
                            #+- or -+, use inverted t
                            ax.plot([x0, x1],[y0, y1] , color = minus_color, linewidth = linewidth, alpha = alpha, zorder = zorder)
                    else :
                        ax.plot([x0, x1],[y0, y1] , color = color, linewidth = linewidth, alpha = alpha, zorder = zorder)
    
    def get_semidual_points(self, itteration):
        '''
        get all the semidual points in a given itteration.
        
        Mostly just a workhorse for get_all_semidual_points
        
        Unlike Euclidean, resoantors aren't stored, so it doesn't put the points on the midpoints of the resonators.
        It uses an indpendnent cylindrical construction. For azimuthal rings past the zeroth you can't really 
        tell the difference by eye, but the zeroth looks bad.
        
        Will fix it for zeroth ring
        '''
        
        xs = numpy.zeros(len(self.points[itteration]) + self.radials[itteration].shape[0])
        ys = numpy.zeros(len(self.points[itteration]) + self.radials[itteration].shape[0])
        
        xs[0:len(self.points[itteration])] = self.radii[itteration]*numpy.cos(self.SDpoints[itteration])
        ys[0:len(self.points[itteration])] = self.radii[itteration]*numpy.sin(self.SDpoints[itteration])
        
        xs[len(self.points[itteration]):] = self.SDradii[itteration]*numpy.cos(self.SDradials[itteration])
        ys[len(self.points[itteration]):] = self.SDradii[itteration]*numpy.sin(self.SDradials[itteration])
        
        #fix the overlay for zeroth itteration
        if itteration == 0:
            zeroRes = self.get_all_resonators(maxItter = 0)
            #take average of end points
            tempXs = (zeroRes[:,0] + zeroRes[:,2])/2.
            tempYs = (zeroRes[:,1] + zeroRes[:,3])/2.
            
            xs[0:len(self.points[itteration])] = tempXs
            ys[0:len(self.points[itteration])] = tempYs
            
        
        return [xs, ys]
    
    def get_all_semidual_points(self, maxItter = -1):
        '''
        return x and y coordinates of all points in the semidual lattice
        (now dominantly a workhorse for generate_semiduals)
        
        As of 7-20-17, this function is falling out of general use because between v3 and a revision of v4
        it stated to store the coordinates it returns as the attributes SDx and SDy,
        and get_all_semidual_points is called at the end of generate_semiduals(), so there
        is no need to call it again and recalculate.
        '''
        if maxItter > self.itter:
            raise ValueError('dont have this many itterations')
        elif maxItter <0:
            maxItter = self.itter
            
        xs = numpy.asarray([])
        ys = numpy.asarray([])
        for i in range(0, maxItter+1):
            [tempx, tempy] = self.get_semidual_points(i)
            xs = numpy.concatenate((xs, tempx))
            ys = numpy.concatenate((ys, tempy))
            
        self.totalSDsites = len(xs)
        self.SDx = xs
        self.SDy = ys
        return [xs,ys]
    
    ######
    #Hamiltonian related methods
    ######
    def generate_Hamiltonian(self, t = 1, maxItter = -1, internalBond = 1000):
        '''
        create the effective tight-binding Hamiltonian
        
        if maxItter< self.itter, it will return the Hamiltonian for 
        a lattice of size maxItter.
        
        Also calculated as stores eigenvectors and eigenvalues for that H
        
        (calling this for maxItter<self.itte will store the eigenvectors and 
        eigenvalues of the subsystem.) Beware.
        
        Will use FW or HW TB coefficients depending on self.modeType
        '''
        
        if maxItter > self.itter:
            raise ValueError('dont have this many itterations')
        elif maxItter <0:
            maxItter = self.itter
        
        self.t = t
        self.internalBond = 1000*self.t
        
        totalSize = 0
        self.startInds = numpy.zeros(self.itter+1).astype('int') #starting index for the points in the ith itteration
        for itt in range(0, self.itter+1):
            self.startInds[itt] = totalSize
            totalSize = totalSize + self.SDcartesian[itt].shape[0]
            
        self.H = numpy.zeros((totalSize, totalSize))
        self.H_HW = numpy.zeros((totalSize*2, totalSize*2))
        
        #loop over the links and fill the Hamiltonian
        for itt in range(0, maxItter+1):
            for link in range(0, self.SDlinks[itt].shape[0]):
                [sourceItt, sourceInd, targetItt, targetInd] = self.SDlinks[itt][link, :]
                [sourceEnd, targetEnd] = self.SDHWlinks[itt][link, 4:]
                sourceEnd = int(sourceEnd)
                targetEnd = int(targetEnd)
                
                source = int(self.startInds[int(sourceItt)] + sourceInd)
                target = int(self.startInds[int(targetItt)] + targetInd)
                
                if self.modeType == 'FW':
                    self.H[source, target] = self.t
                elif self.modeType == 'HW':
                    polarity = sourceEnd^targetEnd #xor of the two ends. Will be one when the two ends are different
                    signum =(-1.)**(polarity) #will be zero when  two ends are same, and minus 1 otherwise
                    self.H[source, target] = self.t * signum
                else:
                    raise ValueError('You screwed around with the mode type. It must be FW or HW.')
                self.H_HW[2*source + sourceEnd, 2*target+targetEnd] = 2*self.t
                
        #fix the bonds between the two ends of the same site
        for site in range(0, totalSize):
            self.H_HW[2*site, 2*site+1] = self.internalBond
            self.H_HW[2*site+1, 2*site] = self.internalBond
                
        self.Es, self.Psis = scipy.linalg.eigh(self.H)
        self.Eorder = numpy.argsort(self.Es)
        
        return
    
    def get_sub_Hamiltonian(self, itteration):
        '''
        return the Hamiltonian for a subsystem (usually of smaller itteration numer)
        
        Does not mess with the full H or eigenvectors/eigenvlaues
        '''
        if itteration < self.itter:
            endPoint= int(self.startInds[itteration+1])
            return self.H[0:endPoint, 0:endPoint]
        else:        
            return self.H
         
    def get_eigs(self, maxItter=-1):
        '''
        returns eigenvectors and eigenvalues for a subsystem (usually of smaller itteration numer)
        
        Does not mess with the full H or eigenvectors/eigenvlaues
        '''
        if maxItter > self.itter:
            raise ValueError('dont have this many itterations')
        elif maxItter <0:
            maxItter = self.itter
            
        Hmat = self.get_sub_Hamiltonian(maxItter)
        
#        Es, Psis = numpy.linalg.eig(Hmat)
        Es, Psis = scipy.linalg.eigh(Hmat)
        Eorder = numpy.argsort(Es)
        
        return [Es, Psis, Eorder]
    
    ##########
    #methods for looking at states and interactions (plus a couple helper/constructor functions)
    #########    
    def get_SDindex(self,num, itt, az = True):
        '''
        get the index location of a semidual point. Point spcified by
        
        number within the itteration points/radials
        itteration number
        az (+True for azimuthal or false for radial)
        
        (useful for making localized states at specific sites)
        '''
        if itt == 0 and az == False:
            raise ValueError('no radials in zeroth itteration')
        
        currInd = 0
        for i in range(0, itt):
            currInd = currInd + len(self.SDpoints[i]) + self.SDradials[i].shape[0]
            
        if az == False:
            currInd = currInd + len(self.SDpoints[itt])
    
        currInd = currInd + num
        
        return currInd
    
    def build_local_state_az(self, site, maxItter= -1):
        '''
        build a state with only one site occupied on an azimuthal ring
        
        site is the number of the site in the maxItter-th azimuthal_ring
        
        if maxItter is negative (or left alone) will work with self.itter
        and automatically populate the site-th site on the outermost azimuthal ring
        '''
        if maxItter > self.itter:
                raise ValueError('dont have this many itterations')
        elif maxItter <0:
            maxItter = self.itter

        state = numpy.zeros(len(self.SDx))*(0+0j)
            
        currentInd = 0
        for itteration in range(0, maxItter+1):
            numPoints = len(self.points[itteration])
            numRadials = self.radials[itteration].shape[0]
            
            #azimuthal points
            if itteration < maxItter:
                pass
            else:
                state[currentInd + site] = 1
            currentInd = currentInd + numPoints    
            
            #radial points
            if itteration !=0:
                currentInd = currentInd + numRadials
                
        return state
    
    def build_local_state(self, site):
        '''
        build a single site state at any location on the lattice.
        
        site is the absolute index coordinate of the lattice site
        (use get_SDindex to obtain this in a halfway sensible fashion)
        '''
        if site >= len(self.SDx):
            raise ValueError('lattice doesnt have this many sites')
            
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
            if self.itter>3:
                self.draw_SDlinks(ax, 3, linewidth = 0.5, color = 'firebrick', zorder = zorder)
            else:
                self.draw_SDlinks(ax, self.itter,linewidth = 0.5, color = 'firebrick', zorder = zorder)
        
        pylab.title(title, fontsize=8)
        ax.xaxis.set_visible(False)
        ax.yaxis.set_visible(False)
        ax.set_aspect('equal')
        return

    def plot_map_state(self, map_vect, ax, title = 'ineraction weight', colorbar = False, plot_links = False, cmap = 'winter', autoscale = False, zorder = 1):
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
            vmin = vals[0]
        else:
            vmax = second_biggest
            vmin = -second_biggest
    
        pylab.sca(ax)
        pylab.scatter(self.SDx, self.SDy,c =  mColors, s = mSizes, marker = 'o', edgecolors = 'k', cmap = cm, vmax = vmax, vmin = vmin, zorder = zorder)
        if colorbar:
            cbar = pylab.colorbar(fraction=0.046, pad=0.04)
            cbar.set_label('interaction energy (AU)', rotation=270)
              
        if plot_links:
            if self.itter>3:
                self.draw_SDlinks(ax, 3, linewidth = 0.5, color = 'firebrick', zorder = zorder)
            else:
                self.draw_SDlinks(ax, self.itter,linewidth = 0.5, color = 'firebrick', zorder = zorder)
        
        pylab.title(title, fontsize=8)
        ax.xaxis.set_visible(False)
        ax.yaxis.set_visible(False)
        ax.set_aspect('equal')
        return

    def get_end_state_plot_points(self,scaleFactor = 0.5, maxItter = -1):
        '''
        find end coordinate locations part way along each resonator so that
        they can be used to plot the field at both ends of the resonator.
        (Will retun all values up to specified itteration. Default is the whole thing)
        
        Scale factor says how far appart the two points will be: +- sclaeFactor.2 of the total length
        
        returns the polt points as collumn matrix
        '''
        if scaleFactor> 1:
            raise ValueError('scale factor too big')
            
        if maxItter > self.itter:
            raise ValueError('dont have this many itterations')
        elif maxItter <0:
            maxItter = self.itter
            
        if maxItter < self.itter:
            size = self.get_SDindex(0,maxItter+1,az = True)
        else:
            size = len(self.SDx)
        plot_points = numpy.zeros((size*2, 2))
        
        resonators = self.get_all_resonators(maxItter)
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
            raise ValueError('You screwed around with the mode type. It must be FW or HW.')
        
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
            if self.itter>3:
                self.draw_SDlinks(ax, 3, linewidth = 0.5, color = 'firebrick', zorder = zorder)
            else:
                self.draw_SDlinks(ax, self.itter,linewidth = 0.5, color = 'firebrick', zorder = zorder)
        
        pylab.title(title, fontsize=8)
        ax.xaxis.set_visible(False)
        ax.yaxis.set_visible(False)
        ax.set_aspect('equal')
        
#        return plotPoints
        return mColors









    
if __name__=="__main__":    
    
    '''
            
    Sample syntax:
        #####
        #loading precalculated layout
        #####
        from LayoutGenerator4 import PlanarLayout
        test = PlanarLayout(file_path = '7gon_3vertex_ 2.pkl')
        
        #####
        #making new layout
        #####
        from LayoutGenerator4 import PlanarLayout
        test = PlanarLayout(gon = 7, vertex = 3,side = 1, radius_method = 'lin')
        test.populate(maxItter = 1)
        
        #####
        #saving computer layout
        #####
        test.save( name = 'filename.pkl') #filename can be a full path, but must have .pkl extension
    '''
    
    def ring_coords(radius, numSteps = 100):
        thetas= numpy.linspace(0,2*numpy.pi, numSteps)
        rs = radius * numpy.ones(numSteps)
        
        xs = rs*numpy.cos(thetas)
        ys = rs*numpy.sin(thetas)
        return [xs, ys]
            
    
    
    
    #########
    ###Make Everything
    ########
    
    ###options
    recalculate = True
#    recalculate = False
    
    makeLatticePlots = True   #plot different itterations of lattice generation
#    makeLatticePlots = False


    doThirdOrder = False
    doFourthOrder = False    #include the fourth order itteration
    
    plotLinks = True      #show the effective model links on the lattice itterations
#    plotLinks = False
    

#    makeJCHplot = False
    makeJCHplot = True

#    testHamiltonian = False
    testHamiltonian = True

#    testEigenvectors = True
    testEigenvectors = False
    
    
    if recalculate:
#        test = PlanarLayout(gon = 7, vertex = 3, side =1, radius_method = 'lin')
#        test = PlanarLayout(gon = 7, vertex = 3, side =1, radius_method = 'exp')
        
#        test = PlanarLayout(gon = 5, vertex = 4, side =1, radius_method = 'lin')
#        test = PlanarLayout(gon = 5, vertex = 4, side =1, radius_method = 'exp')

#        test = PlanarLayout(gon = 4, vertex = 4, side =1, radius_method = 'lin')  #square
#        test = PlanarLayout(gon = 4, vertex = 4, side =1, radius_method = 'exp')

        # test = PlanarLayout(gon = 4, vertex = 5, side =1, radius_method = 'lin')
        
#        test = PlanarLayout(gon = 3, vertex = 6, side =1, radius_method = 'lin')    #triangular
        
        # test = PlanarLayout(gon = 3, vertex = 7, side =1, radius_method = 'lin')


#        test = PlanarLayout(gon = 5, vertex = 3, side =1, radius_method = 'lin')
#        test = PlanarLayout(gon = 6, vertex = 3, side =1, radius_method = 'lin')    #hexagonal
        test = PlanarLayout(gon = 7, vertex = 3, side =1, radius_method = 'lin')
#        test = PlanarLayout(gon = 8, vertex = 3, side =1, radius_method = 'lin')
#        test = PlanarLayout(gon = 9, vertex = 3, side =1, radius_method = 'lin')
#        test = PlanarLayout(gon = 10, vertex = 3, side =1, radius_method = 'lin')



        
#        test = PlanarLayout(gon = 4, vertex = 4, side =1, radius_method = 'lin')
#        test = PlanarLayout(gon = 5, vertex = 4, side =1, radius_method = 'lin')
#        test = PlanarLayout(gon = 6, vertex = 4, side =1, radius_method = 'lin')
#        test = PlanarLayout(gon = 7, vertex = 4, side =1, radius_method = 'lin')


#        test = PlanarLayout(gon = 3, vertex = 5, side =1, radius_method = 'lin') #hinky
#        test = PlanarLayout(gon = 4, vertex = 5, side =1, radius_method = 'lin')
#        test = PlanarLayout(gon = 5, vertex = 5, side =1, radius_method = 'lin')
#        test = PlanarLayout(gon = 6, vertex = 5, side =1, radius_method = 'lin')
#        test = PlanarLayout(gon = 7, vertex = 5, side =1, radius_method = 'lin')

#        test = PlanarLayout(gon = 3, vertex = 7, side =1, radius_method = 'lin')

#        test = PlanarLayout(gon = 3, vertex = 12, side =1, radius_method = 'lin')
#        test = PlanarLayout(gon = 12, vertex = 3, side =1, radius_method = 'lin')

#        test = PlanarLayout(gon = 5, vertex = 3, side =1, radius_method = 'lin')
#        test = PlanarLayout(gon = 4, vertex = 3, side =1, radius_method = 'lin')
        

#        test = PlanarLayout(gon = 7, vertex = 3, side =1, radius_method = 'lin', modeType = 'HW')
#        test = PlanarLayout(gon = 7, vertex = 3, side =1, radius_method = 'lin', modeType = 'FW')
#        test = PlanarLayout(gon = 8, vertex = 3, side =1, radius_method = 'lin', modeType = 'HW')
        
        
        test.itter_generate()
        xs = test.get_xs()
        ys = test.get_ys()
        
        test.itter_generate()
        xs2 = test.get_xs()
        ys2 = test.get_ys()
        
        if doThirdOrder:
            test.itter_generate()
            xs3 = test.get_xs()
            ys3 = test.get_ys()
        
        if doFourthOrder:
            test.itter_generate()
            xs4 = test.get_xs()
            ys4 = test.get_ys()
        
        
        test.generate_semiduals()
        
        test.generate_Hamiltonian()

        
        
        [firstRingX ,firstRingY] = ring_coords(test.radii[0])
        [secondRingX ,secondRingY] = ring_coords(test.radii[1])
        [thirdRingX ,thirdRingY] = ring_coords(test.radii[2])
        if doThirdOrder:
            [fourthRingX ,fourthRingY] = ring_coords(test.radii[3])
        #[fifthRingX ,fifthRingY] = ring_coords(test.radii[4])
        
        #[firstRingX ,firstRingY] = ring_coords(test.radii[0])
        #[secondRingX ,secondRingY] = ring_coords(test.radii[0]*2)
        #[thirdRingX ,thirdRingY] = ring_coords(test.radii[0]*3)
        #[fourthRingX ,fourthRingY] = ring_coords(test.radii[0]*4)
        #[fifthRingX ,fifthRingY] = ring_coords(test.radii[0]*5)
    
    
    
    
    
    ####################################################
    ##Lattice Plots
    ########
    if makeLatticePlots:
        pylab.figure(1)
        pylab.clf()
        ax1 = pylab.subplot(1,1,1)
        pylab.plot(test.x0, test.y0, 'kd')
        
        
        #first itteration
        pylab.plot(xs, ys, 'b.')
        
        [rxs_in, rys_in, rxs_out, rys_out] = test.get_radials(1)
        pylab.plot(rxs_in, rys_in, 'r+')
        pylab.plot(rxs_out, rys_out, 'r+')
        
        test.draw_radials(ax1, 1, color = 'g')
        
        
        pylab.plot(firstRingX, firstRingY, 'c-')
        pylab.plot(secondRingX, secondRingY, 'b-')
        pylab.plot(thirdRingX, thirdRingY, 'k-')
        
        sdx, sdy = test.get_all_semidual_points(maxItter = 1)
        pylab.plot(sdx, sdy,color =  'gold', marker = '+', linestyle = '')
        
        if plotLinks:
            test.draw_SDlinks(ax1, 1)
        
        ax1.set_aspect('equal')
        pylab.title('first itteration')
        pylab.show()
        
        #####################################################
        
        
        pylab.figure(2)
        pylab.clf()
        ax1 = pylab.subplot(1,1,1)
        pylab.plot(test.x0, test.y0, 'kd')
        
        #itterate out to the first ring
        pylab.plot(xs2, ys2, 'b.')
        
        [rxs_in, rys_in, rxs_out, rys_out] = test.get_radials(2)
        pylab.plot(rxs_in, rys_in, 'r+')
        pylab.plot(rxs_out, rys_out, 'r+')
        
        test.draw_radials(ax1, 1, color = 'g')
        test.draw_radials(ax1, 2, color = 'g')
        
        pylab.plot(firstRingX, firstRingY, 'c-')
        pylab.plot(secondRingX, secondRingY, 'b-')
        pylab.plot(thirdRingX, thirdRingY, 'k-')
        
        
        sdx, sdy = test.get_all_semidual_points(maxItter = 2)
        pylab.plot(sdx, sdy,color =  'gold', marker = '+', linestyle = '')
        
        if plotLinks:
            test.draw_SDlinks(ax1, 2)
        
        
        ax1.set_aspect('equal')
        pylab.title('second itteration')
        pylab.show()
        
        ######################################################
        
        if doThirdOrder:
            pylab.figure(3)
            pylab.clf()
            ax1 = pylab.subplot(1,1,1)
            pylab.plot(test.x0, test.y0, 'kd')
            
            #itterate out to the first ring
            pylab.plot(xs3, ys3, 'b.')
            
            [rxs_in, rys_in, rxs_out, rys_out] = test.get_radials(3)
            pylab.plot(rxs_in, rys_in, 'r+')
            pylab.plot(rxs_out, rys_out, 'r+')
            
            test.draw_radials(ax1, 1, color = 'g')
            test.draw_radials(ax1, 2, color = 'g')
            test.draw_radials(ax1, 3, color = 'g')
            
            
            pylab.plot(firstRingX, firstRingY, 'c-')
            pylab.plot(secondRingX, secondRingY, 'b-')
            pylab.plot(thirdRingX, thirdRingY, 'k-')
            pylab.plot(fourthRingX, fourthRingY, 'y-')
            
            sdx, sdy = test.get_all_semidual_points(maxItter = 3)
            pylab.plot(sdx, sdy,color =  'gold', marker = '+', linestyle = '')
            
            #test.draw_SDlinks(ax1, 3, color = 'blueviolet')
            #test.draw_SDlinks(ax1, 3, color = 'firebrick')
            #test.draw_SDlinks(ax1, 3, color = 'gold')
            if plotLinks:
                test.draw_SDlinks(ax1, 3)
            
            ax1.set_aspect('equal')
            pylab.title('third itteration')
            pylab.show()
        
        #####################################################
        
        if doFourthOrder:
            pylab.figure(4)
            pylab.clf()
            ax1 = pylab.subplot(1,1,1)
            pylab.plot(test.x0, test.y0, 'kd')
            
            #itterate out to the first ring
            pylab.plot(xs4, ys4, 'b.')
            
            [rxs_in, rys_in, rxs_out, rys_out] = test.get_radials(3)
            pylab.plot(rxs_in, rys_in, 'r+')
            pylab.plot(rxs_out, rys_out, 'r+')
            
            test.draw_radials(ax1, 1, color = 'g')
            test.draw_radials(ax1, 2, color = 'g')
            test.draw_radials(ax1, 3, color = 'g')
            test.draw_radials(ax1, 4, color = 'g')
            
            
            sdx, sdy = test.get_all_semidual_points(maxItter = 4)
            pylab.plot(sdx, sdy,color =  'gold', marker = '+', linestyle = '')
            
            
            pylab.plot(firstRingX, firstRingY, 'c-')
            pylab.plot(secondRingX, secondRingY, 'b-')
            pylab.plot(thirdRingX, thirdRingY, 'k-')
            pylab.plot(fourthRingX, fourthRingY, 'y-')
            pylab.plot(fifthRingX, fifthRingY, 'deepskyblue')
            
            if plotLinks:
                test.draw_SDlinks(ax1, 4)
            
            ax1.set_aspect('equal')
            pylab.title('fourth itteration')
            pylab.show()
    else:
        pylab.figure(1)
        pylab.clf()
        pylab.figure(2)
        pylab.clf()
        pylab.figure(3)
        pylab.clf()
        pylab.figure(4)
        pylab.clf()
    
    
    
    
    
    ###############################
    ##plot the resonator and model lattices
    ######
    if makeJCHplot:
        pylab.figure(5)
        pylab.clf()
        ax1 = pylab.subplot(1,2,1)
        
        
        #plot the effective model
        if test.itter > 3:
            plotItt = 3
        else:
            plotItt = test.itter
        sdx, sdy = test.get_all_semidual_points(maxItter = plotItt)
        #pylab.plot(sdx, sdy,color =  'gold', marker = '+', linestyle = '')
        #pylab.plot(sdx, sdy,color =  'b', marker = '+', linestyle = '')
        pylab.plot(sdx, sdy,color =  'deepskyblue', marker = '+', linestyle = '', markersize = 10)
        test.draw_SDlinks(ax1, plotItt, linewidth = 1, color = 'firebrick')
        
        #plot the underlying resonator
        alphaFactor = 0.4
        pylab.plot(firstRingX, firstRingY, 'c-', alpha = alphaFactor)
        pylab.plot(secondRingX, secondRingY, 'b-', alpha = alphaFactor)
        pylab.plot(thirdRingX, thirdRingY, 'k-', alpha = alphaFactor)
    #    pylab.plot(fourthRingX, fourthRingY, 'y-', alpha = alphaFactor)
        
        pylab.plot(test.x0, test.y0, 'kd',alpha = alphaFactor)
        if doThirdOrder:
            pylab.plot(xs3, ys3, 'b.',alpha = alphaFactor)
        else:
            pylab.plot(xs2, ys2, 'b.',alpha = alphaFactor)
            
        [rxs_in, rys_in, rxs_out, rys_out] = test.get_radials(plotItt)
        pylab.plot(rxs_in, rys_in, 'r+',alpha = alphaFactor)
        pylab.plot(rxs_out, rys_out, 'r+',alpha = alphaFactor)
        
        for i in range(1, plotItt+1):
            test.draw_radials(ax1, i, color = 'g',alpha = alphaFactor)


        
        ax1.set_aspect('equal')
        pylab.title('circuit QED effective model')
        
        
        
        #####
        ax2 = pylab.subplot(1,2,2)
        pylab.plot(test.x0, test.y0, 'kd')
        
        #itterate out to the first ring
        if doThirdOrder:
            pylab.plot(xs3, ys3, 'b.')
        else:
            pylab.plot(xs2, ys2, 'b.')
        
        [rxs_in, rys_in, rxs_out, rys_out] = test.get_radials(plotItt)
        pylab.plot(rxs_in, rys_in, 'r+')
        pylab.plot(rxs_out, rys_out, 'r+')
        
        for i in range(1, plotItt+1):
            test.draw_radials(ax2, i, color = 'g')
        
        
        pylab.plot(firstRingX, firstRingY, 'c-')
        pylab.plot(secondRingX, secondRingY, 'b-')
        pylab.plot(thirdRingX, thirdRingY, 'k-')
    #    pylab.plot(fourthRingX, fourthRingY, 'y-')
        
        #sdx, sdy = test.get_all_semidual_points(maxItter = 3)
        #pylab.plot(sdx, sdy,color =  'gold', marker = '+', linestyle = '')
        ax2.set_aspect('equal')
        pylab.title('resonator layout')
        
        
        pylab.show()
    else:
        pylab.figure(5)
        pylab.clf()
        
    ##########################################
    
    
    
    
    
    
    
    
    
    
    
    ######################################################
    ##The hamiltonian now
    ##########
    
    
    if testHamiltonian:
        subH = test.get_sub_Hamiltonian(itteration = 1)
        
        pylab.figure(6)
        pylab.clf()
        ax1 = pylab.subplot(1,2,1)
#        pylab.imshow(subH, cmap = 'jet')
#        pylab.imshow(subH, cmap = 'seismic')
#        pylab.imshow(subH, cmap = 'ocean')
#        pylab.imshow(subH, cmap = 'hot', vmin = -1, vmax = 1.5)
        pylab.imshow(subH, cmap = 'winter')
        pylab.title('tight-binding Hamiltonian of the \n zeroth and first itterations')
        
#        ax2 = pylab.subplot(1,2,2)
#        pylab.imshow(test.H, cmap = 'jet')
#        pylab.title(' full tight-binding Hamiltonian')
        
        ax2 = pylab.subplot(1,2,2)
#        pylab.imshow(test.H_HW, cmap = 'jet', vmin = 0, vmax = 2*test.t)
        pylab.imshow(test.H_HW, cmap = 'winter', vmin = 0, vmax = 2*test.t)
        pylab.title(' full tight-binding Hamiltonian')
        
        pylab.show()
        
        
        
        ######################################################
        maxItt = 1
        [Es, Psis, Eorder] = test.get_eigs(maxItter = maxItt)
        
        pylab.figure(7)
        pylab.clf()
        ax1 = pylab.subplot(1,2,1)
        xs = scipy.arange(0,len(Es),1)
        pylab.plot(xs, Es[Eorder], 'b.')
        pylab.title('tight-binding Spectrum of the \n zeroth and first itterations')
        
        ax2 = pylab.subplot(1,2,2)
        xs = scipy.arange(0,len(test.Es),1)
        pylab.plot(xs, test.Es[test.Eorder], 'b.')
        pylab.title(' full tight-binding Spectrum')
        pylab.show()
        
        
            
        
        
        #######################################################
        eigNum = 25
        eigAmps = test.Psis[:,test.Eorder[eigNum]]
        eigProbs = numpy.abs(eigAmps)**2
        mSizes = eigProbs * len(eigProbs)*100
        
        
        pylab.figure(8)
        pylab.clf()
        ax1 = pylab.subplot(1,2,1)
        xs = scipy.arange(0,len(test.Es),1)
        pylab.plot(xs, test.Es[test.Eorder], 'b.')
        
        pylab.plot(xs[eigNum],test.Es[test.Eorder[eigNum]], color = 'firebrick' , marker = '.', markersize = '10' )
        pylab.title(' full tight-binding Spectrum')
        
            
        
        
        ax2 = pylab.subplot(1,2,2)
#==============================================================================
# #        if test.itter > 3:
# #            sdx, sdy = test.get_all_semidual_points(maxItter = 3)
# #        else:
# #            sdx, sdy = test.get_all_semidual_points()
# #        pylab.scatter(sdx, sdy,color =  'goldenrod', s = mSizes, marker = 'o', edgecolors = 'k')
# #        if test.itter > 3:
# #            test.draw_SDlinks(ax2, 3, linewidth = 0.5, color = 'firebrick')
# #        else:
# #            test.draw_SDlinks(ax2, test.itter, linewidth = 0.5, color = 'firebrick')
# #        
# #        titleStr = 'eigenvector weight : ' + str(eigNum)
# #        pylab.title(titleStr)
# #        ax2.set_aspect('equal')
#==============================================================================

        titleStr = 'eigenvector weight : ' + str(eigNum)
        test.plot_layout_state(eigAmps, ax2, title = titleStr, colorbar = True, plot_links = True, cmap = 'Wistia')

#        subItt = 1
#        [Es, Psis, Eorder] = test.get_eigs(maxItter = subItt)
#        eigNum2 = 0
#        eigAmps2 = Psis[:,Eorder[eigNum2]] + 1.0*10**-1
#        titleStr = 'eigenvector weight : ' + str(eigNum)
#        test.plot_layout_state(eigAmps2, ax2, title = titleStr, colorbar = True, plot_links = True, cmap = 'Wistia')
        
        pylab.show()
        
        #test half wave versus full wave
        maxItt = 1
        [Es, Psis, Eorder] = test.get_eigs(maxItter = maxItt)
        
        EsHW, PsisHW = scipy.linalg.eigh(test.H_HW)
        EorderHW = numpy.argsort(EsHW)
        
        typeHW = EsHW[EorderHW[0:len(test.Es)]] + test.internalBond
        typeFW = EsHW[EorderHW[len(test.Es):]] - test.internalBond
        
        
        pylab.figure(11)
        pylab.clf()
        ax1 = pylab.subplot(1,2,1)
        xs = scipy.arange(0,len(Es),1)
        pylab.plot(xs, Es[Eorder], 'r.')
        xs = scipy.arange(0,len(test.Es),1)
        pylab.plot(xs, test.Es[test.Eorder], 'b.')
        pylab.title('tight-binding Spectrum of full wave')
        
        ax2 = pylab.subplot(1,2,2)
        pylab.plot(xs, test.Es[test.Eorder], 'r.')
        pylab.plot(xs, numpy.real(typeHW), 'b.', label = 'half wave')
        pylab.plot(xs, numpy.real(typeFW), marker = '.', color = 'deepskyblue', linestyle = '', label = 'full wave')
        pylab.title(' full tight-binding Spectrum FW v HW')
        ax2.legend(loc = 'upper left')
        pylab.show()
        
        pylab.figure(12)
        pylab.clf()
        ax1 = pylab.subplot(1,2,1)
        maxInd = len(test.Es)-1
        test.plot_end_layout_state(test.Psis[:,maxInd], ax1, title = 'end state local weights', colorbar = True, plot_links = True, cmap = 'Wistia')
        
        ax2 = pylab.subplot(1,2,2)
        test.plot_layout_state(test.Psis[:,maxInd], ax2, title = 'single end weights', colorbar = True, plot_links = True, cmap = 'Wistia')
        
        pylab.show()
        
        
    else:
        pylab.figure(6)
        pylab.clf()
        pylab.figure(7)
        pylab.clf()
        pylab.figure(8)
        pylab.clf()
        
        pylab.figure(11)
        pylab.clf()
    
    
    
    
    
    if testEigenvectors:
        ####temporary test for symmetry of the tight binding Hamiltonian
        ## And of the orthonormality of the computer eigenvectors
        #Es = test.Es
        #H = test.H
        #Psis = test.Psis
        
        H = test.H
        import scipy.linalg
        Es, Psis = scipy.linalg.eigh(H)
        
        #H = test.H + 50*numpy.identity(test.H.shape[0])
        #import scipy.linalg
        ##Es, Psis = scipy.linalg.eigh(H)
        #Es, Psis = numpy.linalg.eig(H)
        
        
        pylab.figure(2)
        pylab.clf()
        pylab.plot(numpy.real(Es), 'b')
        pylab.plot(numpy.imag(Es), 'r')
        pylab.title('real and imaginary parts of eigenvalues')
        pylab.show()
        
        pylab.figure(3)
        pylab.clf()
        ax1 = pylab.subplot(1,2,1)
#        pylab.imshow(H, cmap = 'jet')
        pylab.imshow(H, cmap = 'winter')
        pylab.title('H')
        
        ax2 = pylab.subplot(1,2,2)
#        pylab.imshow(H - numpy.transpose(H), cmap = 'jet')
        pylab.imshow(H - numpy.transpose(H), cmap = 'hot')
        pylab.title('H - transpose(H)')
        
        pylab.show() 
        
        
        #    b = numpy.dot(test.Psis, numpy.transpose(test.Psis))
        #    b = numpy.dot(test.Psis, numpy.transpose(numpy.conj(test.Psis)))
        b = numpy.dot(numpy.transpose(numpy.conj(Psis)), Psis)
        c = numpy.dot(Psis, numpy.linalg.inv(Psis))
        
        for i in range(0, len(Es)):
            vect = Psis[:,i]
            print('eignenvector : ' + str(i))
            print(numpy.linalg.norm(vect))
            print('\n')
            
        #        print numpy.dot(vect, vect)
            print(numpy.dot(numpy.conj(vect), vect))
            print(b[i,i])
            print(c[i,i])
            print('\n\n')
        
        
        pylab.figure(4)
        pylab.clf()
        ax1 = pylab.subplot(1,3,1)
        pylab.imshow(numpy.abs(b), cmap = 'seismic')
        pylab.title('transpose conjugate of Psis *Psi')
        #    pylab.colorbar()
        pylab.colorbar(fraction=0.046, pad=0.04)
        
        ax2 = pylab.subplot(1,3,2)
        pylab.imshow(numpy.abs(c), cmap = 'seismic')
        pylab.title('Psi * inverse of Psis')
        
        ax3 = pylab.subplot(1,3,3)
        pylab.imshow( numpy.abs(numpy.linalg.inv(Psis) - numpy.transpose(numpy.conj(Psis))), cmap = 'seismic')
        pylab.colorbar(fraction=0.046, pad=0.04)
        pylab.title('difference between the two')
        
        pylab.show() 
        
        temp = numpy.where(numpy.abs(b) == numpy.max(numpy.abs(b)))
        ind = temp[0][0]
        vect = test.Psis[:,ind]
        print(numpy.max(numpy.abs(b)))
        print(numpy.linalg.norm(vect))
        print(numpy.dot(vect, vect))
        print(numpy.dot(numpy.conj(vect), vect))
        
        
        
        
        
#    for ind in range(0, len(test.SDx)):
#        print numpy.sum(test.H[:,ind])
