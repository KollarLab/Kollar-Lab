#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 15 12:20:48 2017

@author: kollar2

code to draw hyperbolic tilings in the Poincare disc model

First section calculates a bunch of details about the tiling 
e.g. where the corners of the first polygon are and where
the centers of the circles which form the edges of the first polygon are,
and lattice spacings.

The second section is all the methods necessary for drawing actual pictures of the tilings
and some sample code for making the pictures.
Note: it cannot calculate a kagome-like tiling without also calculating the layout tiling 
underneath at the same time, so the functions automatically do both. If you want them to
appear in different plots, just feed it two different subplots, and it will seperate them
but it has to claculate them together.



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


######
#pick the gon, and vertex
#####

#nn = 7
#vertexNum = 3
#itterations = 3

#nn = 7
#vertexNum = 3
#itterations = 1

nn = 3
vertexNum = 7
itterations = 1

#nn = 5
#vertexNum = 4
#itterations = 3





#####
#set up geometrical factors
######

trigFactor = numpy.tan(numpy.pi/vertexNum)
preFactor = trigFactor**2 + 1
centerR = numpy.sqrt(1/(1-preFactor*(numpy.sin(numpy.pi/nn))**2))
rr = numpy.sqrt(numpy.abs(preFactor))*numpy.sin(numpy.pi/nn)*centerR



#####
#set up arraus and plot points
######
thetas = numpy.linspace(0,2*numpy.pi, 1000)

borderx = numpy.cos(thetas)
bordery = numpy.sin(thetas)


########
#most of these points aren't used much anymore
center1x = centerR*numpy.cos(numpy.pi/nn)
center1y = centerR*numpy.sin(numpy.pi/nn)

center2x = centerR*numpy.cos(numpy.pi/nn)
center2y = -centerR*numpy.sin(numpy.pi/nn)

intersectionx = centerR*numpy.cos(numpy.pi/nn) - numpy.sqrt(rr**2 - (centerR*numpy.sin(numpy.pi/nn))**2)
intersectiony = 0

intersectionx2 = centerR*numpy.cos(numpy.pi/nn) - numpy.sqrt(3) * centerR*numpy.sin(numpy.pi/nn) 
intersectiony2 = 0

center3x = centerR
center3y = 0
#center4x = numpy.cos(2*2*numpy.pi/nn)*centerR
#center4y = numpy.sin(2*2*numpy.pi/nn)*centerR
center4x = numpy.cos(2*numpy.pi/nn)*centerR
center4y = numpy.sin(2*numpy.pi/nn)*centerR
circle3xs = rr*numpy.cos(thetas) + center3x
circle3ys = rr*numpy.sin(thetas) + center3y
circle4xs = rr*numpy.cos(thetas) + center4x
circle4ys = rr*numpy.sin(thetas) + center4y

borderIntersectionx = (centerR**2 - rr**2 + 1)/(2*centerR)

kagomePointx = centerR -rr
kagomePointy = 0
############


#######
#first corner of the central polygon
#########
vertex1x = centerR*numpy.cos(numpy.pi/nn) - trigFactor * centerR*numpy.sin(numpy.pi/nn) 
vertex1y = 0
vertex1 = vertex1x + 1j*vertex1y

vertex2x = numpy.cos(2*numpy.pi/nn)*vertex1x
vertex2y = numpy.sin(2*numpy.pi/nn)*vertex1x
vertex2 = vertex2x + 1j*vertex2y

kagomeVertex1x = centerR - numpy.sqrt(trigFactor**2 +1)*numpy.sin(numpy.pi/nn)*centerR
kagomeVertex1y = 0
kagomeVertex1 = kagomeVertex1x + 1j*kagomeVertex1y

kagomeVertex2x = numpy.cos(2*numpy.pi/nn)*kagomeVertex1x
kagomeVertex2y = numpy.sin(2*numpy.pi/nn)*kagomeVertex1x
kagomeVertex2 = kagomeVertex2x + 1j*kagomeVertex2y



#####
#print curvature results
#####
argg = (vertex1-vertex2)/(1-vertex1*numpy.conj(vertex2))
kagomeArgg = (kagomeVertex1-kagomeVertex2)/(1-kagomeVertex1*numpy.conj(kagomeVertex2))

latticeConstant = 2*numpy.arctanh(numpy.abs(argg))
kagomeLatticeConstant = 2*numpy.arctanh(numpy.abs(kagomeArgg))

print str(nn) + '-gon:, ' + str(vertexNum) + ' vertex:'
print 'lattice constant = ' + str(latticeConstant)

print ''

print str(nn) + '-gon:, ' + str(vertexNum) + ' vertex, medial: '
print 'meadial lattice constant = ' + str(kagomeLatticeConstant)


#####
#check against Coxeter Formulas
#####

phi = numpy.arccosh(numpy.cos(numpy.pi/nn)/numpy.sin(numpy.pi/vertexNum))
chi = numpy.arccosh((1/numpy.tan(numpy.pi/nn))*(1/numpy.tan(numpy.pi/vertexNum)))
psi = numpy.arccosh(numpy.cos(numpy.pi/vertexNum)/numpy.sin(numpy.pi/nn))

print ''
print phi*2
print chi
print psi








#####
#draw a tiling
#####

def hyperline(a, b, numEdgePoints = 20):
    '''
    Draw a geodesic line in the Poincare disc model between two points a and b
    a and b are complex coordinates for the two end points
    
    Returns the line as a complex valued vector
    '''
    lamd=((1-a*numpy.conj(b))/(numpy.abs(1-a*numpy.conj(b))))*numpy.sqrt((1-a*numpy.conj(a))/(1-b*numpy.conj(b)))
        
    z=numpy.linspace(-1,1,numEdgePoints)
        
    z=(z*(a-lamd*b) - (a+lamd*b))/(z*(1-lamd)-(1+lamd))
    return z

def draw_gon(axis, linecorners, goncolor = 'r', linewidth = 1):
    '''
    takes as arguments the subplot in which to draw, and the complex 
    coordinates of the corners of the gon, plus some draw options
    
    Will itterate through the sides of the polygon and physically draw each one
    '''
    for ind in range(len(linecorners)):
        Line = hyperline(linecorners[ind], linecorners[numpy.mod(ind+1, nn)])
        axis.plot(numpy.real(Line), numpy.imag(Line), goncolor, linewidth = linewidth)
    return
    
def store_gon(linecorners):
    '''
    Alternative method to draw_gon so that you can store all the points in the
    edges and save them to draw later
    
    takes as input the complex coordinates of all the corners of the polygon
    
    returns all the edges as vectors, concatenated together
    '''
    for ind in range(len(linecorners)):
        Line = hyperline(linecorners[ind], linecorners[numpy.mod(ind+1, nn)])
        if ind ==0:
            gonlines = Line
        else:
            gonlines = numpy.concatenate((gonlines, Line))
    return gonlines

def reflect_point(z_target, line_point1, line_point2):
    '''
    Computes a proper conformal, metric consistent reflection about a geodesic of the 
    Poincare disc model.
    
    Takes as input the original point (in complex format), and the complex
    coordinates for two points on the lircle of reflection. (Typically
    these will be two corners of a polygon that has already been drawn/claculated)
    
    returns the complex coordinate of the reflected point
    '''
    a = line_point1
    b = line_point2
    part1 = b*numpy.conj(a)-a*numpy.conj(b)+(a-b+b*numpy.conj(b)*a-a*b*numpy.conj(a))*numpy.conj(z_target)
    part2 = (numpy.conj(a)-numpy.conj(b)+b*numpy.conj(b)*numpy.conj(a)-a*numpy.conj(a)*numpy.conj(b)+(a*numpy.conj(b)-numpy.conj(a)*b)*numpy.conj(z_target))
    z_out = part1/part2
    
    return z_out

def reflect_gon_on_edge(corners ,line_point1, line_point2):
    '''
    function to reflect an entire polygon properly.
    
    Takes are argument an array of the corners of the polygon (complex format),
    and two points on the lircle of reflection
    
    returns an array of the corners (vertices) of the reflected polygon
    '''
    zs = numpy.zeros(len(corners))*(1 + 1j)
    for ind in range(0, len(corners)):
        zs[ind] = reflect_point(corners[ind], line_point1, line_point2)
    return zs
    
def itterate_gons(corners, kagomeCorners, axis, kagomeAxis, numItt = 1, kagomeColor = 'b', color = 'r', linewidth = 2):
    '''
    Function to itterate through and create the polygons of the tiling.
    NOTE: THIS VERSION IS OLD!!
    It is the original version that draws as is goes and in therefore extremely slow.
    
    takes as input the corners of the fundamental polygon and the corners of it's medial polygon
    Also takes the two subplots on which to draw and a bunch of draw options.
    
    NumItt is the number of recursive itterations.
    
    Will recursively itterate through and draw as it goes.
    '''
    if numItt  > 0:
#        print numItt
        for ind in range(0, nn):
            newGon = reflect_gon_on_edge(corners, corners[ind], corners[numpy.mod(ind+1, nn)] )
            newKagomeGon = reflect_gon_on_edge(kagomeCorners, corners[ind], corners[numpy.mod(ind+1, nn)] )
            draw_gon(axis, newGon, goncolor = color,linewidth = linewidth)
            draw_gon(kagomeAxis, newKagomeGon, goncolor = kagomeColor, linewidth = linewidth)
            
            #recursively itterate
            itterate_gons(newGon, newKagomeGon, axis, kagomeAxis, numItt-1 , kagomeColor, color, linewidth)
    return

def itterate_and_store_gons(corners, kagomeCorners , lines, kagomeLines, numItt = 1, method = 'NaNseparated'):
    '''
    Function to itterate through and create the polygons of the tiling.
    NOTE: This is the new version
    It's default method is to store all the edges in a single vector by smooshing
    them together with a NaN seperating each closed polygon.
    the alternate method is to not put the NaN separator. This requires a different method of plotting, 
    which is slower, but if you actually want the edges as output, this is probably 
    slightly easier to process
    
    takes as input the corners of the fundamental polygon and the corners of it's medial polygon
    Also takes the two subplots on which to draw and a bunch of draw options.
    
    NEEDS as input the edges of the fundamental polygon and the fundamental polygon of the medial.
    It needs to concatenate to these arrays.
    
    NumItt is the number of recursive itterations.
    
    Will recursively itterate through and return the giant vector of all the gons.
    
    '''
    if numItt  > 0:
#        print numItt
        for ind in range(0, nn):
            newGon = reflect_gon_on_edge(corners, corners[ind], corners[numpy.mod(ind+1, nn)] )
            newKagomeGon = reflect_gon_on_edge(kagomeCorners, corners[ind], corners[numpy.mod(ind+1, nn)] )
            gonLines = store_gon(newGon)
            kagomeGonLines = store_gon(newKagomeGon)
            
            if method == 'NaNseparated':
                lines = numpy.concatenate((lines, gonLines, [numpy.NaN]))
                kagomeLines = numpy.concatenate((kagomeLines, kagomeGonLines, [numpy.NaN]))
            else:
                lines = numpy.concatenate((lines, gonLines))
                kagomeLines = numpy.concatenate((kagomeLines, kagomeGonLines))
                
#            print 'side num = ' +str(ind) + ' , itteration number = ' + str(numItt)
#            print lines.size
            
            #recursively itterate
            [lines, kagomeLines] = itterate_and_store_gons(newGon, newKagomeGon,lines, kagomeLines, numItt-1, method = method)
    return [lines, kagomeLines]


def draw_tiling(corners, kagomeCorners, axis, kagomeAxis, numItt = 1, kagomeColor = 'b', color = 'r', linewidth = 2, method = 'fast'):
    '''
    method to actually draw the tilings.
    
    Takes as input the corners of the fundamental polygons for the regular tiling, and it's medial.
    Then it calls the itterate function and recursively builds both.
    
    Three operating methods:
        draw_while: is the oldest way that draws as it itterates. This is paralytically slow
        nice_lines: draws after, but doesn't use NaN seperation, so it needs to draw each polygon
                    seperately and it's only slightly less paralytic
        fast: draws after and uses NaNseparation, so the plotting is fast
    
    It draws the tiling on the two sublots specified and returns the huge arrays of all the edges.
    '''
    if method == 'draw_while':
        #draw starting gon
        draw_gon(axis, corners, goncolor = color,linewidth = linewidth)
        draw_gon(kagomeAxis, kagomeCorners, goncolor = kagomeColor, linewidth = linewidth)
        
        t1 = time.time()
        itterate_gons(corners, kagomeCorners, axis, kagomeAxis, numItt , kagomeColor, color, linewidth)
        t_el = time.time()-t1
        print 'calculaton time = ' + str(t_el)
        return
    else:
        lines = store_gon(corners)
        kagomeLines = store_gon(kagomeCorners)    
    
        if method == 'nice_lines':
            t1 = time.time()
            [lines, kagomeLines] = itterate_and_store_gons(corners, kagomeCorners,lines, kagomeLines, numItt, method = 'notseperated')
            t_el = time.time()-t1
            print 'calculaton time = ' + str(t_el)
            
            print 'plotting...'
            startdex = 0
            while startdex< len(lines):
                segment = lines[startdex:startdex+numEdgePoints]
                axis.plot(numpy.real(segment), numpy.imag(segment), color, linewidth = linewidth)
                startdex = startdex+numEdgePoints
            startdex = 0
            while startdex< len(kagomeLines):
                segment = kagomeLines[startdex:startdex+numEdgePoints]
                kagomeAxis.plot(numpy.real(segment), numpy.imag(segment), kagomeColor, linewidth = linewidth)
                startdex = startdex+numEdgePoints
            return [lines, kagomeLines]
        
        else:# method == 'fast':
            t1 = time.time()
            [lines, kagomeLines] = itterate_and_store_gons(corners, kagomeCorners,lines, kagomeLines, numItt, method = 'NaNseparated')
            t_el = time.time()-t1
            print 'calculaton time = ' + str(t_el)
            
            print 'plotting...'
            axis.plot(numpy.real(lines), numpy.imag(lines), color, linewidth = linewidth)
            kagomeAxis.plot(numpy.real(kagomeLines), numpy.imag(kagomeLines), kagomeColor, linewidth = linewidth)
            return [lines, kagomeLines]

#corners in the first polygon
angles = scipy.arange(0,2*numpy.pi, 2*numpy.pi/nn)
cornersx = numpy.cos(angles) * vertex1x
cornersy = numpy.sin(angles) * vertex1x
corners = cornersx + 1j*cornersy

kagomeAngles = angles + numpy.pi/nn
kagomeCornersx = numpy.cos(kagomeAngles) * kagomeVertex1x
kagomeCornersy = numpy.sin(kagomeAngles) * kagomeVertex1x
kagomeCorners = kagomeCornersx + 1j*kagomeCornersy


#####
#testing the draw functions
#####

pylab.figure(4)
pylab.clf()
ax1 = pylab.subplot(1,1,1)
ax1.set_aspect('equal')
pylab.plot(borderx, bordery, 'cornflowerblue')
for ind in range(0,nn):    
    #plot the graphene vertices
    vertx = numpy.cos(ind*2*numpy.pi/nn) * vertex1x
    verty = numpy.sin(ind*2*numpy.pi/nn) * vertex1x
    pylab.plot(vertx, verty, color = 'cornflowerblue', marker = 'd')
    
    kagomeVertx = numpy.cos(ind*2*numpy.pi/nn + numpy.pi/nn) * kagomeVertex1x
    kagomeVerty = numpy.sin(ind*2*numpy.pi/nn + numpy.pi/nn) * kagomeVertex1x
    pylab.plot(kagomeVertx,kagomeVerty, color = 'darkgoldenrod', marker = 'd')

#draw_gon(ax1, corners, goncolor = 'r')
#draw_gon(ax1, kagomeCorners, goncolor = 'b')
#
#temp = reflect_gon_on_edge(corners, corners[0], corners[1])
#draw_gon(ax1, temp, 'r')
#
#temp3 = reflect_gon_on_edge(kagomeCorners, corners[0], corners[1])
#draw_gon(ax1, temp3, 'b')

draw_tiling(corners, kagomeCorners, ax1, ax1, numItt = 2)
ax1.axis('off')

pylab.show


#####
#simplest way to draw tiling, two subplots of same figure
#####
#
#fig = pylab.figure(5)
#pylab.clf()
#titleStr = str(nn) + '-gon:, ' + str(vertexNum) + ' vertex:'
#saveStr = str(nn) + 'gon' + str(vertexNum) + 'vertex_' +str(itterations) + 'itterations'
#pylab.suptitle(titleStr)
#ax1 = pylab.subplot(1,2,1)
#ax1.set_aspect('equal')
#ax1.axis('off')
#pylab.plot(borderx, bordery, 'cornflowerblue')
#pylab.title('graphene')
#
#ax2 = pylab.subplot(1,2,2)
#ax2.set_aspect('equal')
#ax2.axis('off')
#pylab.plot(borderx, bordery, 'cornflowerblue')
#pylab.title('kagome')
#
#fig.set_size_inches(11, 5)
#
#
#t = time.time()
##draw_tiling(corners, kagomeCorners, ax1, ax2, numItt = itterations, kagomeColor = 'mediumblue', color = 'goldenrod', linewidth = 1)
##draw_tiling(corners, kagomeCorners, ax1, ax2, numItt = itterations, kagomeColor = 'deepskyblue', color = 'mediumblue', linewidth = 1.5)
##draw_tiling(corners, kagomeCorners, ax1, ax2, numItt = itterations, kagomeColor = 'mediumblue', color = 'darkgoldenrod', linewidth = 1.5)
##draw_tiling(corners, kagomeCorners, ax1, ax2, numItt = itterations, kagomeColor = 'mediumblue', color = 'firebrick', linewidth = 1.5)
#draw_tiling(corners, kagomeCorners, ax1, ax2, numItt = itterations, kagomeColor = 'dodgerblue', color = 'mediumblue', linewidth = 1.5)
#pylab.show()
#t_elapsed = time.time()-t
#
#print ' '
#print ' '
#print 'method = fast , itterations = ' + str(itterations) 
#print 'elapsed time = ' +str(t_elapsed)
#print ' '
#
#savePath = r'/Users/kollar2/Documents/HouckLab/HyperbolicPlanning/PoincarePlots/BB/'
#####savePath = r'/Users/kollar2/Documents/HouckLab/HyperbolicPlanning/PoincarePlots/'
##fig.savefig(savePath + saveStr + '.png' , dpi = 300)
####fig.savefig(saveStr + '.png' , dpi = 300)




#####
#alternate saving method: two seperate figures
#####

fig_layout = pylab.figure(14)
pylab.clf()
ax1 = pylab.subplot(1,1,1)
ax1.set_aspect('equal')
ax1.axis('off')
fig_layout.tight_layout()
pylab.plot(borderx, bordery, 'cornflowerblue')

fig_effective = pylab.figure(15)
pylab.clf()
ax2 = pylab.subplot(1,1,1)
ax2.set_aspect('equal')
ax2.axis('off')
fig_effective.tight_layout()
pylab.plot(borderx, bordery, 'cornflowerblue')

draw_tiling(corners, kagomeCorners, ax1, ax2, numItt = itterations, kagomeColor = 'dodgerblue', color = 'mediumblue', linewidth = 1.5)



fig_layout.set_size_inches(5, 5)
fig_effective.set_size_inches(5, 5)


saveStr = str(nn) + 'gon' + str(vertexNum) + 'vertex_' +str(itterations) + 'itterations'
savePath = r'/Users/kollar2/Documents/HouckLab/HyperbolicPlanning/PoincarePlots/BB/'
#fig_layout.savefig(savePath + saveStr + '_layout.png' , dpi = 300)
#fig_effective.savefig(savePath + saveStr + '_effective.png' , dpi = 300)









######
#optional figures
######

#fig = pylab.figure(9)
#pylab.clf()
#titleStr = str(nn) + '-gon:, ' + str(vertexNum) + ' vertex:'
#saveStr = str(nn) + 'gon' + str(vertexNum) + 'vertex_' +str(itterations) + 'itterations'
#pylab.suptitle(titleStr)
#ax = pylab.subplot(1,1,1)
#pylab.plot(borderx, bordery, 'cornflowerblue')
##draw_tiling(corners, kagomeCorners, ax, ax, numItt = itterations, kagomeColor = 'deepskyblue', color = 'mediumblue', linewidth = 1.5)
#draw_tiling(corners, kagomeCorners, ax, ax, numItt = itterations, kagomeColor = 'dodgerblue', color = 'mediumblue', linewidth = 1.5)
##draw_tiling(corners, kagomeCorners, ax, ax, numItt = itterations, kagomeColor = 'mediumblue', color = 'firebrick', linewidth = 1.5)
##draw_tiling(corners, kagomeCorners, ax, ax, numItt = itterations, kagomeColor = 'mediumblue', color = 'darkgoldenrod', linewidth = 1.5)
#
#ax.set_aspect('equal')
#ax.axis('off')
#
#pylab.show()




#fig = pylab.figure(6)
#pylab.clf()
#titleStr = str(nn) + '-gon:, ' + str(vertexNum) + ' vertex:'
#saveStr = str(nn) + 'gon' + str(vertexNum) + 'vertex_' +str(itterations) + 'itterations'
#pylab.suptitle(titleStr)
#ax1 = pylab.subplot(1,2,1)
#ax1.set_aspect('equal')
#ax1.axis('off')
#pylab.plot(borderx, bordery, 'cornflowerblue')
#pylab.title('both')
#
#ax2 = pylab.subplot(1,2,2)
#ax2.set_aspect('equal')
#ax2.axis('off')
#pylab.plot(borderx, bordery, 'cornflowerblue')
#pylab.title('layout')
#
#fig.set_size_inches(11, 5)
#
#pylab.figure(7)
#ax = pylab.subplot(1,1,1)
#
#draw_tiling(corners, kagomeCorners, ax1, ax1, numItt = itterations, kagomeColor = 'deepskyblue', color = 'mediumblue', linewidth = 1)
#draw_tiling(corners, kagomeCorners, ax2, ax, numItt = itterations, kagomeColor = 'deepskyblue', color = 'mediumblue', linewidth = 1)
#
##fig.savefig(r'/Users/kollar2/Desktop/7_3.png', dpi = 200, transparent=True)
##fig.savefig(r'/Users/kollar2/Desktop/3_7.png', dpi = 200, transparent=True)




###
#doodle
RR = 1/kagomeLatticeConstant
KK = -1/RR**2
print 'Gaussian Curvature = ' + str(KK)
