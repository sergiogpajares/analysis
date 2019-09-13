# -*- coding:utf8 -*-
#
#analysis module
#Author: Sergio García Pajares
#last update: 04-09-2019

#=======================================================================
#==== INFO =============================================================
#=======================================================================
'''
This is an analysis module developed for TEI.

    Contains:
     - Regresion tools that gets the regresion coeficients
     - Plot special funtion that manage errorbars and regresion
     - Some physical constants
     
    
    Dependencies:
     - Numpy
     - Matplotlib
     - Scipy


Author: Sergio García Pajares
Last update: 16-03-2019

'''
__all__ = [
'ponderated_mean','data_plot','linear_regresion','linear_ponderated_regresion',
'linear_origin_regresion','series_ponderated_mean','skip_value','funplot',
'fun3plot',
#'c','e','e0','u0'
]

__version__ = '2.1.0'


'''
Versions history
   2.1.0 
 ---------------------------
   - ponderated_mean: added
   - series_ponderated_mean: added

   2.0.1 
 ---------------------------
   - linear_origin_regresion: added

   
   2.0.0 
 ---------------------------
   - general: inclusion of dynamic lambda use in regresion funtions and
              data_plot
   
   - data_plot: linear coeficientes has been substituided by a regresion
                funtion so, now plot can be used with any function. New
                np parameter added.
   - data_plot: no longer plots grid
   - data_plot: ecolor parameter is now set 'k' by default
   
   - linear_ponderated: debuged 
   
   - data_multi_plots: disapears
   - data_sin_plot: disapears
    
   Previous 
 ---------------------------
    Not recorded

'''

#=======================================================================
#==== IMPORTS ==========================================================
#=======================================================================
from numpy import *
import numpy as np #this dual import is a temporal bug fix before main implementation. It's necesary due to mix of funcs
                   #from main update and before it.
import matplotlib.pyplot as plt

from scipy.stats import itemfreq

#=======================================================================
#==== DATA =============================================================
#=======================================================================

#####################
# --- PHYSICAL CONSTANTS ---
#
#
#c=299792458
#'''
#speed of light in vacuum (m/s)
#'''
#
#e=1.6021766208E-19
#'''
#electron charge magnitude (C)
#'''
#
#e0=8.854187817E-12
#'''
#permittivity of free space (N/m)
#'''
#
#u0=12.566370614E-7
#'''
#permeability of free space (N/A²)
#'''
#
#####################
# --- UNITS MANAGEMENT ---
#
#
#The tranform everything to SI unitis
#'''
#m=1
#
#c = m/1E2
#m = m/1E3
#u = m/1E6
#n = m/1E9  
#'''
#=======================================================================
#==== FUNCTIONS ========================================================
#=======================================================================
#####################
# --- STATISTICS ---
#
def ponderated_mean (x,dx):
    '''
    Calculates ponderated mean
    
        PARAMETERS:
            x,  1D-array-like: x values
            dx, 1D-array-like: x error values
        
        RETURNS:
          _  _
         (x,dx ) tuple
            _
            x,  number: mean
             _
            dx, number: error of the mean
               
    '''
    x=asarray(x,dtype=float)
    dx=asarray(dx,dtype=float)
    
    w=1/(dx**2)
    sw=sum(w)
    swx=sum(w*x)
    
    return((swx/sw),(1/sqrt(sw)))
    
def series_ponderated_mean(x,y,dy):
    '''
    Calculates the ponderated mean of a secuence
    
        PARAMETERS:
            x,  1d-array like: x values
            y,  1d-array like: y values in which ponderated mean is 
                               going to be calculated
            dy, 1d-array like: precision of the y values
        
        RETURNS:
                      _    _
          (x-unique , y , dy )
            
             x-unique, 1d-array: unique values in the original data sheet
             _
             y, 1d-array: mean of the elements with same x value
              _
             dy, 1d-array: error of the mean of the same x elements
        
        
        EXAMPLE 
        
         GIVEN:              RETURNED:
         --------------                  _        _
         | x   y   dy |      X-unique    y       dy
         |------------|  
         | 1  5.1  .2 |         1       5.02    11.01   
         | 1  4.9  .1 |
         |------------|         2       0.04     0.08
         | 2  10   .3 |
         |------------|
         | 1  5.2  .1 |
         |------------|
         | 2  11   .1 |
         | 2  12   .2 |
         | 2   8   .5 |
         |------------|
         | 1  5.0 .05 |
         --------------

    '''
    
    x=asarray(x)
    y=asarray(y)
    dy=asarray(dy)
    
    assert len(x.shape) == 1, "x must be a 1d array like object" 
        #unique flatten arrays if axis is not specified
    assert len(y.shape)  == 1, "y must be a 1d array like object"
    assert len(dy.shape) == 1, "dy must be a 1d array like object"
    assert size(x) == size (y), "x and y mus have the same number of elemnts"
    assert size(y) == size(y), "y and dy must have the same number of elements"
    
    
    x_values  = unique(x) #get unique values of x
    
    y_values  = ones_like(x_values,dtype=float) # create an array in which store y means
    dy_values = ones_like(x_values,dtype=float) # create an array in which store dy of means
    
    i = 0 #initialice #it's the counter of x-values elements
    
    for value in x_values:
        indices = where (x == value) #get indices of x original array
        y_values[i] , dy_values[i] = ponderated_mean(y[indices],dy[indices]) 
            #calculate for chosen values
        i += 1 #next x_value
    
    
    return x_values , y_values , dy_values

#####################
# --- REGRESION ---
#

def linear_regresion(x,y,f=False):
    '''
    Calculates the linear regresion.
    
        PARAMETERS:
            x, 1D-array-like: x points
            y, 1D-array-like: y points
            
            OPTIONAL
            f=False, bool or None: do you want a regresion funtion f(x) 
                                   returned?
                                    True: Returns f
                                    False: Don't return f
                                    
                                    None: Only returns f.
        
        
        RETURNS:
        
        Depends on f optional parameter
            
            if f == False
                (a,b,da,db) tuple
            
            if f == True
                (a,b,da,db,f) tuple
                
            if f == None
                (f)
                    
         For  y = a x + b 
            
           a,  float: a coeficient
           b,  float: b coeficient
           da, float: a error
           db, float: b error
           
           f, funtion: f(x)=ax+b
         
    '''
    #--- preparing inputs ---------------
    x=asarray(x) #turn them into arrays
    y=asarray(y)
    
    x=skip_value(x) #skip None
    y=skip_value(y)
    
    #--- checking -----------------------
    assert len(x.shape)== 1, 'x must be a vector'
    assert len(y.shape)== 1, 'y must be a vector'
    assert size(x)==size(y), 'x and y must have the same number of elements'
    
    #--- calculate twice used values ----
    N=size(x) #number of elements
    
    sx=sum(x) #x sumation
    sx2=sum(x*x) #x square sumation
    
    sy=sum(y) #y sumation
    
    sxy=sum(x*y) # xy sumation
        
    delta=float(N*sx2-sx**2) #common denominator for both paramenteres
    
    
    #--- getting linear coeficients -----
    #ax+b
    a=(N*sxy-(sx*sy))/(delta)
    b=((sx2*sy)-(sx*sxy))/(delta)
    
    #--- getting error ------------------
    sigmay=sqrt((1./(N-2))*sum((y-b-a*x)**2))
    
    da=sigmay*sqrt(N/delta)
    db=sigmay*sqrt(sx2/delta)
    
    if f==False: return(a,b,da,db) #normal return
    
    elif f==None: 
        f=lambda x: a*x+b #Define reggresion func
        return(f)
    
    else: #f==True
        f=lambda x: a*x+b #Define regresion func
        return(a,b,da,db,f)
    
def linear_ponderated_regresion(x,y,dy,f=False):
    '''
    Calculates the linear ponderated regresion.
    
        PARAMETERS:
            x, 1D-array-like: x points
            y, 1D-array-like: y points
            dy, 1D-array-like: y error values
            
            OPTIONAL
            f=False, bool or None: do you want a regresion funtion f(x) 
                                   returned?
                                    True: Returns f
                                    False: Don't return f
                                    
                                    None: Only returns f.
        
        RETURNS:
           
        Depends on f optional parameter
            
            if f == False
                (a,b,da,db) tuple
            
            if f == True
                (a,b,da,db,f) tuple
                
            if f == None
                (f)
                    
         For  y = a x + b 
            
           a,  float: a coeficient
           b,  float: b coeficient
           da, float: a error
           db, float: b error
           
           f, funtion: f(x)=ax+b
         
    '''
    #--- preparing inputs ---------------
    x=asarray(x) #turn them into arrays
    y=asarray(y)
    dy=asarray(dy)
    
    x=skip_value(x) #skip None
    y=skip_value(y)
    dy=skip_value(dy)
    
    #--- checking -----------------------
    assert len(x.shape)== 1, 'x must be a vector'
    assert len(y.shape)== 1, 'y must be a vector'
    assert len(dy.shape)== 1, 'dy must be a vector'
    
    
    assert size(x)==size(y), 'x and y must have the same number of elements'
    assert size(y)==size(dy), 'y and dy must have the same number of elements'
    
    #--- calculate twice used values ----
    N=size(x) #number of elements
    w=1/(dy**2)
    
    sw=sum(w)
    
    swx=sum(w*x) #x sumation
    swx2=sum(w*x*x) #x square sumation
    
    wy=w*y
    swy=sum(wy) #y sumation
    
    swxy=sum(w*x*y) # xy sumation
        
    delta=float(sw*swx2-(swx)**2) #common denominator for both paramenteres
    
    
    #--- getting linear coeficients -----
    #ax+b
    a=(sw*swxy-(swx*swy))/(delta)
    b=(swx2*swy-(swx*swxy))/(delta)
        
    #--- getting error ------------------
    da=sqrt(sw/delta)
    db=sqrt(swx2/delta)
    
    
    if f==False: return(a,b,da,db) #normal return
    
    elif f==None: 
        f=lambda x: a*x+b #Define reggresion func
        return(f)
    
    else: #f==True
        f=lambda x: a*x+b #Define regresion func
        return(a,b,da,db,f)


def linear_origin_regresion(x,y,f=False):
    '''
    Calculates the linear regresion.
    
        PARAMETERS:
            x, 1D-array-like: x points
            y, 1D-array-like: y points
            
            OPTIONAL
            f=False, bool or None: do you want a regresion funtion f(x) 
                                   returned?
                                    True: Returns f
                                    False: Don't return f
                                    
                                    None: Only returns f.
        
        
        RETURNS:
        
        Depends on f optional parameter
            
            if f== False
                (a,b,da,db) tuple
            
            if f== True
                (a,b,da,db,f) tuple
                
            if f== None
                (f)
                    
         For  y = a x
            
           a,  float: a coeficient
           da, float: a error
           
           f, funtion: f(x)=ax+b
         
    '''
    #--- preparing inputs ---------------
    x=asarray(x) #turn them into arrays
    y=asarray(y)
    
    x=skip_value(x) #skip None
    y=skip_value(y)
    
    #--- checking -----------------------
    assert len(x.shape)== 1, 'x must be a vector'
    assert len(y.shape)== 1, 'y must be a vector'
    assert size(x)==size(y), 'x and y must have the same number of elements'
    
    #--- calculate twice used values ----
    N=size(x) #number of elements
    
    sx2=sum(x*x) #x square sumation
    
    sxy=sum(x*y) # xy sumation
    
    
    #--- getting linear coeficients -----
    #ax+b
    a=sxy/sx2
    
    #--- getting error ------------------
    sigmay=sqrt((1./(N-1))*sum((y-a*x)**2))
    
    da=sigmay/sqrt(sx2)
    
    if f==False: return(a,da) #normal return
    
    elif f==None: 
        f=lambda x: a*x #Define regresion func
        return(f)
    
    else: #f==True
        f=lambda x: a*x #Define regresion func
        return(a,da,f)
        
        

_regresiones={
  0                  : linear_regresion,
  'linear'           : linear_regresion,
  1                  : linear_ponderated_regresion,
  'linear_ponderated': linear_ponderated_regresion,
  2                  : linear_origin_regresion,
  'linear_origin'    : linear_origin_regresion
}
'''
This dict is used by data_plot function to get regresion from
the different regresion functions. It's a private var
'''




#####################
# --- PLOTTING ---
#
def funplot (f,xmin,xmax,n=100,fmt='',legend='',title='',xlabel='',
ylabel='',label='',adjust=False):
    '''
    Plots an |R --> |R fuction 
    
        PARAMETERS:
            f, function: |R-->|R function
            xmin, number: lower limit in x axis
            xmax, number: upper limit in x axis
            
            OPTIONAL
            n, int: number of x divisions for x. 
            fmt, str: line format
            legend, str: name of function for legend
            xlabel, str: xlabel
            ylabel, str: tlabel
            adjust, bool: for adjusting axis to xmin and xmax limits.
        
        RETURNS
            A list of Line2D objects representing the plotted data.
            It's generated by matplotlib.pyplot.plot() func.
            
    '''
    #plotting
    x  = np.linspace(xmin,xmax,n) #create x
    y  = f(x) #eval func. Values will be used later.
    
    if fmt != '': #manage fmt provided or not by user
        gr = plt.plot(x,y,fmt,label=label) #plot
    else:
        gr = plt.plot(x,y,    label=label) #plot
    
    #personalization
    if legend != '' : plt.legend() #create legend
    if title != '' : plt.title('u'+title)
    if xlabel != '' : plt.xlabel('u'+xlabel)
    if ylabel != '' : plt.ylabel(ylabel)
    if adjust == True : ptl.axis([xmin,xmax,np.min(y),np.max(y)])
    
    return gr 

def fun3plot (F,x0,x1,y0,y1,z0,z1,n=10, ax=None, normalize=True):
    '''
    Plot a 3D quiver of an F: |R³ ------> |R³ function.
    
        PARAMETERs:
            F, callable: funtion for plot. Must be F(x,y,z) = (v1,v2,v3)
            
            x0, number: x min bound
            x1, number: x min bound
            
            y0, number: x min bound
            y1, number: x max bound
            
            z0, number: x min bound
            z1, number: x max bound
            
            n=10, integer: number of arrows per axis
            
            ax = None: axis in which you want to make the plot.
                       by default it takes current axis
            normalize=True, bool: True normalize arrows and gives it's 
                                  module using color.
        
        RETURNS
            <mpl_toolkits.mplot3d.art3d.Line3DCollection object>
            
    '''
    # get limits
    assert isinstance(n,int) , "n is the number of arrows per axis. It must be an integer." 
    x = np.linspace(x0,x1,n,dtype=float)
    y = np.linspace(y0,y1,n,dtype=float)
    z = np.linspace(z0,z1,n,dtype=float)
    
    # prepare broadcast
    x.shape = ( n, 1, 1)
    y.shape = ( 1, n, 1)
    z.shape = ( 1, 1, n)
    
    # create grid
    x,y,z = np.meshgrid(x,y,z)
    
    # eval funtion
    F_x , F_y , F_z = F(x,y,z)
    
    #in case axis not provided
    if ax == None: ax = plt.gca() 
    
    # plot
    if normalize:
        plot = ax.quiver(x, y, z, F_x*.1, F_y*.1, F_z*.1, length=0.5, normalize = True)#, colors=mpl.cm.jet)
        #add cm here
        #plt.gcf().colorbar(plot, shrink=0.85)
        return plot
    else:
        return ax.quiver(x, y, z, F_x*.1, F_y*.1, F_z*.1)


def data_plot(x,y,fmt='bo',fmtr='b-',dx='',dy='',ecolor='k',label='',xlabel='',ylabel='',ms=3,regresion=None,np=300,extrapole=False,adjust=False):
    '''
    Plots a pair of data. It can also plot its error bars and
    linear regresion. Regresion type can be specified (see below).
    
        PARAMETERS:
            x, 1d array-like: x
            y, 1d array-like: y
            
            OPTIONAL
            fmt='bo', str: format for the points
            label='', str: label for data
            dx='', 1d array-like or number: error values for x
            dy='', 1d array-like or number: error values for y
            ecolor='k', str: format for the color bar
            ms=3, number: marker size
            
            -----
            y=ax+b
            regresion=None, : Do we want to draw linear regresion?
                                Values can be: 
                                
              =====================================================
              |               REGRESION POSIBILITIES              |
              |===================================================|
              |             KEY          |          TYPE          |
              |---+-----------------------------------------------|
              | 0 | linear               | linear                 |
              | 1 | linear_ponderated    | linear ponderated      |
              | 2 | linear_origin        | linear crossing origin |
              |   |                      |                        |
              |--------------------------+------------------------|
              | None  (or True)          | linear or linear pon-  |
              |                          | derated depending on   |
              |                          | dy type (number or     |
              |                          | array like)            |
              |--------------------------+------------------------|
              | False                    | No regresion drawn     |
              -----------------------------------------------------
              | Function f(x) specififed by the user              |
              =====================================================
                                   
            fmtr='b-', str: format for regresion line
            extrapole=False, bool: Do we want regresion crossing
            np=300, int: Number of points used for regresion
            
            -----
            adjust=False, bool: It adjust graphic axis to data
        
        RETURNS: 
            if regresion is drawn:
                tuple (gr,regr)
                    gr:   points draw object
                    regr: regresion draw object
                        
            if regresion isn't drawn:
                gr: points draw object
            
            
    '''
    x = asarray(x,dtype=float) #work with arrays
    y = asarray(y,dtype=float)
    
    x = skip_value(x) #skip none values
    y = skip_value(y)
    
    assert size(x) == size(y), "x and y must have the same number of elements"
    
    
    #---ERROR BARS---
    if type(dx)!= str or type(dy)!= str : #Check if one of them is introduced
        grb=plt.errorbar(x,y,xerr=dx,yerr=dy,fmt=fmt,ms=ms,ecolor=ecolor) #isn't labeled
    #---NORMAL PLOTTING---
    gr=plt.plot(x,y,fmt,label=label,ms=ms)
    '''
    For avoinding weird legend as temporal normal plotting must be
    overplot on errorbar plotting. Only this second plot will be
    labeled so legend will not be modified.
    '''
    
    #---SETTING UP FORMAT---
    if xlabel != '': plt.xlabel(xlabel)
    if ylabel != '': plt.ylabel(ylabel)
    if adjust : plt.axis([min(x),max(x),min(y),max(y)])
    
    
    
    rtype=type(regresion) #it's used serveral times
    #---REGRESION----
    '''
        SERGIO DEL FUTURO, PLANTEATE CAMBIAR ESTA PUTA LOCURA DE CONDICIONALES
        DE SÁBDO A LAS 12PM :( POR UN DICCIONARIO DE SUBRUTINAS QUE 
        ACTÚE EN FUNCIÓN DEL TIPO Y SI ES TRUE O FALSE jeje
    '''
    
    
    if regresion == False and rtype!=int: #check if user want regresion
        pass
    else:
        #---------- GET REG FUNNTION
        dytype=type(dy) #it's used several times
        if regresion == None or regresion == True and rtype!=int: #True for backward working
                #user wants regresion but haven't specified which one
                #automatic choosing depending on if dy values
                #are a single number or an array-like object (linear or
                #ponderated_linear)
            
            if dytype==ndarray or dytype==list or dytype==tuple:
                fr=linear_ponderated_regresion(x,y,dy,None)
            
            elif dytype==str: #not specified
                fr= linear_regresion(x,y,None)
            else: #it's a number
                fr= linear_regresion(x,y,None)
        
        elif isinstance(regresion, (int, long, float, complex)): #it's a number
            if _regresiones[regresion]==linear_ponderated_regresion:
                assert type(dy)==ndarray or type(dy)==list or type(dy)==tuple,\
                "If linear ponderated regresion is choosen, dy must be specified as an array-like object. Otherwise use normal linear regresion" 
                assert size(dy)==size(y), "y and dy must have the same number of elements"
                
                fr=linear_ponderated_regresion(x,y,dy,None)
            else:
                fr=_regresiones[regresion](x,y,None)
            '''
            Asking func dict. As far as my imagination isn't enough and
            linear_ponderated takes an aditionar parameter, the first
            if checks isn't linear_poderated.
            '''
        
        else:# rtype == function:  #user specifies the function
            fr=regresion
            
            
        
        
        
        
        #---------- EXTRAPOLE OPTIONS
        if extrapole: #True
            xrmin,xrmax,_,_=plt.axis() #get actual axis limits
            plt.gca().autoscale(enable=False, axis='both')
        else: #False
            xrmin=min(x)
            xrmax=max(x)
        
        
        
        #---------- GETTING X POINTS
        assert type(np)==int, "np must be an int"
        xr= linspace(xrmin,xrmax,np)
        
        #---------- PLOT REG   
        
        rgr=plt.plot(xr,fr(xr),fmtr,alpha=0.3) #plot regresion 
        '''
        alpha parameter makeS it transparent so it looks better.
        '''
        return gr,rgr
    

    return gr








############

def skip_value (x,value=None):
    '''
    Returns an array from an array-like object without some specified
    values.
    
    It's a copy of the ifc function skip_value
    
    This func will be updated soon as it is non vectorize and it's use
    may ralentize the exec of your program.
    
        PARAMENTERS:
            x, narray-like: array in which we want to skip some key
                            values
            
            value, str: value we want to skip. (None by default)
        
        RETURNS:
            narray: copy (not slice) of elements without the choosen 
                    values
    '''
    x = asarray(x)
    shape = x.shape #storing original shape
    x.shape = (x.size) #reshaping
    l=[] #empty list, declare
    for i in range(x.size):
        if x[i] != value : l.append(i)
    
    return x[l]
