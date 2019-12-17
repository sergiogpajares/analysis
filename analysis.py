# -*- coding:utf8 -*-
#
#  ------> AUTHOR INFO <------------------------------------------------
#
#  analysis module
#  Author: Sergio García Pajares
#  Mail me  'uo272591@uniovi.es' or 'sergiogarciapajares@gmail.com'
#  last update: 17-12-2019
#
#  ------> LICENSE <----------------------------------------------------
#  
#  Copyright 2019 Sergio G. Pajares <sergio@sergio-linux>
#  
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#  
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#  
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#  MA 02110-1301, USA.

#=======================================================================
#==== INFO =============================================================
#=======================================================================
'''
This is an analysis module developed for TEI (Experimetal Techniques I)
an mandatory course of Physics Degree at University of Oviedo.

    Contains:
     - Regresion tools that gets the regresion coeficients
     - Some analytic basic function to work with raw labdata
     - Plot special funtion that manage errorbars and regresion
     - Some OS manage functions to automatize some common processes
     - Some physical constants
     
    
    Dependencies:
     - NumPy
     - SciPy
     - Matplotlib
     
     - OS


Author: Sergio García Pajares
Last update: 17-12-2019

[WARNING] This version is an old version and API will change soon to 
          implement POO. 

'''
__all__ = [
'ponderated_mean','data_plot','linear_regresion','linear_ponderated_regresion',
'linear_origin_regresion','series_ponderated_mean','skip_value','funplot',
'newdirectory','fun3plot',
'c','e','e0','u0'
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
import numpy as np 
import matplotlib.pyplot as plt
#from ifc import skip_value
from scipy.stats import itemfreq

#=======================================================================
#==== DATA =============================================================
#=======================================================================

#####################
# --- PHYSICAL CONSTANTS ---
#

c=299792458
'''
speed of light in vacuum (m/s)
'''

e=1.6021766208E-19
'''
electron charge magnitude (C)
'''

e0=8.854187817E-12
'''
permittivity of free space (N/m)
'''

u0=12.566370614E-7
'''
permeability of free space (N/A²)
'''

#####################
# --- UNITS MANAGEMENT ---
#

#The tranform everything to SI unitis
'''
m=1

c = m/1E2
m = m/1E3
u = m/1E6
n = m/1E9  
'''
#=======================================================================
#==== FUNCTIONS ========================================================
#=======================================================================
#####################
# --- OPERATIVE SYSTEM ---
#
def newdirectory (directory_name, cwd = None):
    '''
    This func is thought to create new directories considering that info
    could be overwriten if we choose an existing path. For example we have
    used this dic in a privious execution. It will ask the user if the path
    exist. If it doesn't will create it.
    
    This func should work either on Linux or Windows.
    
    
        PARAMETERS
            directory_name, str: Is the name of the new directory
            cwd = None, str: Is the path in which we want to create the
                             new directory. If not provided, the current
                             directory at which we are executing the file
                             will be used
                             
        RETUNRS
            newpath, str: Is the full path to the new directory
    '''
    #--- Check inputs are valid ----------------------------------------
    assert isinstance(directory_name,str), "directory_name must be an string"
    if cwd == None: cwd = os.getcwd()
    else: 
        assert isinstance(cwd,str), "The working directory must be an string"
        assert os.path.exists(cwd), "The working directory must exist"
    
    #--- Create the directory ------------------------------------------
    newpath = os.path.join(cwd,directory_name)
    if os.path.isdir(newpath): #check if the directory we want to create
                               #already exists
        AskAgain = True
        while AskAgain: #ask user if he wants to use the same directory
            answer = input("%s already exists. If you continu all \
the info in %s could be overwritten. Do you want to continue? (y/N): "
%(directory_name,directory_name))
            if answer.upper() == 'Y':
                AskAgain = False
                return(newpath)
            elif answer == 'N' or answer == '':
                exit(0) #end execution
    else: #Directory doesn't exist create it
        os.mkdir(newpath)
        return(newpath)


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
    x=np.asarray(x,dtype=float)
    dx=np.asarray(dx,dtype=float)
    
    w=1/(dx**2)
    sw=np.sum(w)
    swx=np.sum(w*x)
    
    return((swx/sw),(1/np.sqrt(sw)))
    
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
    
    x=np.asarray(x)
    y=np.asarray(y)
    dy=np.asarray(dy)
    
    assert np.len(x.shape) == 1, "x must be a 1d array like object" 
        #unique flatten arrays if axis is not specified
    assert np.len(y.shape)  == 1, "y must be a 1d array like object"
    assert np.len(dy.shape) == 1, "dy must be a 1d array like object"
    assert np.size(x) == np.size (y), "x and y mus have the same number of elemnts"
    assert np.size(y) == np.size(y), "y and dy must have the same number of elements"
    
    
    x_values  = np.unique(x) #get unique values of x
    
    y_values  = np.ones_like(x_values,dtype=float) # create an array in which store y means
    dy_values = np.ones_like(x_values,dtype=float) # create an array in which store dy of means
    
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
    x=np.asarray(x) #turn them into arrays
    y=np.asarray(y)
    
    x=np.skip_value(x) #skip None
    y=np.skip_value(y)
    
    #--- checking -----------------------
    assert np.len(x.shape) == 1, 'x must be a vector'
    assert np.len(y.shape) == 1, 'y must be a vector'
    assert np.size(x) == np.size(y), 'x and y must have the same number of elements'
    
    #--- calculate twice used values ----
    N=np.size(x) #number of elements
    
    sx=np.sum(x) #x sumation
    sx2=np.sum(x*x) #x square sumation
    
    sy=np.sum(y) #y sumation
    
    sxy=np.sum(x*y) # xy sumation
        
    delta=float(N*sx2-sx**2) #common denominator for both paramenteres
    
    
    #--- getting linear coeficients -----
    #ax+b
    a=(N*sxy-(sx*sy))/(delta)
    b=((sx2*sy)-(sx*sxy))/(delta)
    
    #--- getting error ------------------
    sigmay=np.sqrt((1./(N-2))*np.sum((y-b-a*x)**2))
    
    da=sigmay*np.sqrt(N/delta)
    db=sigmay*np.sqrt(sx2/delta)
    
    if f == False: return(a,b,da,db) #normal return
    
    elif f == None: 
        f = lambda x: a*x+b #Define reggresion func
        return(f)
    
    else: #f==True
        f = lambda x: a*x+b #Define regresion func
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
    x=np.asarray(x) #turn them into arrays
    y=np.asarray(y)
    dy=np.asarray(dy)
    
    x=np.skip_value(x) #skip None
    y=np.skip_value(y)
    dy=np.skip_value(dy)
    
    #--- checking -----------------------
    assert np.len(x.shape)== 1, 'x must be a vector'
    assert np.len(y.shape)== 1, 'y must be a vector'
    assert np.len(dy.shape)== 1, 'dy must be a vector'
    
    
    assert np.size(x)==np.size(y), 'x and y must have the same number of elements'
    assert np.size(y)==np.size(dy), 'y and dy must have the same number of elements'
    
    #--- calculate twice used values ----
    N=np.size(x) #number of elements
    w=1/(dy**2)
    
    sw=np.sum(w)
    
    swx=np.sum(w*x) #x sumation
    swx2=np.sum(w*x*x) #x square sumation
    
    wy=w*y
    swy=np.sum(wy) #y sumation
    
    swxy=np.sum(w*x*y) # xy sumation
        
    delta=float(sw*swx2-(swx)**2) #common denominator for both paramenteres
    
    
    #--- getting linear coeficients -----
    #ax+b
    a=(sw*swxy-(swx*swy))/(delta)
    b=(swx2*swy-(swx*swxy))/(delta)
        
    #--- getting error ------------------
    da=np.sqrt(sw/delta)
    db=np.sqrt(swx2/delta)
    
    
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
    x=np.asarray(x) #turn them into arrays
    y=np.asarray(y)
    
    x=np.skip_value(x) #skip None
    y=np.skip_value(y)
    
    #--- checking -----------------------
    assert np.len(x.shape)== 1, 'x must be a vector'
    assert np.len(y.shape)== 1, 'y must be a vector'
    assert np.size(x)==np.size(y), 'x and y must have the same number of elements'
    
    #--- calculate twice used values ----
    N=np.size(x) #number of elements
    
    sx2=np.sum(x*x) #x square sumation
    
    sxy=np.sum(x*y) # xy sumation
    
    
    #--- getting linear coeficients -----
    #ax+b
    a=sxy/sx2
    
    #--- getting error ------------------
    sigmay=np.sqrt((1./(N-1))*np.sum((y-a*x)**2))
    
    da=sigmay/np.sqrt(sx2)
    
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
    x = np.asarray(x,dtype=float) #work with arrays
    y = np.asarray(y,dtype=float)
    
    x = np.skip_value(x) #skip none values
    y = skip_value(y)
    
    assert np.size(x) == np.size(y), "x and y must have the same number of elements"
    
    
    #---ERROR BARS---
    if type(dx)!= str or type(dy)!= str : #Check if one of them is introduced
        grb=plt.errorbar(x,y,xerr=dx,yerr=dy,fmt=fmt,ms=ms,ecolor=ecolor,capsize=5,elinewidth=1,markeredgewidth=1) #isn't labeled
                           # capsize, elinewidth, markeredgewidth formats the error bar
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
        xr= np.linspace(xrmin,xrmax,np)
        
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
    x = np.asarray(x)
    shape = x.shape #storing original shape
    x.shape = (x.size) #reshaping
    l=[] #empty list, declare
    for i in range(x.size):
        if x[i] != value : l.append(i)
    
    return x[l]
