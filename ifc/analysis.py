#Analysis submodule
#Author: Sergio Garcia Pajares
#documentation language: engish
#Last update: 2019-02-22
#version: 0.0
#
#  version history
#
#




'''
This submodule contains the tools prepragramed for data analysis both
in IFC and TEI for academical propouses at University of Oviedo.

Last update: 22-02-2019
Version: 0.0
'''

#=======================================================================
#======   IMPORTS   ====================================================
#=======================================================================

from numpy import *
import matplotlib.pyplot as plt


#=======================================================================
#======   DEFINITIONS   ================================================
#=======================================================================

def funplot (f,xmin,xmax,n=100,fmt='b-',legend='',title='',xlabel='',ylabel='',adjust=False):
    #SIN ACABAR 
    '''
    Plots an R-->R fuction 
    
        PARAMETERS:
            f, function: R-->R function
            xmin, number: lower limit in x axis
            xmax, number: upper limit in x axis
            
            OPTIONAL
            n, int: number of x divisions for x. 
            fmt, str: line format
            legend, str: name of function for legend
            xlabel, str: xlabel
            ylabel, str: tlabel
            adjust, bool: for adjusting axis to xmin and xmax limits.
    '''
    #plotting
    x=np.linspace(xmin,xmax,n) #create x
    y=f(x) #calculate y
    gr=plt.plot(x,y,fmt,label=leyenda) #plot and evaluate
    
    #personalization
    if legend != '' : plt.legend() #create legend
    if title != '' : plt.title(u+title)
    if xlabel != '' : plt.xlabel(u+xlabel)
    if ylabel != '' : plt.ylabel(ylabel)
    if adjust==True : ptl.axis([min(x),max(x),min(y),max(y)])
    
    return gr 

def linear_regresion(x,y):
    '''
    Calculates the linear regresion.
    
        PARAMETERS:
            x, 1D-array-like: x points
            y, 1D-array-like: y points
        
        RETURNS:
            (a,b,da,db) tuple
        
            For  y = a x + b 
        
            a,  float: a coeficient
            b,  float: b coeficient
            da, float: a error
            db, float: b error
         
    '''
    #--- preparing inputs ---------------
    x=np.asarray(x) #turn them into arrays
    y=np.asarray(y)
    
    #--- checking -----------------------
    assert len(x.np.shape)== 1, 'x must be a vector'
    assert len(y.np.shape)== 1, 'y must be a vector'
    assert np.size(x)==np.size(y), 'x and y must have the same number of elements'
    
    #--- calculate twice used values ----
    N=size(x) #number of elements
    
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
    sigmay=np.sqrt(1./(N-2))*np.sum((y-b-a*x)**2)
    
    da=sigmay*np.sqrt(N/delta)
    db=sigmay*np.sqrt(sx2/delta)
    
    return(a,b,da,db)

def data_quick_plot(x,y,fmt='bo',label='',dx=None,dy=None,regresion=False,fmtr='b-',a=None,b=None,adjust=False):
    '''
    Plots a pair of data. It can also plot its error bars and linear
    regresion.
    
        PARAMETERS:
            x, 1d array-like: x
            y, 1d array-like: y
            
            OPTIONAL
            fmt='bo', str: format for the points
            label='', str: label for data
            dx=None, 1d array-like: error values for x
            dy=None, 1d array-like: error values for y
            
            -----
            y=ax+b
            regresion=False, bool: Do we want to draw linear regresion?
                                   True = yes   False = no
            fmtr='b-', str: format for regresion line
            a=None, number: a coeficient. If not specified is calculated
            b=None, number: b coeficient. If not specified is calculated
            
            -----
            adjust=False, bool: It adjust graphic axis to data
        
        RETURNS:
            if not regresion drawed:
                points draw object
            
            if regresion drawed:
                tuple (gr,regr)
                    gr:   points draw object
                    regr: regresion draw object
    '''
    #input reading
    x=asarray(x)
    y=asarray(y)
    
    x=skip_value(x) #getting out of None values
    y=skip_value(y)
     
    #normal plotting
    gr=plt.plot(x,y,fmt,label=label)
    
    #error bars plotting
    if dx!=None and dy!= None : gr=plt.err(x,y,xerr=dx,yerr=dy,label=label)
    
    #adjust
    if adjust == True: plt.axis([np.min(x),np.max(x),np.min(y),np.max(y)])
        
    #regresion plotting
    if regresion == True: #represent linear regresion
        if a==None or b==None: a,b,da,db=linear_regresion(x,y) 
            #if coeficientes are not specify, calculate them
        
        fregr= lambda x,a,b : a*x+b #define regresion function
        xregr= array([np.min(x),xp.max(x)])
        
        regr=plt.plot(x,fregr(x,a,b),fmtr) #plot regresion 
        
        return gr,regr

    return gr


def skip_value (x,value=None):
    '''
    Returns an array from an array-like object without some specified
    values
    
        PARAMENTERS:
            x, narray-like: array in which we want to skip some key
                            values
            
            value, str: value we want to skip. (None by default)
        
        RETURNS:
            narray: copy of elements without the choosen values
    '''
    x=asarray(x)
    shape=x.shape #storing original shape
    x.shape=(x.size) #rehaping
    l=[] #empty list
    for i in range(x.size):
        if x[i]!= value : l.append(i)
    return x[l]
