# -*- coding:utf8 -*-
#
#  ------> AUTHOR INFO <------------------------------------------------
#
#  analysis module
#  Author: Sergio García Pajares
#  Mail me  'sergio.garcia.pajares @alumnos.uva.es' or 
#           'sergiogarciapajares@gmail.com'
#  last update: 19-09-2020
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
This is an analysis module was originally developed for TEI (Experimetal
Techniques I), an mandatory course of Physics Degree at University of
Oviedo. But it's development was kept as a hobby for his author. 

    Contains:
     - Regresion tools that gets the fitting coeficients
     - Some analytic basic function to work with raw labdata
     - Plot special funtion that manage errorbars and fitting
     - Some OS manage functions to automatize some common processes
     - Some physical constants
     
    
    Dependencies:
     - NumPy
     - SciPy
     - Matplotlib
     
     - OS


Author: Sergio García Pajares
Last update: 17-12-2019
Copyright: GNU General Public License either version 2 of the License,
or (at your option) any later version.

[INFO]    Examples supposes that analysis has been imported as ana.

'''

__all__ = ['DataPlot','Fit','funplot','fun3plot','ponderated_mean','series_ponderated_mean',
            'newDirectory','autoPathRenamer','skip_value','legend','xrad','yrad','setLatex']

__version__ = '3.0.0'

'''

#=======================================================================
#==== IMPORTS ==========================================================
#=======================================================================
import numpy as np 
import matplotlib
import matplotlib.pyplot as plt
from scipy.stats import itemfreq
from scipy.optimize import curve_fit
import os

#=======================================================================
#==== DATA =============================================================
#=======================================================================


#=======================================================================
#==== FUNCTIONS ========================================================
#=======================================================================
#####################
# --- OPERATIVE SYSTEM ---
#
def newDirectory (directory_name, cwd = None):
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
            answer = input("%s already exists. If you continue all \
the info in %s could be overwritten. Do you want to continue? (y/N): "
%(directory_name,directory_name))
            
            answer = answer.upper()
            if answer == 'Y':
                AskAgain = False
                return(newpath)
                
            elif answer == 'N' or answer == '':
                AskAgain = False
            
            else:
                print("\nSorry, I couldn't understand you")
                
                # --- Change dir option ---
                AskAgain2 = True #second question
                while AskAgain2:
                    answer = input("Do you want to enter a new directory name? (Y/n): ").upper()
                    if answer == 'Y' or answer == '':
                        AskAgain2 = False
                        directory_name = input("\nPlease introduce the new directory name\n  ")
                        newDirectory(directory_name,cwd)
                        
                    elif answer == 'N':
                        exit(0) #end execution
                    
                    else:
                        print("\nSorry, I couldn't understand you")
    
    else: #Directory doesn't exist create it
        os.mkdir(newpath)
        return(newpath)

def autoPathRenamer(path,log=False):
    '''
    Adds an  autonumber to the end of a file to path in case
    the file already exist. In case the file has an extension
    format is expected to be .aaa at the end of the filename.

        PARAMETERS:
            path, str: path to the file
            log = False, bool: print or not in terminal info
                
        
        RETURNS:
            newpath, str
    '''
    if os.path.isfile(path):
        
        if log: print('[WARN] file <{}> already exist'.format(path))
        counter=1

        if path[-4] == '.': #case there's extension
            extension = path[-4:]
            path = path[:-4]+str(counter)+extension
            while os.path.isfile(path):
                counter += 1
                digits = np.ceil(np.log10(counter))
                if log: print('[WARN] file <{}> already exist'.format(path))
                path = path[:-(digits+4)]+str(counter)+extension
            
        else:# case there is no extension
            path = path[:-1]+str(counter)+extension
            while os.path.isfile(path):
                counter += 1
                digits = np.ceil(np.log10(counter))
                if log: print('[WARN] file <{}> already exist'.format(path))
                path = path[:-digits]+str(counter)+extension
        
        if log: print("[INFO] Path has been changed to <{}>".format(path))
    
    return path

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
    
    assert len(x.shape) == 1, "x must be a 1d array like object" 
        #unique flatten arrays if axis is not specified
    assert len(y.shape)  == 1, "y must be a 1d array like object"
    assert len(dy.shape) == 1, "dy must be a 1d array like object"
    assert np.size(x) == np.size (y), "x and y mus have the same number of elemnts"
    assert np.size(y) == np.size(y), "y and dy must have the same number of elements"
    
    
    x_values  = np.unique(x) #get unique values of x
    
    y_values  = np.ones_like(x_values,dtype=float) # create an array in which store y means
    dy_values = np.ones_like(x_values,dtype=float) # create an array in which store dy of means
    
    i = 0 #initialice #it's the counter of x-values elements
    
    for value in x_values:
        indices = np.where (x == value) #get indices of x original array
        y_values[i] , dy_values[i] = ponderated_mean(y[indices],dy[indices]) 
            #calculate for chosen values
        i += 1 #next x_value
    
    
    return x_values , y_values , dy_values
    
    


#####################
# --- REGRESION ---
#

def linear_regresion(x,y,**kwargs):
    '''
    Calculates the linear regresion.
    
        PARAMETERS:
            x, 1D-array-like: x points
            y, 1D-array-like: y points
            
        RETURNS:
            For  y = a x + b 
            
                a,  float: a coeficient
                b,  float: b coeficient
                da, float: a error
                db, float: b error
                
                f, funtion: f(x)=ax+b
        
            ( (a,b) , (da,db) , r2 , f ) tuple
        
        DETAILS:
        For more details about it's meaning and
        calculation see: Introducción al Análisis
        de errores, Johon R. Taylor (Reverté 2014)
                    
         
    '''
    #--- preparing inputs ---------------
    x=np.asarray(x) #turn them into arrays
    y=np.asarray(y)
    
    #x=skip_value(x) #skip None
    #y=skip_value(y)
    
    #--- checking -----------------------
    assert len(x.shape) == 1, 'x must be a vector'
    assert len(y.shape) == 1, 'y must be a vector'
    assert np.size(x) == np.size(y), 'x and y must have the same number of elements'
    
    #--- calculate twice used values ----
    N=np.size(x) #number of elements
    
    sx=np.sum(x) #x sumation
    sx2=np.sum(x*x) #x square sumation
    
    sy=np.sum(y) #y sumation
    sy2=np.sum(y*y) #y square sumation
    sxy=np.sum(x*y) # xy sumation
        
    delta=float(N*sx2-sx**2) #common denominator for both paramenteres
    
    
    #--- getting linear coeficients -----
    #ax+b
    a=(N*sxy-(sx*sy))/(delta)
    b=((sx2*sy)-(sx*sxy))/(delta)
    
    #--- getting error ------------------
    sigmay=np.sqrt((1./(N-2))*np.sum((y-b-a*x)**2))
    
    r2 = ( sxy - (sx*sy/N) )**2 / ( ( sx2-(sx**2/N))*(sy2 - (sy**2/N)) ) #correlation squared coeficient
    
    da=sigmay*np.sqrt(N/delta)
    db=sigmay*np.sqrt(sx2/delta)
    
    f = lambda x: a*x+b #Define regresion func
    return( np.array([a,b]) ,np.array([da,db]) ,r2, f )
    
def linear_ponderated_regresion(x,y,dy,**kwargs):
    '''
    Calculates the linear ponderated regresion.
    
        PARAMETERS:
            x, 1D-array-like: x points
            y, 1D-array-like: y points
            dy, 1D-array-like: y error values
            
         RETURNS:
                    
        For  y = a x + b 
            
           a,  float: a coeficient
           b,  float: b coeficient
           da, float: a error
           db, float: b error
           
           f, funtion: f(x)=ax+b
        
        ( (a,b) , (da,db) , r2 , f )
         
    '''
    #--- preparing inputs ---------------
    x  = np.asarray(x) #turn them into arrays
    y  = np.asarray(y)
    dy = np.asarray(dy)
    
    #x  = skip_value(x) #skip None
    #y  = skip_value(y)
    #dy = skip_value(dy)
    
    #--- checking -----------------------
    assert len(x.shape) == 1, 'x must be a vector'
    assert len(y.shape) == 1, 'y must be a vector'
    assert len(dy.shape) == 1, 'dy must be a vector'
    
    
    assert np.size(x)==np.size(y), 'x and y must have the same number of elements'
    assert np.size(y)==np.size(dy), 'y and dy must have the same number of elements'
    
    #--- calculate twice used values ----
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
    
    r2 = None #<----------------------------------------------------------------------------------------------------
    
    f=lambda x: a*x+b #Define regresion func
    return( np.array([a,b]) , np.array([da,db]), r2 , f )


def linear_origin_regresion(x,y,**kwargs):
    '''
    Calculates the linear regresion.
    
        PARAMETERS:
            x, 1D-array-like: x points
            y, 1D-array-like: y points
            
            
        RETURNS:
                    
        For  y = a x
            
           a,  float: a coeficient
           da, float: a error
           
           f, funtion: f(x)=ax+b
           
        ( (a) , (da) , r , f )
         
    '''
    #--- preparing inputs ---------------
    x=np.asarray(x) #turn them into arrays
    y=np.asarray(y)
    
    #x=skip_value(x) #skip None
    #y=skip_value(y)
    
    #--- checking -----------------------
    assert len(x.shape)== 1, 'x must be a vector'
    assert len(y.shape)== 1, 'y must be a vector'
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
    r2 = None #<-----------------------------------------------------------------------------------------------------
    
    f=lambda x: a*x #Define regresion func
    return( np.array([a]) , np.array([da]) , r2 , f)

def auto_linear(x,y,dx,dy):
    '''
    Choose which regresion between linear and linear ponderated is 
    needed considering the type of dy.
    
    If dy is an array like 
    
        PARAMETERS:
            x, 1D-array-like: x points
            y, 1D-array-like: y points
            dx, number or array-like: x error
            dy, number or array-like: y error
                    
        RETURNS:
        For  y = a x + b 
            
           a,  float: a coeficient
           b,  float: b coeficient
           da, float: a error
           db, float: b error
           
           f, funtion: f(x)=ax+b
        
        ( (a,b) , (da,db) , r2 , f )
        
         
    '''
    if len(dy)==1 or len(dy)==0: #It's a number
        print("[INFO] linear regresion used")
        return linear_regresion(x,y)
    else: #It's an array like object
        print("[INFO] linear ponderated regresion used")
        return linear_ponderated_regresion(x,y,dy)
    


def custom_fitting (x,y,dy,func):
    '''
    Automated calling of scipy.optimeze.curve_fit
    '''
    if dy == []: dy = None #solve some implementation issue due to diference in DataPlot dy=[]
                           #default argument and curve_fitting sigma = None expected
    p, pcov, infodict, _, _ = curve_fit(func, x, y,sigma=dy, full_output=True)
        #pcov is the covariance matrix of the parameters
    dp = np.sqrt(np.diag(pcov)) # errores estándar de los parámetros

    return (p,dp,1-(infodict['fvec']**2).sum()/((y-y.mean())**2).sum(),lambda t: func(t,*p))

def noFit(**kwargs):
    '''
    No fitting function to allow use of ref=False in DataPlot func
    '''
    return ((),(),None,None)


def sinusoidal(x,y,dy):

    return custom_fitting(x,y,dy,lambda t, A, phi, B: A*np.sin(x + phi) + B)



# ======================== FITTING CLASS ============================= #
class Fit (object):
    '''
    This class aims to make easier fitting problems solution.
    It's a wrapper around all fitting functions for easier use.
    '''
    regfuncs={
        #This dict is used by fit class builder to calc regresion from
        #the different regresion functions. Add here a regresion function
        #and in fit.dict and fit.help and it will be implemented.
        False              : (noFit,('No fitting',[],[])),
        1                  : (auto_linear,('f(x)=ax+b',['a','b'] , ['da','db'])),
        'auto_linear'      : (auto_linear,('f(x)=ax+b, auto linear',['a','b'] , ['da','db'])),
        2                  : (linear_regresion, ('f(x)=ax+b',['a','b'] , ['da','db'])),
        'linear'           : (linear_regresion, ('f(x)=ax+b, linear',['a','b'] , ['da','db'])),
        3                  : (linear_ponderated_regresion,('f(x)=ax+b',['a','b'] , ['da','db'])),
        'linear_ponderated': (linear_ponderated_regresion,('f(x)=ax+b, linear ponderated',['a','b'] , ['da','db'])),
        4                  : (linear_origin_regresion,('f(x)=ax',['a'],['da'])),
        'linear_origin'    : (linear_origin_regresion,('f(x)=ax',['a'],['da'])),
        5                  : (sinusoidal,('f(x)=Asin(x+phi)+B',['A','phi','B'],['dA','dphi','dB'])),
        'sinusoidal'       : (sinusoidal,('f(x)=Asin(x+phi)+B',['A','phi','B'],['dA','dphi','dB']))
    }

    def __init__ (self,x,y,dx=[],dy=[],reg=True):
        '''
            PARAMETERS  
                x, array-like:
                y, array-like:
                
                OPTIONAL
                
                dx, number or array-like:
                dy, number or array-like:
                
                reg: type of fitting (see below)
                
                
            ATTRIBUTES
                p, numpy.array: array containing the fitting parameters
                dp, numpy.array: array containing the estimated error on
                                 fitting paramenters
                r2, float: squared correlation coeficient R²
                f, lambda.function: function f(x) representing the
                                    fitted curve.
                type: reg type provided by user
                data, list [x,y,dx,dy]: original data provided by the user
                dict, dictionary: This dictionary provides a user-friendly
                      way to acces all info stored by the user. Including
                      parameters in a human redable way.

                      KEYS:
                        x : x original data
                        dx: x error original data
                        y : y original data
                        dy: y error original data
                        f : fitted lambda.function
                        r2: squared correlation coeficient R²

                        Also all the keys of the parameters depending
                        on the type of reg. To check such parameters
                        see the regfuncs dictionary

            
            REGRESION POSIBILITIES
            
              =====================================================
              |               REGRESION POSIBILITIES              |
              |===================================================|
              |             KEY          |          TYPE          |
              |---+-----------------------------------------------|
              | 0 | False                | No regresion drawn     |
              |---+----------------------+------------------------|
              | 1 | auto_linear          | linear or linear pon-  |
              |   |                      | derated depending on   |
              |   |                      | dy type (number or     |
              |   |                      | array like)            |
              |---+----------------------+------------------------|
              | 2 | linear               | linear                 |
              |---+----------------------+------------------------|
              | 3 | linear_ponderated    | linear ponderated      |
              |---+----------------------+------------------------|
              | 4 | linear_origin        | linear crossing origin |
              |---+----------------------+------------------------|
              | 5 | sinusoidal           | sinusoidal             |
              |---+----------------------+------------------------|
              |   |                      |                        |
              +--------------------------+------------------------|
              | Function f(x) specififed by the user              |
              =====================================================

              Note: In case of user specified function format must
              be  f(x,params) where params are the value we want to
              estimate
            
            EXAMPLES
                >>> import numpy, analysis
                >>> x = [1,2,3,4]
                >>> y = [2.8,3.9,6.1,7.9]
                >>> 
                >>> dx = [.02,.15,.10,.05]
                >>> dy = [.9,.15,.2,.3]
                >>> 
                >>> myfit = analysis.fit(x,y,dx,dy,reg='linear_ponderated')
                >>> print(myfit.p) #show parameters
                       numpy.array([ 2.0101029  , -0.05116932 ])
                >>> print(myfit.dp) #show error on parameters
                       numpy.array([ 0.14936723 ,  0.39837133 ])
                >>> print(myfit.dict['a']," +- ",myfit.dict['da'])
                       2.0101029 +- 0.14936723
        '''
        if reg == False or reg == 'no':
            #second assertion is to avoid interpreting 0 as False
            self.p , self.dp , self.r2, self.f = None, None, None, None
            self.type = reg
            self.data = np.array([x,y,dx,dy])
            
        #--- Prepare data ---
        x = np.asarray(x,dtype=float) #work with arrays
        y = np.asarray(y,dtype=float)
    
        assert np.size(x) == np.size(y), "x and y must have the same number of elements"
        
        
        #--- Get regresion ---
        if isinstance(reg,(int,str)):
            self.p , self.dp , self.r2 ,self.f = self.regfuncs[reg][0](x=x,y=y,dx=dx,dy=dy)
        else: #Case of custom function instead of reg type
            self.p , self.dp , self.r2 ,self.f = custom_fitting(x,y,dy,func=reg)
        
         # <- WHERE ->
         #  self.p: params of fitting
         #  self.dp: error of params of fitting
         #  self.r2 square of correlation coeficient
         #  self.f lambda func containing the fitted function
        
        #--- Customize instance ---
        self.type = reg #reg type
        self.data = [x,y,dx,dy] #original data provided

        #--- Dictionary ---
        self.dict = {
            'x'   : self.data[0],
            'y'   : self.data[1],
            'dx'  : self.data[2],
            'dy'  : self.data[3],
            'type': self.type, #chage to human readable
            'r2'  : self.r2,
            'f'   : self.f,
        }

        #Defining parameters and error of parameters in terms of parameters name
        #stored in regfuncs dictinary
        human_names = self.regfuncs.get(self.type)#returns None if key is not avaible,
                                                  # the case of user specified function
        if human_names == None: #case user define function
            for i in np.arange(len(self.p)):
                self.dict.update(
                    {
                        'p'+str(i+1)  : self.p[i],
                        'dp'+str(i+1): self.p[i]
                    }
                )
        else:
            human_names = human_names[1] #choose tuple of names
            # params
            for i in np.arange(len(human_names[1])):
                self.dict.update( {human_names[1][i] : self.p[i]} )

            #errors
            for i in np.arange(len(human_names[2])):
                self.dict.update( {human_names[2][i]:self.dp[i]} )
        
    def help (self):
        '''
        Prints info about the params
        '''

        try:
            print(self.regfuncs[self.type][1][0]) #reimplement using, key is not in dicctionary
        except:
            print("Reg function was especified manually")

    def GetType (self):
        '''
        Return human readable type of regresion
        '''



#####################
# --- PLOTTING ---
#
def setLatex(mode='pdf'):
    '''
    Wrapper aroud matplotlib to use LaTeX rendering configuration
    
    IT MUST BE CALLED BEFORE ANY PLOTTING! I recommend to call it
    just after importing matplotlib or matplotlib.pyplot

    PARAMETERS
        mode=pdf: it might be 'pdf' or 'png' depending on what
           plt.savefig() we will like to use.

           Note: pdf mode will also allow you to use plt.show()

    '''

    if mode=='pdf':
        matplotlib.rcParams['text.usetex'] = True
        matplotlib.rcParams['backend'] = 'pgf'
    elif mode=='png':
        matplotlib.rcParams['text.usetex'] = True
        matplotlib.rcParams['backend'] = 'TkAgg'
    else: raise ValueError("type must be 'pdf' or 'png'")

def funplot (f,xmin,xmax,n=100,fmt='',legend='',title='',xlabel='',
ylabel='',label='',adjust=False,**aditional_plot_params):
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

            adtional_plot_params: all params accepted by
                                  matplotlib.pyplot.plot()
        
        RETURNS
            A list of Line2D objects representing the plotted data.
            It's generated by matplotlib.pyplot.plot() func.
            
    '''
    #plotting
    x  = np.linspace(xmin,xmax,n) #create x
    y  = f(x) #eval func. Values will be used later.
    
    if fmt != '': #manage fmt provided or not by user
        gr = plt.plot(x,y,fmt,label=label,**aditional_plot_params) #plot
    else:
        gr = plt.plot(x,y,    label=label,**aditional_plot_params) #plot
    
    #personalization
    if legend != '' : plt.legend() #create legend
    if title != '' : plt.title('u'+title)
    if xlabel != '' : plt.xlabel('u'+xlabel)
    if ylabel != '' : plt.ylabel(ylabel)
    if adjust == True : plt.axis([xmin,xmax,np.min(y),np.max(y)])
    
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

class DataPlot (Fit):
    '''
    This class aims to make easier solving data plotting with 
    (or without) fitting and it's plotting. It's used is focused on
    pairs of data (x ~ y). It can also manage errorbar plotting.
    '''
    def __init__ (self,x,y,fmt='none',fmtr='none',reg=False,dx=[],dy=[],n=300,ax=None,**aditional_params): #quitar ax
        '''
        

            PARAMETERS
                x, array-like: x data to plot and fit
                y, array-like: y data to plot and fit
                fmt='none', str: format of points using 
                                 matploplotlib.pyplot.plot()
                fmtr='none',str: format of regresion line using
                                 matploplotlib.pyplot.plot()
                reg=False, number str or callable: type of regresion
                        For types see below.
                dx=[], number or array-like: standart errors
                        of x-values drawing it's error bar
                dy=[], number or array-like: standart errors
                        of y-values for it's use in fitting and
                        drawing of errorbars
                n=300, int: number of points used to draw the fitted
                            function
                ax=None, matplotlib.axes: axes in wich points and fitted
                                          function is drawn.
                                          Default value will use
                                          matplotlib.pyplot.gca()
                **aditional_params: every parameter accpeted by 
                                    matplotlib.pyplot.plot or 
                                    matplotlib.pyplot.plot depending on
                                    errorbars provided or not.
            
            ATRIBUTES
                p, numpy.array: array containing the fitting parameters
                dp, numpy.array: array containing the estimated error on
                                 fitting paramenters
                r2, float: squared correlation coeficient R²
                f, lambda.function: function f(x) representing the
                                    fitted curve.
                ax: axes in wich plot is made
                gr: points plot
                rgr: fitted curve plot
                type: reg type provided by user
                data, list [x,y,dx,dy]: original data provided by the user
                dict, dictionary: This dictionary provides a user-friendly
                      way to acces all info stored by the user. Including
                      parameters in a human redable way.

                      KEYS:
                        x : x original data
                        dx: x error original data
                        y : y original data
                        dy: y error original data
                        f : fitted lambda.function
                        r2: squared correlation coeficient R²

                        Also all the keys of the parameters depending
                        on the type of reg. To check such parameters
                        see the regfuncs dictionary

                        Note: in case user defined reg function the key
                        of paramters will be  p1 , p2 , ... ,pn and 
                        dp1 , dp2 , ... , dpn

            
            REGRESION POSIBILITIES
            
              =====================================================
              |               REGRESION POSIBILITIES              |
              |===================================================|
              |             KEY          |          TYPE          |
              |---+-----------------------------------------------|
              | 0 | False                | No regresion drawn     |
              |---+----------------------+------------------------|
              | 1 | auto_linear          | linear or linear pon-  |
              |   |                      | derated depending on   |
              |   |                      | dy type (number or     |
              |   |                      | array like)            |
              |---+----------------------+------------------------|
              | 2 | linear               | linear                 |
              |---+----------------------+------------------------|
              | 3 | linear_ponderated    | linear ponderated      |
              |---+----------------------+------------------------|
              | 4 | linear_origin        | linear crossing origin |
              |---+----------------------+------------------------|
              | 5 | sinusoidal           | sinusoidal             |
              |---+----------------------+------------------------|
              |   |                      |                        |
              +--------------------------+------------------------|
              | Function f(x) specififed by the user              |
              =====================================================

              Note: In case of user specified function format must
              be  f(x,params) where params are the value we want to
              estimate

                

                
        '''
        
        # -- PREPARE INPUT --
        x = np.asarray(x,dtype=float) #work with arrays
        y = np.asarray(y,dtype=float)

        assert np.size(x) == np.size(y), "x and y must have the same number of elements"
        if reg in (1,'auto_linear',2,'linear',3,'linear_ponderated',4,'linear_origin'): n = 2 #simplify in linear cases
        
        
        if ax == None: self.ax = plt.gca()
        else: self.ax = ax
        #else: assert isinstance(ax,matplotlib.axes) , "ax is not a valid axis"
        
        #change default matplotlib parmaters in case not porvided
        nondefaulterrorparams = [('ms'           , 3 ), #marker size
                                 ('ecolor'       ,'k'), #errorbar color
                                 ('capsize'        ,5), #errorbar cap size
                                 ('elinewidth'     ,1), #errorbar line width
                                 ('markeredgewidth',1)] #error bar cap width
        
        # -- POINTS PLOTTING --
        if len(dx)!=0  or len(dy)!= 0 : #Check if one of them is introduced
            for key , value in nondefaulterrorparams:
                if not key in aditional_params: # if a new value has not 
                                                # been provided add it
                    aditional_params.update( { key : value } )
        
            self.gr = plt.errorbar(x,y,xerr=dx,yerr=dy,fmt=fmt,**aditional_params)
            
        else:
            self.gr = plt.plot(x,y,fmt,**aditional_params)
        
        # -- REGRESION --
        super().__init__(x,y,dx,dy,reg)
                #make regresion from parent class and get all
                #parameters set up for regresion

        sample = np.linspace(np.min(x),np.max(x),n,endpoint=True)
        self.rgr = plt.plot(sample,self.f(sample),fmtr,alpha=.3)
            
def xrad(ax=None):#multiples=np.pi/3,ax=None):
    '''
    Change x axis thicks to radians in terms of pi

        PARAMETERS:
            ax=None: axes in wich change the labels
            Default parameters uses matplotlib.pyplot.gca()
    '''
    if ax == None: ax=plt.gca()
    ax.xaxis.set_major_locator(plt.MultipleLocator(np.pi / 2))
    ax.xaxis.set_minor_locator(plt.MultipleLocator(np.pi / 8))
    ax.xaxis.set_major_formatter(plt.FuncFormatter(format_func))   
        
        
def yrad(ax=None):#multiples=np.pi/3,ax=None):
    '''
    Change y axis thicks to radians in terms of pi

        PARAMETERS:
            ax=None: axes in wich change the labels
            Default parameters uses matplotlib.pyplot.gca()
    '''
    if ax == None: ax=plt.gca()
    ax.yaxis.set_major_locator(plt.MultipleLocator(np.pi / 2))
    ax.yaxis.set_minor_locator(plt.MultipleLocator(np.pi / 8))
    ax.yaxis.set_major_formatter(plt.FuncFormatter(format_func))

def format_func(value, tick_number):
    '''
    Auxiliar function for yrad and xrad
    '''
    # find number of multiples of pi/2
    N = int(round(2 * value / np.pi))
    if N == 0:
        return "0"
    elif N == 1:
        return r"$\pi/2$"
    elif N == 2:
        return r"$\pi$"
    elif N % 2 > 0:
        return r"$%d\pi/2$"%(N)
    else:
        return r"$%d\pi$"%(N // 2)

def legend(ax=None):
        '''
        Show legend without displaying error bars

            PARAMETERS:
                ax, axes: axes in wich to draw legend
                    If not specified matplotlib.pyplot.gca()
                    is used.
            RETURNS:
                legend object
        '''
        #remove error bar from legend
        if ax == None: ax = plt.gca()
        handles, labels = ax.get_legend_handles_labels()
        handles = [handle[0] for handle in handles]

        return ax.legend(handles, labels)


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
    x.shape = (x.size) #reshaping
    l=[] #empty list, declare
    for i in range(x.size):
        if x[i] != value : l.append(i)
    

    return x[l]

