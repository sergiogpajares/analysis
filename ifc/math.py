from numpy import *
'''
last update: 10-03-2019
'''

__all__=['segundo_grado','finite_diff','integral',
'integral_funciones_lineales','matriz_nxn_1an2','maxmin','mr','norm',
'raiz_cubica','segundo_grado','volumen_esfera','volumen_cilindro_recto',
'volumen_paralelepipedo_recto']

def segundo_grado(a,b,c):
    '''
    Parametros (a,b,c) ax2+bx+c=0.
    
    Devuleve las soluciones reales o complejas de la ecuacion de
    segundo grado en la forma (x1,x2). Si es una ecuacion de primer 
    grado devuleve la tupla (x1,None)
    '''
    '''
    >>>> CORREGIR <<<<
    Es preferible que siempre devuelva el mismo numero de valores la 
    funcion por lo que se introduce el None, pero no es del todo reco-
    mendable. Otra de las posiblidades es asignar la solucion,
    a una varible como una tupla de dos soluciones. Y usar en la impre-
    sion un if que decida si es o no un tuple.
    
    Por otro lado, es recomendable usar el minimo numero de returns.
    '''
    
    #ECUACION PRIMER GRADO
    if a==0:
        x1=-c/b
        x2=None

    #ECUACION SEGUNDO GRADO
    else :
        #CALCULO DISCRIMINANTE
        d=b**2-4*a*b
	
        #CALCULO SOLUCION REAL
        if d>=0:
            x1=(-b+d**.5)/(2*a)
            x2=(-b-d**.5)/(2*a)
	
        #CALCULO SOLUCION COMPLEJA
        else :
            x1=(-b+complex(d)**.5)/(2*a)
            x2=(-b-complex(d)**.5)/(2*a)

    return x1,x2

def finite_diff (f, x, h = 0.0001):
    '''
    Aproxima la primera y segunda derivada
    de f en x; h tiene un valor por defecto 0.0001
    '''

    fm, f0, fp = f(x-h), f(x), f(x+h)
    df, ddf = (fp-fm)/(2.*h), (fp-2.*f0+fm)/h**2.
    return df, ddf

def mr (funcion_posicion,t):
    '''
    Devuelve (posicion, velocidad, aceleracion) de un movil
    en movimiento rectilineo en un instante t de tiempo. Se pasa
    como parametro (funcion de posicion en funcion de t, tiempo).
    '''
    
    v,a=finite_diff(f,t) #Derivamos numericamente la posicion, la pri-
                         #mera derivada corresponde a la velocidad
                         #la segunda a la aceleracion.

    x=funcion_posicion(t) #evalua la funcion posicion en t
    return  x,v,a

def maxmin (f, a, b, n=1001):
    '''
    Parametros:
    (funcion, extremo izq intervalo, extremo drch intervalo,
      *numero de pasos)
      
    Devuleve:
    (max abs(f(x)), min abs(f(x)) (si hay dos identicos el primero)
    
    Lo calcula subdividiendo en n divisiones el intervalo, evaluando
    en todos los puntos y devolviendo el mayor de los puntos.
    '''
    
    h=(b-a)/float(n-1) #el 'paso'
    X=[a+i+h for i in range(n)] #lista de valores x a evaluar
    Y=[f(x) for x in X] #lista de valores de f(x)
    return max(Y), min(Y)

def integral_funciones_lineales(f,a,b):
    '''
    integral(f, a, b)
    
    Devuelve el valor de la integral de f entre a y b.
    
    Usa primera aproximacion en trapecio entre el principio y el final.
       b             b - a  
      / f(x) dx =~  ------- ( f(a)+f(b) ) 
     a                 2    
    '''
    return (b-a)/2*(f(a)+f(b))

def integral(f,a,b,n=10000):
    '''
    integral(f,a,b,*n) 
    
    Devuleve la integral de f(x)dx entre a y b calculando n puntos
    intermedios. n por defecto= 10 000
    '''
    X=[]#Lista puntos X a completar en bucle
    deltaX=(b-a)/n #Paso del intervalo
    
    for i in range(n+1):#Puntos x en los que evaluar la funcion
        X.append(a) #se reutiliza la variable a. Pasa a significar
                    #valor Xi
        a+=deltaX
        
    Y=[f(i) for i in X] #Puntos y de la funcion
    I=0 #valor de la integral
    for i in range(1,n+1):#Sumatorio de todos los trapecios.(n+1)=len(X)
        I+=(X[i]-X[i-1])*(Y[i-1]+Y[i])*.5 #Area de cada trapecio.
        '''
        El area de cada trapecio (entre X1,X2) se calcula como:
        I+=  Y2(X2-X1)    -   .5(Y2-Y1)(X2-X1)
        Area rectangulo      Area triangulo
        
        Y2 |                Y2 |   __          Y2 |
           |    /|             |  |  |            |   /|
        Y1 |   / |       =  Y1 |  |  |      -  Y1 |  /_|
           |___|_|_____        |__|__|____        |________
              X1 X2               X1 X2             X1 X2
                
            Sacando factor comun y operando se obtiene:
            I+=(X2-X1)(Y1+Y2)*.5
        '''
        '''
        OTRA POSIBILIDAD (mas optimizada):
        h=(a+b)/n #Paso del interavlo
        k=0 #
        for in in range(1,n): #sumatorio desde 1 hasta n
            k+=f(a+i*h)
        return (h*(f(a)+f(b))/2)+k)
        '''
    del X, Y #borramos listas inecesarias
    return I

def matriz_nxn_1an2 (n):
    '''
    parametros: (n9
    
    Devuelve una matriz n x n con valores desde 1,---,n^2.
    
    Ejemplo n=3       Desarrollo           
            
    |  1  2  3  |   |   1    2   ---  n  |
    |  4  5  6  |   |  n+1  n+2  ---  2n |  
    |  7  8  9  |   | 2n+1  2n+2 --- n^2 | 
    

    '''
    #Definimos lista de filas de la matriz
    
    return [range(i*n+1,(i+1)*n+1) for i in range(n)]
    '''
    La matriz es una lista de filas, que a su vez son listas. 
        
    Se deduce una formula que exprese los valores de cada fila en 
    fucion del numero de fila i. Los valores de la fila i
    valen:|(i*n)+1 --- (i+1)*n| donde i=0,---,n-1
        
    Ejemplo n=3       Desarrollo             Termino general
        
    |  1  2  3  |   |   1    2   ---  n  |     fila i
    |  4  5  6  |   |  n+1  n+2  --- n+n |  | in+1 in+2 --- (i+1)n | 
    |  7  8  9  |   | 2n+1  2n+2 --- 3n  | 
        
    Notese que las filas se empiezan a contar por el cero hasta la
    fila n-1 (total de n filas).
    
    La fila se genera en lista comprimida siguiendo este termino
    general y se anade a la lista. El +1 que se anade al cierre 
    de la lista es para incluir el ultimo elemento (intervalo 
    abierto).
    
    El bucle recorre desde la primera fila (i=0), hasta la ultima
    (i=n-1).
    '''

def volumen_esfera(r):
    '''
    Parametros: (r)
    Devuelve el volumen de la esfera de radio r
    '''
    return r*r*r*pi*4/3.

def volumen_cilindro_recto(r,h):
    '''
    Parametros: (r,h)
    Devuelve el volumen del cilindro recto de base con radio r 
    y altura h.
    '''
    return r*r*h*pi

def volumen_paralelepipedo_recto(a,b,c):
    '''
    Parametros(a,b,c)
    Devuelve el volumen de un paralalepipedo de dimensiones a,b,c.
    '''
    return a*b*c

def raiz_cubica(x,EPS=1.00E-06):
    '''
    Parametros (x, EPS)
    
    Devuelve la raiz cubica de un numero x con una precision mayor que
    EPS. Si EPS no se especifica se tomara un valor por defecto de
    1.00E-06
    '''
    x1=1. #valor de x1 inicial
    coeficiente=1./3 #Usar la maxima precision del ordenador
    x2=10. #Declaracion de variable
    while abs(x2-x1)>EPS:
        x1=x2 #next loop
        x2=coeficiente*(2*x1+(x/(x1*x1)))
    return x2

def euclidean_distance (A,B):
    '''
    Calaculates the ecludean distance between two n-dimensional vectors.
    Vectors can be calculated together if they are gathered as:
    
      A=[[a1,a2, ... ,an], <-vec 1   B=[[y1,y2, ... ,yn],
         [b1,b2, ... ,bn], <-vec 2      [z1,z2, ... ,zn],
           .                              .
           .                              .
           .                              .
    
    PARAMETERS:
        A: narray or list,  eucludian coordinates of A
        B: narray or list,  euclidian coordinates of B
    
    RETURNS:
        array,  euclidean distance between A and B
         array([distance([a1,...,an],[y1,...,yn])],
               [distance([b1,...,bn],[z1,...,zn])],
                 .
                 .
                 .
    '''
    A=asarray(A) #if the input is not an array, convert it into one
    B=asarray(B)
    
    assert len(A.shape)< 3,"A can't be bigger than a matrix of vectors"
    assert len(B.shape)< 3,"B can't be bigger than a matrix of vectors"
    
    #ADAPTING FORMAT
    if len(A.shape)==1 and len(B.shape)==1:
        pass
        
    else:# len(A.shape)==2 and len(B.shape)==1:
        '''
        While comparing a vector with a matrix of vector we need to
        turn the vector into a 1 row matrix.
        '''
        if len(A.shape)==1:
            A.shape=(1,size(A))#turning into matrix
        else:# len(B.shape)==1:
            B.shape=(1,size(B))#turning into matrix

        assert len(A[0,:])==len(B[0,:]),"vectors inside the matrix must \
have the same number of dimendions"
    
    #GETTING DISTANCE
    distance=sqrt(sum(((A-B)**2),axis=1)) #axis for vertical sum #make calculus
    distance.shape=(size(distance),1) #reshape for vertical output
    return  distance

def norm(A):
    '''
    Returns the norm of an n-dimensional vector.
    
    PARAMETERS:
        A: 1d-array or list,  vector
    
    RETURNS:
        number,  vector's norm
    
    '''
    A=asarray(A) #if the input is not array, convert it into one
    
    assert len(A.shape) == 1, 'A must be an 1d-array or list'
    
    return sqrt(sum((A)**2))


