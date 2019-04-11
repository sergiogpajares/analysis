def media_aritmetica(conjunto):
    '''
    parametros(conjunto) donde conjunto es una lista o tupla con todos
    los datos.
    
    Devuelve la media artimetica del conjunto de datos.
    '''
    assert type(conjunto)==list or type(conjunto)==tuple,\
    'El argumento debe ser una lista o una tupla'
    
    return sum(conjunto)/float(len(conjunto))

def media_geometrica(conjunto):
    '''
    parametros(conjunto) donde conjunto es una lista o tupla con todos
    los datos.
    
    Devuelve la media geometrica del conjunto de datos.
    '''
    assert type(conjunto)==list or type(conjunto)==tuple,\
    'El argumento debe ser una lista o una tupla'
    producto=1
    for i in cojunto:
        producto*=i
    return (producto)^(1./len(conjunto))
    
def varianza(x):
    '''
    parametros(conjunto) donde conjunto es una lista o tupla con todos
    los datos.
    
    Devuelve la varianza del conjunto de datos.
                _                   _
               \       _  2        \     2
           2   /_ ( xi-x )  fi  =  /_ (xi) fi      _2
      sigma = ----------------     ------------  - x
                     N                  N
    '''
    assert type(x)==list or type(x)==tuple,\
    'El argumento debe ser una lista o una tupla'
    
    media=media_aritmetica(x)
    sumatorio=0
    for i in x:
        sumatorio+=(i*i)
    return float(sumatorio)/len(x) - media*media

def desviacion_tipica(x):
    '''
    parametros(x) donde x es una lista o tupla con todos los datos.
    
    Devuelve la varianza del conjunto de datos.
                  __________________          _____________________
                 /  _                        /  _
                /  \       _  2             /  \      2
      sigma =  /   /_ ( xi-x )             /   /_ (xi)         _ 2
              /   ----------------    =   /   ------------  -  x 
            \/           N              \/          N
    '''
    return sqrt(varianza(x))

def covarianza(x,y):
    '''
    parametros (x,y) donde x,y son dos listas o tuplas de pares de datos
    en la forma x=[x1,x2,...,xn] y=[y1,y2,...,yn]
    
    Devuleve la covarianza de los pares de datos:
                 _                    _
                \      _     _       \
                /_ (xi-x)(yj-y)      /_ (xi yj)    _ _
    sigma x,y = ----------------  =  ----------- - x y
                        N                N
    '''
    assert type(x)==list or type(x)==tuple,\
    'El primer argumento debe ser una lista o una tupla'
    assert type(x)==list or type(x)==tuple,\
    'El segundo argumento debe ser una lista o una tupla'
    assert len(x)==len(y),\
    'Ambas listas deben tener el mismo numero de elementos'
    
    xmedia=media_aritmetica(x)
    ymedia=media_aritmetica(y)
    sumatorio=0 #inicializamos variable
    for i in range(len(x)): #x,y de igual dimension por construccion
        sumatorio+=(x[i]*y[i])
    return ((float(sumatorio)/len(x))-(xmedia*ymedia))

def regresion_lineal(x,y):
    '''
    parametros(x,y) donde x,y son dos listas o tuplas de pares de datos
    
    Devuelve una tupla (a,b,r2) de la recta de mejor ajuste en la forma
     y = ax + b, r2: coeficiente correlacion r2
     
                                             _         _
                                            \         \            _ 2
          covarianza       _     _    2   b /_ yi + a /_ xi yi - N y
      a = ----------- ; b= y - a x ; r = ----_--------------------------
           varianza                         \    2      _ 2
                                            /_ yi   - N y
              
    '''
    #datos valores
    xmedia=media_aritmetica(x)
    ymedia=media_aritmetica(y)
    N=len(x) #x,y igual dimension por construccion.
    
    #valores ajuste
    a=float(covarianza(x,y))/(varianza(x))
    b=ymedia-a*xmedia
    
    #calculo r2
    sumay=sum(y)
    sumay2=sum([i*i for i in y])
    sumaxy=sum([x[i]*y[i] for i in range(N)])
    Ny2=N*ymedia*ymedia
    
    r2=(b*sumay+a*sumaxy-Ny2)/(sumay2-Ny2)
    
    
    return a, b, r2

def coeficiente_correlacion_pearson(x,y):
    '''
    parametros(x,y) donde x,y son dos listas o tuplas de pares de datos
    
    Devuelve r, coeficiente de correlacion de pearson de los pares de x
    e y
    '''
    return (covarianza(x,y))/(varianza(x)*varianza(y))


