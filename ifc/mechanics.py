from numpy import sqrt

def momento_inercia_cilindro_macizo(m,r):
    '''
    parametros (m,r)
    
    Devuleve en kg m2 el momento de inercia de un cilindro macizo de 
    masa m (kg) y radio r (m) que gira entorno a su eje longitudinal. 
    Unidades SI.
    '''
    return .5*m*r*r

def momento_inercia_esfera_maciza(m,r):
    '''
    parametros (m,r)
    
    Devuleve el momento de inercia en kg m2 de una esfera maciza de 
    masa m (kg) y radio r (m).  Unidades SI.
    '''
    return .4*m*r*r

def momento_inercia_varilla(m,l):
    '''
    parametros (m,l)
    
    Devuelve el momento de inercia en kg m2 de una varilla de masa m 
    (kg)y longitud que rota respecto a un eje perpendicular que pasa 
    por su centro. Unidades SI.
    '''
    return 1./12*m*l*l

def posicion_vertical_pelota(y0,v0,t,g=9.81):
    '''
    parametros (y0,v0,t,*g=9.81)
    
    Devuelve la posicion de un movil en mrua que parte de y0 con 
    velocidad v0 en la vertical en el instante de tiempo t. Usar 
    unidades SI. Valor por defecto de la gravedad 9.81 m/s2.
    '''
    return y0+v0*t-t*t*g*.5

def densidad_media(m,v):
    '''
    Parametros(m,v)
    
    Devuelve la densidas media de un cuerpo de masa m y volumen v, en
    unidades de los valores introduciodos.
    '''
    return float(m)/v

def mrua_tiempo(x,x0,v0,a,t0=0):
    '''
    parametros(x,x0,v0,a,t0)
    
    Devuelve el intante de tiempo en el que se encuentra un movil en 
    movimiento rectilineo uniformemente acelerado que parte en el intan-
    te t0 de x0 con velocidad v0 y aceleracion a. Todas las unidades
    en SI. Si t0, no es especificado se asume t0=0.
    '''
    return t0+((-v0-sqrt(v0*v0-2.*a*(x0-x)))/a)

def mrua_velocidad_dando_x(x,x0,v0,a,t0=0):
    '''
    parametros(x,x0,v0,a,t0)
    
    Devuelve la velocidad v que lleva un movil en movimiento rectilineo
    uniformemente acelerado que parte en el intan te t0 de x0 con velo-
    cidad v0 y aceleracion a. Todas las unidades en SI. Si t0, no es es-
    pecificado se asume t0=0.
    '''
    return v0+a*mrua_tiempo(x,x0,v0,a,t0)

def mrua_velocidad(t,a,t0=0,v0=0):
    '''
    Devuelve la velocidad de un movil en el intante t que sigue mrua con
    aceleracion a, velocidad inicial x0 y instante inicial t0. Si t0 y 
    v0 no se especifican toman el valor 0. Unidades SI.
    '''    
    return v0+a*(t-t0)

def mrua_poscion(t,a,x0=0,v0=0,t0=0):
    '''
    Devuelve la posicion de un movil en el instante t que sigue mrua con
    aceleracion a, posicion inicial x0 velocidad inicial v0 e instante
    incial t0. Si, x0,v0 o t0 no se especifican toman valor 0. Unidades
    SI
    '''
    dt=t-t0
    return x0+(v0*dt)+(.5*a*(t**2))

def km2ms(v):
    '''
    Converts velocity from km/h to m/s. It's a ufunc.
    
    PARAMETERS:
        v: number, velocity in km/h
    RETURNS:
        velicity in m/s
    '''
    return v/3.6

def ms2km(v):
    '''
    Converts velocity from m/s to km/h. It's a ufunc.
    
    PARAMETERS:
        v: number, velocity in m/s
    RETURNS:
        velicity in km/h
    '''
    return v*3.6
