#Recopilation module of funtions used on IFC
#
#encoding: ASCII
#documentation language: spanish
#Author: Sergio Garcia Pajares
#Last update: 2019-03-01  21:55

'''
Recolpilacion de funciones empleadas en IFC durante el periodo 
2018-2019.

Ultima actualizacion: 01-03-2019
'''
#== MODULO PRINCIPAL ===================================================

#DATA:
pi= 3.141592653589793
e = 2.718281828459045

#IMPORTS



#FUNCTIONS:
def horizontal_line():
    '''
    Imprime por pantalla una linea horizontal de ancho 81.
    '''
    print('------------------------------------------------------------\
--------------------')
   
def read_csv_rows(file_name,delimitator=',',comment='#',empty=None):
    '''
    Returns a list of list of rows of date (making them an eval) with
    shape:
    
    
      CSV:                     RETURN:
     a1,b1,...,z1    [ [a1,b1,...,z1],[a2,b2,...,z2],...[an,bn,...,zn] ]
     a2,b2,...,z2
      |  |  |   |
     an,bn,...,zn
     
     If there is an empty value, list will contain a None value.
     It check's the file existence
     
     
     PARAMETERS:
        file_name, str: name or path to the csv file
        
        OPTIONAL
        delimitator, str: delimitation character
        comment, str: comment marker. It must be at the begining of the
                      line.
        empty, str: if there is an empty value this parameters specifies
                    which value will fit it's position in the list. 
    
    
    '''
    #-- comprobaciones a los parametros ------
    assert type(file_name)==str, 'el nombre debe ser una cadena de carateres'
    assert type(delimitator)==str, 'el delimitador debe ser una cadena de caracteres'
    
    #-- lectura de lineas -------
    try:
        archivo=open(file_name,'r')
    except:
        print('El fichero <%s> no existe'%(file_name))
        exit(1)
    
    raw=archivo.readlines()
    archivo.close()
    
    L=[] #Declare result list
    #-- splitting rows -----
    for i in range(len(raw)):
        if raw[i][0]!=comment: #ignore comment lines
            fila=raw[i].split(delimitator) #separamos columnas
            for j in range(len(fila)): 
                try:#writting none in case of not been able to eval
                    fila[j]=eval(fila[j]) #make eval
                except:
                    fila[j]=empty
            L.append(fila)
    
    del raw #borramos lista de proceso interno.
    return L #devolvemos la lista de listas de filas

    
def read_csv_columnas_primitivo(file_name,delimitator=','):
    '''
    Devuelve una lista de listas de columnas de datos que tuviese el csv
    en la forma:
    
      CSV:                     FUNCION:
     a1,b1,...,z1    [ [a1,a2,...,an],[b1,b2,...,bn],...[z1,z2,...,zn] ]
     a2,b2,...,z2
      |  |  |   |
     an,bn,...,zn
     
    Internamente, primero realiza lectura por filas y despues reorganiza
    a columas. La funcion ignorara todas las lineas que empiecen por '#'.
    '''
    #-- usamos lectura por filas en primera instacia ----
    filas=read_csv_filas(file_name,delimitator)
    #-- separamos en columnas ----
    
    z=range(len(filas[0])) #z es la lista de indices de elementos en
                           #cada fila, i.e. numero de columnas.
    
    columnas=[] 
    for i in z:              #creamos una lista de listas vacias con
        columnas.append([])  #el numero de elementos en cada fila, 
                             #i.e. numero de comunas.
    
    for i in range(len(filas)):
        for j in z:
            columnas[j].append(filas[i][j])
    del z, filas #borramos listas del proceso intermedio.
    
    return columnas
    
#def read_csv_columnas(file_name,delimitator=','):
#    '''
#    Devuelve una lista de listas de columnas de datos que tuviese el csv
#    en la forma:
#    
#      CSV:                     FUNCION:
#     a1,b1,...,z1    [ [a1,a2,...,an],[b1,b2,...,bn],...[z1,z2,...,zn] ]
#     a2,b2,...,z2
#      |  |  |   |
#     an,bn,...,zn
#     
#    Internamente, primero realiza lectura por filas y despues reorganiza
#    a columas. La funcion ignorara todas las lineas que empiecen por '#'.
#    '''
#    
#    archivo=open(file_name,'r')
#    linea=archivo.readline()
#    linea=linea.split(delimitador) #separamos en filas
#    columnas=len(filas[0])*[[],]
#    while linea != '': #Leer todas la lineas
#        linea=linea.split(delimitador) #separamos en filas
#        if linea[0]==('#' or '\n'): continue
#        
#        
#    #inacabado
    


def existencia_fichero(file_name):
    '''
    Devuelve True si existe y False si no lo hace, el fichero cuyo path 
    (o nombre) pasamos como argumento.
    '''
    assert type(file_name)==str,\
    'El path debe ser una cadena de caracteres'
    
    try:
        f=open(file_name,'r')
        f.close()
    except:
        return False
    else:
        return True


def overwrite(file_name):
    '''
    Checks if a file already existis. In case it does, ask the user if
    they want to replace it or not. If they want, then let's the program
    runs. Otherwise, exits ends the program
    
    PARAMETERS:
        file_name, str: name or path to the file
    RETURNS:
        nothing
        
        However, prints the replacing or not question
        
    '''
    if existencia_fichero(file_name):
        while True:
            q=raw_input("\n%s already exist. Do you want to replace it? (y/n) "%(file_name))
            
            if q == 'y' or q == 'Y':
                break #stop loop and keep executing
            elif q == 'n' or q == 'N':
                exit(0) #we dont want to overwrite, stop program
            else:
                pass #not proper answer, repeat loop
