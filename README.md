# Parallel-k-means

#compile

gcc -Wall -fopenmp -o omp.out main_omp.c utils/image_io.c utils/segmentation_omp.c -lm

#run

.\omp.out -k 4 -t 4 imgs/large.jpg


## Documentación de Algoritmo K-means 🚀

# Metodo (Worker): 📋

En esta clase se inicializan las variables de worker, las cuales serán compartidas a todos los nodos correspondientes. 

# Metodo (squared_norm): 📋

Aqui se realiza el cálculo de la distancia euclidiana.

```
double Worker::squared_norm(Point p1, Point p2) --> recibe como entrada 2 puntos a considerar.

sum += pow(p1.values[j] - p2.values[j], 2.0); --> Calculo de la distancia, se realiza una sumatoria para todos los valores.

```

### FileManager 📋

La clase FileManager se un archivo csv que contiene todos los puntos que representa los pixeles de una imagen y se define 
la cantidad de clusters y la maxima cantidad de iteraciones (Esto se realiza solo en el nodo 0).

```
void FileManager::createPointsFile()

std::vector<u_int8_t> pixel = imgOriginal.getPixel(j,i); --> se extraen las coordenadas de cada pixel que representa a la imagen,
se hace de manera iterativa, de acuerdo a las dimensiones de la imagen.

myfile << (int)pixel[n] ; --> se guarda todos estos puntos en el documento csv para usarlos durante la clasificación k-means en los nodos.

```








