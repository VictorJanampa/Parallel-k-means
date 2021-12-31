# Parallel-k-means

#compile

gcc -Wall -fopenmp -o omp.out main_omp.c utils/image_io.c utils/segmentation_omp.c -lm

#run

.\omp.out -k 4 -t 4 imgs/large.jpg


## Documentación de Algoritmo K-means 🚀

#Metodo (Worker): 📋

En esta clase se inicializan las variables de worker, las cuales serán compartidas a todos los nodos correspondientes. 

#Metodo (squared_norm): 📋

Aqui se realiza el cálculo de la distancia euclidiana.

```
double Worker::squared_norm(Point p1, Point p2) --> recibe como entrada 2 puntos a considerar.

sum += pow(p1.values[j] - p2.values[j], 2.0); -->Calculo de la distancia, se realiza una sumatoria para todos los valores.

```

### Pre-requisitos 📋


