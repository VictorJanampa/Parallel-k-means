# Parallel-k-means

#compile

gcc -Wall -fopenmp -o omp.out main_omp.c utils/image_io.c utils/segmentation_omp.c -lm

#run

.\omp.out -k 4 -t 4 imgs/large.jpg




Documentación de Algoritmo K-means

Metodo Worker:  En esta clase se inicializan las variables de worker, las cuales serán compartidas a todos los nodos correspondientes. 
