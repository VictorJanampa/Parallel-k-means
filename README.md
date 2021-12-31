# Parallel-k-means

#compile

mpic++ -fopenmp -o main main.cpp Worker.cpp jpeg.cpp FileManager.cpp -ljpeg

#run

mpirun -np 4 ./main 4 10 data/test.jpg


## Documentación de Algoritmo K-means 🚀

### Worker: 📋

En esta clase se inicializan las variables de worker, las cuales serán compartidas a todos los nodos correspondientes. 

### Squared_norm: 📋

Aqui se realiza el cálculo de la distancia euclidiana.

```
double Worker::squared_norm(Point p1, Point p2) --> recibe como entrada 2 puntos a considerar.

sum += pow(p1.values[j] - p2.values[j], 2.0); --> Calculo de la distancia, se realiza una sumatoria para todos los valores.

```

### CreateDataset (FileManager) 📋

La clase FileManager se un archivo csv que contiene todos los puntos que representa los pixeles de una imagen y se define 
la cantidad de clusters y la maxima cantidad de iteraciones (Esto se realiza solo en el nodo 0).

```
void Worker::createDataset(std::string FileImg) --> se definen las variables cantidad de clusters y la maxima cantidad de iteraciones
  if (rank == 0) --> (Esto se realiza solo en el nodo 0)

```

```
void FileManager::createPointsFile()

std::vector<u_int8_t> pixel = imgOriginal.getPixel(j,i); --> se extraen las coordenadas de cada pixel que representa a la imagen,
se hace de manera iterativa, de acuerdo a las dimensiones de la imagen.

myfile << (int)pixel[n] ; --> se guarda todos estos puntos en el documento csv para usarlos durante la clasificación k-means en los nodos.

```

### ReadDataset 📋

Permite leer el archivo csv que contiene los puntos de la imagen, esto solo se realiza en el nodo 0. 

```
while (getline(infile, line, '\n')) 
if (count == 0) -->lectura del archivo en el rank 0.

else
{ Point point; point.id = num; point.size = total_values;...  --> para el resto de ranks se obtiene las coordenadas de los puntos de la imagen.

```

### ScatterDataset 📋

Realiza un Scatter de los puntos totales  (numPoints / numNodes)  a cada nodo de acuerdo a la cantidad de procesadores 
y se realiza un Bcast de la maxima cantidad de iteraciones, y la cantidad de puntos totales.

```
MPI_Scatterv(dataset.data(), pointsPerNode, datasetDisp, pointType, localDataset.data(), num_local_points,pointType, 0, comm);
--> Realiza el scater sobre todos los nodos.

MPI_Bcast(&total_values, 1, MPI_INT, 0, comm); --> envia la dimensión de los puntos a cada nodo.

```

### ExtractCluster 📋

Crea los clusters y realiza  un Bcast de los centroides y la cantidad de centroides.

```
MPI_Bcast(&K, 1, MPI_INT, 0, comm); --> Bcaster del número de cluster

MPI_Bcast(clusters.data(), K, pointType, 0, comm); --> Envia los valores de los centroides.
```

### GetIdNearestCluster 📋

Permite obtener el cluster de cada punto para esto se utilizo la distancia euclidiana hacia los centroides

```
for (int k = 1; k < K; k++) -->se obtiene la distancia euclidiana del resto de clusters
    {
        sum = 0.0;
        double dist;

        sum = squared_norm(clusters[k], p);
        dist = sqrt(sum);

        if (dist < min_dist)
        {
            min_dist = dist;
            idCluster = k;
        }
}
```

### Run 📋

Asigna el cluster de todos los puntos que tienen un nodo y hace un all reduce para calcular el nuevo  valor de los centroides.

```
#pragma omp parallel for shared(memCounter)
    for (int i = 0; i < localDataset.size(); i++) --> Se hace un for compartido para actualizar menCounter 


#pragma omp atomic update
            memCounter[new_mem]++; --> se declara esta operacion como atomica para evitar bloquear todo el vector, mutex de una posision en el Array


MPI_Allreduce(memCounter, resMemCounter, K, MPI_INT, MPI_SUM, comm); --> Recogemos el numero de puntos que pertenecen a cada cluster


#pragma omp parallel for  shared(reduceArr) 
    for (int i = 0; i < K; i++)
    {
        for (int j = 0; j < total_values; j++) --> Se hace un for compartido para actualizar los datos
        

MPI_Allreduce(reduceArr, reduceResults, K * total_values, MPI_DOUBLE, MPI_SUM,comm); --> Recogemos los elementos de todos los clusters


MPI_Allreduce(&notChanged, &globalNotChanged, 1, MPI_INT, MPI_SUM, comm); --> recogemos el valor de cambio
```

### UpdateLocalSum 📋

Permite realizar la sumatoria local  de las distancias de cada cluster para facilitar la obtencion de los centroides.

```
for (int i = 0; i < localDataset.size(); i++)
    {
        for (int j = 0; j < total_values; j++)
        {
            localSum[memberships[i]].values[j] += localDataset[i].values[j];
        }
    }
```

### ComputeGlobalMembership 📋

Une todos los puntos distribuidos, es decir, obtener el cluster de todos los puntos.

```
MPI_Reduce(&localMem, &globalMember, numPoints, MPI_INT, MPI_SUM, 0, comm);
```

### WriteClusterMembership 📋

Crea un archivo csv que contiene los datos de cada cluster y crea la imagen segmentada


