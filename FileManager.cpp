#include "FileManager.h"
#include <fstream>
#include "jpeg.h"

FileManager::FileManager( int numClusters, int maxIteration, std::string filename,std::string fileimg) {
    marengo::jpeg::Image imgOriginal(fileimg);
    this->numPoints = imgOriginal.getWidth() * imgOriginal.getHeight();
    this->pointDimension = imgOriginal.getPixelSize();
    this->numClusters = numClusters;
    this->maxIteration = maxIteration;
    this->filename = filename;
    this->fileimg = fileimg;
}

void FileManager::createPointsFile(){
    srand(time(NULL));
    std::ofstream myfile;
    marengo::jpeg::Image imgOriginal(fileimg);
    myfile.open("data/" + filename + ".csv");
    myfile << pointDimension << "," << numClusters << "," << maxIteration << "\n";
    for(int i = 0; i < imgOriginal.getHeight(); i++){
        for(int j=0 ;j < imgOriginal.getWidth(); j++){
            std::vector<u_int8_t> pixel = imgOriginal.getPixel(j,i); /*se extraen las coordenadas de cada pixel que representa a la imagen*/
            for(int n = 0;n < pointDimension;n++){ /*se hace de manera iterativa, de acuerdo a las dimensiones de la imagen*/
                myfile << (int)pixel[n] ; /*se guarda todos estos puntos en el documento csv para usarlos durante la clasificación k-means en los nodos*/
                if(n < pointDimension-1){
                    myfile << ",";
                }
            }
            myfile << "\n";

        }

    }
    myfile.close();
}
