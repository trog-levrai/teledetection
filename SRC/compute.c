#include <stdlib.h>
#include <stdio.h>
#include <gtk/gtk.h>
#include <math.h>

#define K 8
#define EPSILON 0.00001

/*******************************************************
IL EST FORMELLEMENT INTERDIT DE CHANGER LE PROTOTYPE
DES FONCTIONS
*******************************************************/


/*---------------------------------------
  Proto:


  But:

  Entrees:
          --->le tableau des valeurs des pixels de l'image d'origine
	      (les lignes sont mises les unes à la suite des autres)
	  --->le nombre de lignes de l'image,
	  --->le nombre de colonnes de l'image,
          --->le tableau des valeurs des pixels de l'image resultat
	      (les lignes sont mises les unes à la suite des autres)


  Sortie:

  Rem:

  Voir aussi:

  ---------------------------------------*/

//Returns the squared distance between a pixel and a mean
double get_dist(guchar* pixel, double* mean) {
  double ans = 0.;
  for (size_t i = 0; i < 5; ++i) {
    double aux = (double)pixel[i] - mean[i];
    ans += aux * aux;
  }
  return ans;
}

//Returns the nearest mean to the pixel at (i,j)
unsigned char get_mean(guchar* pixel, double** means) {
  unsigned char ans = 0;
  float min = 1024.;
  for (size_t k = 0; k < K; ++k) {
    float diff = get_dist(pixel, means[k]);
    if (min > diff) {
      min = diff;
      ans = k;
    }
  }
  return ans;
}

void update_avg(guchar* pixel, double* mean) {
  for (size_t i = 0; i < 5; ++i)
    mean[i] += pixel[i];
}

void average_mean(size_t coeff, double* mean) {
  for (size_t i = 0; i < 5; ++i)
    mean[i] /= coeff > 0 ? coeff : 1;
}

double** compute_new_means(unsigned char** mapping, guchar** image, int nbLine, int nbCol) {
  double** means = malloc(sizeof (double*) * K);
  for (size_t i = 0; i < K; ++i)
    means[i] = calloc(sizeof (double), 5);

  size_t mean_size [K];
  for (size_t i = 0; i < K; ++i)
    mean_size[i] = 0;

  for (size_t i = 0; i < nbCol; ++i) {
    for (size_t j = 0; j < nbLine; ++j) {
      update_avg(image[j * nbCol + i], means[mapping[i][j]]);
      ++mean_size[mapping[i][j]];
    }
  }

  for (size_t i = 0; i < K; ++i) {
    average_mean(mean_size[i], means[i]);
    printf("%zu ", mean_size[i]);
  }
  printf("\n");
  return means;
}

double get_delta(double** means, double** new_means) {
  double ans = 0.;
  for (size_t i = 0; i < K; ++i) {
    double diff = 0.;
    for (size_t j = 0; j < 5; ++j)
      diff += (means[i][j] - new_means[i][j]) * (means[i][j] - new_means[i][j]);
    if (diff > ans)
      ans = diff;
  }
  return ans;
}

static int sort_lambda(const void* e1, const void* e2) {
  return *(unsigned char*)e1 - *(unsigned char*)e2;
}

double* get_init_mean(unsigned char amp, unsigned char min) {
  double* mean = malloc(sizeof (double) * 5);
  for (size_t i = 0; i < 5; ++i) {
    double val = (double)min + (((double)amp / (double)K) * (double)i);
    val += (double) amp / (2. * (double)K);
    mean[0] = val;
    mean[1] = val;
    mean[2] = val;
    mean[3] = val;
    mean[4] = val;
  }
  return mean;
}

void ComputeImage(guchar *pucImaOrig,
		  int NbLine,
		  int NbCol,
		  guchar *pucImaRes)
{
  int iNbPixelsTotal, iNumPix;
  int iNumChannel, iNbChannels=3; /* on travaille sur des images couleurs*/
  guchar ucMeanPix;

  printf("Segmentation de l'image.... A vous!\n");

  //Greyscale of sorted neighbours vectors
  guchar** image = malloc(sizeof (guchar*) * NbCol * NbLine);
  guchar min = 255;
  guchar max = 0;

  for (size_t i = 0; i < NbCol * NbLine; ++i) {
    guchar val = (guchar)(pucImaOrig[i * 3] + pucImaOrig[i * 3 + 1]  + pucImaOrig[i * 3 + 2]) / 3.;
    if (val < min)
      min = val;
    if (val > max)
      max = val;
    image[i] = calloc(sizeof (guchar), 5);
    image[i][0] = (guchar)(pucImaOrig[i * 3] + pucImaOrig[i * 3 + 1]  + pucImaOrig[i * 3 + 2]) / 3.;
    if (i + 1 < NbCol * NbLine)
      image[i][1] = (guchar)(pucImaOrig[(i+1) * 3] + pucImaOrig[(i+1) * 3 + 1]  + pucImaOrig[(i+1) * 3 + 2]) / 3.;
    if (i - 1 > 0)
      image[i][2] = (guchar)(pucImaOrig[(i-1) * 3] + pucImaOrig[(i-1) * 3 + 1]  + pucImaOrig[(i-1) * 3 + 2]) / 3.;
    if (i - NbCol > 0)
      image[i][3] = (guchar)(pucImaOrig[(i - NbCol) * 3] + pucImaOrig[(i - NbCol) * 3 + 1]  + pucImaOrig[(i - NbCol) * 3 + 2]) / 3.;
    if (i + NbCol < NbCol * NbLine)
      image[i][4] = (guchar)(pucImaOrig[(i + NbCol) * 3] + pucImaOrig[(i + NbCol) * 3 + 1]  + pucImaOrig[(i + NbCol) * 3 + 2]) / 3.;
    qsort(image[i], 5, sizeof (guchar), sort_lambda);
  }

  //Initialisation des K poids de manière régulière
  double** means = malloc(sizeof (double*) * K);
  for (size_t i = 0; i < K; ++i) {
    means[i] = get_init_mean(max - min, min);
  }

  //Initialisation du mapping pixel -> classe
  unsigned char** mapping = malloc(sizeof (unsigned char*) * NbCol);
  for (size_t i = 0; i < NbCol; ++i) {
    mapping[i] = malloc(sizeof (unsigned char) * NbLine);
    for (size_t j = 0; j < NbLine; ++j)
      mapping[i][j] = get_mean(image[j * NbCol + i], means);
  }

  char changed = 1;
  int counter = 0;
  while (changed) {
    //Compute new means
    double** new_means = compute_new_means(mapping, image, NbLine, NbCol);
    if (get_delta(means, new_means) < EPSILON)
      changed = 0;
    free(means);
    means = new_means;

    //Compute new mapping
    for (size_t i = 0; i < NbCol; ++i) {
      for (size_t j = 0; j < NbLine; ++j)
        mapping[i][j] = get_mean(image[j * NbCol + i], means);
    }
    ++counter;
  }
  printf("Counter is %d\n", counter);

  iNbPixelsTotal = NbCol * NbLine;
  for (iNumPix=0;
      iNumPix<iNbPixelsTotal*iNbChannels;
      iNumPix += iNbChannels){
    //[>moyenne sur les composantes RVB <]
    ucMeanPix=(unsigned char)
        ((
          *(pucImaOrig+iNumPix) +
          *(pucImaOrig+iNumPix+1) +
          *(pucImaOrig+iNumPix+2))/3);
    ucMeanPix = mapping[(iNumPix / 3) % NbCol][(iNumPix / 3) / NbCol] == K - 1 ? 255 : 0;

    //[> sauvegarde du resultat <]
    for (iNumChannel=0;
        iNumChannel<iNbChannels;
        iNumChannel++)
      *(pucImaRes+iNumPix+iNumChannel)= ucMeanPix;
  }
/*  for (size_t i = 0; i < NbCol; ++i) {
    for (size_t j = 0; j < NbLine; ++j) {
      guchar val = mapping[i][j] == K - 1 ? 255 : 0;
      pucImaRes[i + NbCol * j] = val;
      pucImaRes[i + NbCol * j + 1] = val;
      pucImaRes[i + NbCol * j + 2] = val;
    }
  }*/

  free(means);
  free(image);
  for (size_t i = 0; i < NbCol; ++i)
    free(mapping[i]);
  free(mapping);
}
