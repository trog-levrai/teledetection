# include <gtk/gtk.h>
# include <stdio.h>
# include <stdlib.h>
# include <string.h>
# include <time.h>

# include "compute.h"
# include "math.h"

GdkPixbuf** get_images(int* nb_images) {
  DIR           *d;
  struct dirent *dir;
  *nb_images = 0;
  d = opendir("./");

  GdkPixbuf** out = NULL;

  if (d) {

    size_t size = 0;

    // Read all images in directory
    while ((dir = readdir(d)) != NULL) {
      if (strcmp(dir->d_name, ".") && strcmp(dir->d_name, "..")) {
        size++;
        out = realloc(out, size * sizeof (GdkPixbuf*));

        char file[1000];
        strcpy(file, "./");
        strcat(file, "/");
        strcat(file, dir->d_name);

        GtkWidget *image  = gtk_image_new_from_file(file);
        if (image == NULL) {
          printf("%s, Wrong file.\n", file);
          continue;
        }

        GdkPixbuf* pGdkPixbuf = gtk_image_get_pixbuf(GTK_IMAGE(image));
        if (pGdkPixbuf == NULL) {
          printf("%s, Wrong file.\n", file);
          continue;
        }

        // Greyscale all images

        int NbCol  = gdk_pixbuf_get_width(pGdkPixbuf);
        int NbLine = gdk_pixbuf_get_height(pGdkPixbuf);

        GdkPixbuf* pGdkPixbufImaRes = gdk_pixbuf_copy(pGdkPixbuf);

        guchar* pucImaOrig = gdk_pixbuf_get_pixels(pGdkPixbuf);
        guchar* pucImaRes   = gdk_pixbuf_get_pixels(pGdkPixbufImaRes);

        ComputeImage(pucImaOrig, NbLine, NbCol, pucImaRes);

        out[size - 1] = pGdkPixbufImaRes;

        size_t count = 0;
        for (size_t i = 0; i < NbLine * NbCol; ++i) {
          if (pucImaRes[i * 3] > 0)
            ++count;
        }

        float r = (float)count / (float)(NbCol * NbLine);
        printf("%s: %f\tpercents of clouds founds\n", file, r * 100);

        *nb_images = *nb_images + 1;
      }
    }
    closedir(d);
  }

  return out;
}

int main(int argc, char* argv[])
{
  int nb = 0;
  get_images(&nb);

  return 0;
}
