#define JPEG_INTERNALS
extern "C"{
#include <jpeglib.h>
#include <jerror.h>
#include <jconfig.h>
#include <jmorecfg.h>
#include <jdct.h>
}
#include <setjmp.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <math.h>

class Image {
public:
  int readImage(const char *filename);
  void writePGMimage(const char *filename, short *matDCT);
  unsigned int getRows();
  unsigned int getCols();
  unsigned char* getRawImg();
  unsigned int getElemQuantTbl(const unsigned int num);
  void getQuantMat(const char *filename, float *quant);
private:
  jpeg_decompress_struct cinfo;
  JSAMPROW rawImg;
};
