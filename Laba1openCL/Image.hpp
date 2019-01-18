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
  int readColorImage(const char *filename);
  void writePGMimage(const char *filename, unsigned char *matDCT);
  void writeJPEGimage(const char *filename);
  unsigned int getRows();
  unsigned int getCols();
  unsigned char* getRawImg();
  unsigned char* getYRawImg();
  unsigned char* getCbRawImg();
  unsigned char* getCrRawImg();
  void setYRawImg(unsigned char* Yimg);
  void setCbRawImg(unsigned char* Cbimg);
  void setCrRawImg(unsigned char* Crimg);
  unsigned int getElemQuantTbl(const unsigned int num);
  unsigned int getElemQuantTblTwo(const unsigned int num);
  unsigned int getElemQuantTblThree(const unsigned int num);
  void getQuantMat(const char *filename, float *quant);
private:
  jpeg_decompress_struct cinfo;
  JSAMPROW rawImg;
  JSAMPROW rawImgY;
  JSAMPROW rawImgCb;
  JSAMPROW rawImgCr;
  int mySize = 0;
  int originalRows = 0;
  int originalCols = 0;
};
