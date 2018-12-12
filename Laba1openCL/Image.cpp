#include "Image.hpp"
using namespace std;
#define DEQUANTIZE(coef,quantval)  (((FAST_FLOAT) (coef)) * (quantval))
#define RANGE_CENTER  CENTERJSAMPLE

void my_jpeg_idct_float (FLOAT_MULT_TYPE * quantptr,
                         JCOEFPTR coef_block,
                         float* output_buf);

struct my_error_mgr {
  struct jpeg_error_mgr pub;  /* "public" fields */
  
  jmp_buf setjmp_buffer;  /* for return to caller */
};

typedef struct my_error_mgr * my_error_ptr;

/*
 * Here's the routine that will replace the standard error_exit method:
 */

METHODDEF(void)
my_error_exit(j_common_ptr cinfo)
{
  /* cinfo->err really points to a my_error_mgr struct, so coerce pointer */
  my_error_ptr myerr = (my_error_ptr)cinfo->err;
  
  /* Always display the message. */
  /* We could postpone this until after returning, if we chose. */
  (*cinfo->err->output_message) (cinfo);
  
  /* Return control to the setjmp point */
  longjmp(myerr->setjmp_buffer, 1);
}

int Image::readImage(const char *filename)
{
  struct my_error_mgr jerr;
  FILE *infile;
  JSAMPARRAY buffer;
  int row_stride;
  if ((infile = fopen(filename, "rb")) == NULL) {
    fprintf(stderr, "can't open %s\n", filename);
    return 0;
  }
  /* Step 1: allocate and initialize JPEG decompression object */
  /* We set up the normal JPEG error routines, then override error_exit. */
  cinfo.err = jpeg_std_error(&jerr.pub);
  jerr.pub.error_exit = my_error_exit;
  /* Establish the setjmp return context for my_error_exit to use. */
  if (setjmp(jerr.setjmp_buffer)) {
    /* If we get here, the JPEG code has signaled an error.
     * We need to clean up the JPEG object, close the input file, and return.
     */
    jpeg_destroy_decompress(&cinfo);
    fclose(infile);
    return 0;
  }
  /* Now we can initialize the JPEG decompression object. */
  jpeg_create_decompress(&cinfo);
  /* Step 2: specify data source (eg, a file) */
  jpeg_stdio_src(&cinfo, infile);
  /* Step 3: read file parameters with jpeg_read_header() */
  (void)jpeg_read_header(&cinfo, TRUE);
  /* Step 4: set parameters for decompression */
  /* In this example, we don't need to change any of the defaults set by
   * jpeg_read_header(), so we do nothing here.
   */
  cinfo.out_color_space = JCS_GRAYSCALE;
  /* Step 5: Start decompressor */
  (void)jpeg_start_decompress(&cinfo);
  /* JSAMPLEs per row in output buffer */
  row_stride = cinfo.output_width * cinfo.output_components;
  /* Make a one-row-high sample array that will go away when done with image */
  buffer = (*cinfo.mem->alloc_sarray)
  ((j_common_ptr)&cinfo, JPOOL_IMAGE, row_stride, 1);
  rawImg = new JSAMPLE [cinfo.image_height*cinfo.image_width];
  /* Step 6: while (scan lines remain to be read) */
  /*           jpeg_read_scanlines(...); */
  
  /* Here we use the library's state variable cinfo.output_scanline as the
   * loop counter, so that we don't have to keep track ourselves.
   */
  int index = 0;
  while (cinfo.output_scanline < cinfo.output_height) {
    /* jpeg_read_scanlines expects an array of pointers to scanlines.
     * Here the array is only one element long, but you could ask for
     * more than one scanline at a time if that's more convenient.
     */
    (void)jpeg_read_scanlines(&cinfo, buffer, 1);
    for (size_t i = 0; i < row_stride; i += cinfo.output_components) {
      //rawImgMy[index++] = buffer[0][i];
      rawImg[index++] = buffer[0][i];
    }
    /* Assume put_scanline_someplace wants a pointer and sample count. */
    //put_scanline_someplace(buffer[0], row_stride);
  }
  /* Step 7: Finish decompression */
  
  (void)jpeg_finish_decompress(&cinfo);
  /* We can ignore the return value since suspension is not possible
   * with the stdio data source.
   */
  fclose(infile);
  return 0;
}

unsigned int Image::getRows()
{
  return  cinfo.image_height;
}
unsigned int Image::getCols()
{
  return cinfo.image_width;
}
unsigned char* Image::getRawImg()
{
  return rawImg;
}

unsigned int Image::getElemQuantTbl(const unsigned int num)
{
  return cinfo.quant_tbl_ptrs[0]->quantval[num];
}

void Image::writePGMimage(const char *filename, unsigned char *matDCT){
  FILE *file = fopen(filename, "wb");
  //cout<<"Width: "<<cinfo.image_width<<" Heigh: "<<cinfo.image_height<<endl;
  fprintf(file, "P5\n%i %i\n255\n", cinfo.image_width, cinfo.image_height);
  fwrite(matDCT, sizeof(unsigned char), getRows() * getCols(), file);
  /*
  for (int i = 0; i<cinfo.image_height; i++){
    for (int j = 0; j<cinfo.image_width; j++){
      fwrite(&matDCT[i*cinfo.image_width+j], 1, 1, file);
    }
    
  }*/
  fclose(file);
}

void clearZigzag(const char* zigzahDQT, unsigned int* DQT, const int row, const int col){
  int curY = -1;
  int curX = 2;
  int y = 0;
  int x = 0;
  DQT[0] = static_cast<unsigned char>(zigzahDQT[0]);
  for (int i = 1; i<col; i++){
    for (int j = 0; j<(i+1); j++){
      if (x<(col-1)){
        x++;
      } else {
        y++;
        x = 0;
      }
      if((i+1)%2 == 0){
        curY++;
        curX--;
      } else {
        curY--;
        curX++;
      }
      DQT[y*col+x] = static_cast<unsigned char>(zigzahDQT[curY*col+curX]);
    }
    if((i+1)%2 == 0){
      curY+=2;
      curX-=1;
    } else {
      curY-=1;
      curX+=2;
    }
  }
  curY = col;
  curX = 0;
  for (int i = (col-2); i>0; i--){
    for (int j = 0; j<=i; j++){
      if (x<(col-1)){
        x++;
      } else {
        y++;
        x = 0;
      }
      if((i+1)%2 == 0){
        curY++;
        curX--;
      } else {
        curY--;
        curX++;
      }
      DQT[y*col+x] = static_cast<unsigned char>(zigzahDQT[curY*col+curX]);
    }
    if((i+1)%2 == 0){
      curY = col;
    } else {
      curX = col;
    }
  }
  DQT[col*row-1] = static_cast<unsigned char>(zigzahDQT[col*row-1]);
}

void Image::getQuantMat(const char *filename, float *quant)
{
  ifstream fin (filename, ios::binary);
  char buff;
  char check[2];
  char jpgStart[2] = {'\xff', '\xd8'};
  bool flag = true;
  fin.read(check, 2);
  for (int i = 0; i<2; i++) {
    cout<<check[i]<<endl;
    if (check[i] != jpgStart [i]) {
      flag = false;
      break;
    }
  }
  if (!flag) {
    cout<<"Это не JPG"<<endl;
    return;
  }
  
  char lenght[2];
  int x;
  fin.read(check, 2);
  while (check[0] != '\xee' && check[1] != '\xdb'){
    fin.read(lenght, 2);
    x = (lenght[0]*256) + lenght[1] - 2;
    fin.seekg(x, fin.cur);
    fin.read(check, 2);
    if (fin.eof()){
      return;
    }
  }
  
  fin.read(lenght, 2);
  fin>>buff;
  x = (lenght[0]*256) + lenght[1] - 3;
  char *DQT = new char [64];
  unsigned int *intDQT = new unsigned int[64];
  fin.read(DQT, 64);
  for (int i = 0; i<64; i++)
    intDQT[i] = 0;
  fin.close();
  clearZigzag(DQT, intDQT, 8, 8);
  for (int i = 0; i<DCTSIZE2; i++){
    quant[i] = intDQT[i];
  }
}
