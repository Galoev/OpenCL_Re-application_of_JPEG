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

int Image::readColorImage(const char *filename)
{
  /* We use our private extension JPEG error handler.
   * Note that this struct must live as long as the main JPEG parameter
   * struct, to avoid dangling-pointer problems.
   */
  struct my_error_mgr jerr;
  /* More stuff */
  FILE * infile;    /* source file */
  JSAMPARRAY buffer;    /* Output row buffer */
  int row_stride;    /* physical row width in output buffer */
  
  /* In this example we want to open the input file before doing anything else,
   * so that the setjmp() error recovery below can assume the file is open.
   * VERY IMPORTANT: use "b" option to fopen() if you are on a machine that
   * requires it in order to read binary files.
   */
  
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
  
  (void) jpeg_read_header(&cinfo, TRUE);
  /* We can ignore the return value from jpeg_read_header since
   *   (a) suspension is not possible with the stdio data source, and
   *   (b) we passed TRUE to reject a tables-only JPEG file as an error.
   * See libjpeg.txt for more info.
   */
  
  originalCols = cinfo.image_width;
  originalRows = cinfo.image_height;
  
  if (cinfo.image_width % 8 != 0) {
    cinfo.image_width += (8 - cinfo.image_width % 8);
  }
  
  if (cinfo.image_height % 8 != 0) {
    cinfo.image_height += (8 - cinfo.image_height % 8);
  }
  
  /* Step 4: set parameters for decompression */
  
  /* In this example, we don't need to change any of the defaults set by
   * jpeg_read_header(), so we do nothing here.
   */
  cinfo.out_color_space = JCS_YCbCr;
  
  /* Step 5: Start decompressor */
  
  (void) jpeg_start_decompress(&cinfo);
  /* We can ignore the return value since suspension is not possible
   * with the stdio data source.
   */
  
  /* We may need to do some setup of our own at this point before reading
   * the data.  After jpeg_start_decompress() we have the correct scaled
   * output image dimensions available, as well as the output colormap
   * if we asked for color quantization.
   * In this example, we need to make an output work buffer of the right size.
   */
  /* JSAMPLEs per row in output buffer */
  row_stride = cinfo.output_width * cinfo.output_components;
  /* Make a one-row-high sample array that will go away when done with image */
  buffer = (*cinfo.mem->alloc_sarray)
  ((j_common_ptr) &cinfo, JPOOL_IMAGE, row_stride, 1);
  rawImgY = new JSAMPLE [cinfo.image_height*cinfo.image_width];
  rawImgCb = new JSAMPLE [cinfo.image_height*cinfo.image_width];
  rawImgCr = new JSAMPLE [cinfo.image_height*cinfo.image_width];
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
    (void) jpeg_read_scanlines(&cinfo, buffer, 1);
    for (size_t i = 0; i < row_stride; i += cinfo.output_components) {
      //rawImgMy[index++] = buffer[0][i];
      rawImgY[index] = buffer[0][i];
      rawImgCb[index] = buffer[0][i+1];
      rawImgCr[index++] = buffer[0][i+2];
    }
    /* Assume put_scanline_someplace wants a pointer and sample count. */
  }
  
  /* Step 7: Finish decompression */
  
  (void) jpeg_finish_decompress(&cinfo);
  /* We can ignore the return value since suspension is not possible
   * with the stdio data source.
   */
  
  /* After finish_decompress, we can close the input file.
   * Here we postpone it until after no more JPEG errors are possible,
   * so as to simplify the setjmp error logic above.  (Actually, I don't
   * think that jpeg_destroy can do an error exit, but why assume anything...)
   */
  fclose(infile);
  
  /* At this point you may want to check to see whether any corrupt-data
   * warnings occurred (test whether jerr.pub.num_warnings is nonzero).
   */
  
  /* And we're done! */
  return 1;
}

void Image::writeJPEGimage(const char *filename)
{
  /* This struct contains the JPEG compression parameters and pointers to
   * working space (which is allocated as needed by the JPEG library).
   * It is possible to have several such structures, representing multiple
   * compression/decompression processes, in existence at once.  We refer
   * to any one struct (and its associated working data) as a "JPEG object".
   */
  struct jpeg_compress_struct comp_cinfo;
  /* This struct represents a JPEG error handler.  It is declared separately
   * because applications often want to supply a specialized error handler
   * (see the second half of this file for an example).  But here we just
   * take the easy way out and use the standard error handler, which will
   * print a message on stderr and call exit() if compression fails.
   * Note that this struct must live as long as the main JPEG parameter
   * struct, to avoid dangling-pointer problems.
   */
  struct jpeg_error_mgr jerr;
  /* More stuff */
  FILE * outfile;    /* target file */
  JSAMPROW row_pointer[1];  /* pointer to JSAMPLE row[s] */
  int row_stride;    /* physical row width in image buffer */
  
  /* Step 1: allocate and initialize JPEG compression object */
  
  /* We have to set up the error handler first, in case the initialization
   * step fails.  (Unlikely, but it could happen if you are out of memory.)
   * This routine fills in the contents of struct jerr, and returns jerr's
   * address which we place into the link field in cinfo.
   */
  comp_cinfo.err = jpeg_std_error(&jerr);
  /* Now we can initialize the JPEG compression object. */
  jpeg_create_compress(&comp_cinfo);
  
  /* Step 2: specify data destination (eg, a file) */
  /* Note: steps 2 and 3 can be done in either order. */
  
  /* Here we use the library-supplied code to send compressed data to a
   * stdio stream.  You can also write your own code to do something else.
   * VERY IMPORTANT: use "b" option to fopen() if you are on a machine that
   * requires it in order to write binary files.
   */
  if ((outfile = fopen(filename, "wb")) == NULL) {
    fprintf(stderr, "can't open %s\n", filename);
    exit(1);
  }
  jpeg_stdio_dest(&comp_cinfo, outfile);
  
  /* Step 3: set parameters for compression */
  
  /* First we supply a description of the input image.
   * Four fields of the cinfo struct must be filled in:
   */
  comp_cinfo.image_width = originalCols;   /* image width and height, in pixels */
  comp_cinfo.image_height = originalRows;
  comp_cinfo.jpeg_width = originalCols;
  comp_cinfo.jpeg_height = originalRows;
  comp_cinfo.input_components = 3;    /* # of color components per pixel */
  comp_cinfo.in_color_space = JCS_YCbCr;   /* colorspace of input image */
  /* Now use the library's routine to set default compression parameters.
   * (You must set at least cinfo.in_color_space before calling this,
   * since the defaults depend on the source color space.)
   */
  jpeg_set_defaults(&comp_cinfo);
  /* Now you can set any non-default parameters you wish to.
   * Here we just illustrate the use of quality (quantization table) scaling:
   */
  jpeg_set_quality(&comp_cinfo, 100, TRUE /* limit to baseline-JPEG values */);
  
  /* Step 4: Start compressor */
  
  /* TRUE ensures that we will write a complete interchange-JPEG file.
   * Pass TRUE unless you are very sure of what you're doing.
   */
  jpeg_start_compress(&comp_cinfo, TRUE);
  
  /* Step 5: while (scan lines remain to be written) */
  /*           jpeg_write_scanlines(...); */
  
  /* Here we use the library's state variable cinfo.next_scanline as the
   * loop counter, so that we don't have to keep track ourselves.
   * To keep things simple, we pass one scanline per call; you can pass
   * more if you wish, though.
   */
  int tmpWidth = getCols();
  row_stride = tmpWidth * 3;  /* JSAMPLEs per row in image_buffer */
  while (comp_cinfo.next_scanline < comp_cinfo.image_height) {
    /* jpeg_write_scanlines expects an array of pointers to scanlines.
     * Here the array is only one element long, but you could pass
     * more than one scanline at a time if that's more convenient.
     */
    JSAMPROW tmp = new JSAMPLE [row_stride];
    int index = comp_cinfo.next_scanline;
    for (int i = 0; i < row_stride; i += 3){
      tmp[i] = rawImgY[index*tmpWidth+(i/3)];
      tmp[i+1] = rawImgCb[index*tmpWidth+(i/3)];
      tmp[i+2] = rawImgCr[index*tmpWidth+(i/3)];
    }
    row_pointer[0] = tmp;
    (void) jpeg_write_scanlines(&comp_cinfo, row_pointer, 1);
  }
  
  /* Step 6: Finish compression */
  
  jpeg_finish_compress(&comp_cinfo);
  /* After finish_compress, we can close the output file. */
  fclose(outfile);
  
  /* Step 7: release JPEG compression object */
  
  /* This is an important step since it will release a good deal of memory. */
  jpeg_destroy_compress(&comp_cinfo);
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

unsigned char* Image::getYRawImg()
{
  return rawImgY;
}
unsigned char* Image::getCbRawImg()
{
  return rawImgCb;
}
unsigned char* Image::getCrRawImg()
{
  return rawImgCr;
}

unsigned int Image::getElemQuantTbl(const unsigned int num)
{
  return cinfo.quant_tbl_ptrs[0]->quantval[num];
}

unsigned int Image::getElemQuantTblTwo(const unsigned int num)
{
  return cinfo.quant_tbl_ptrs[0]->quantval[num];
}

unsigned int Image::getElemQuantTblThree(const unsigned int num)
{
  return cinfo.quant_tbl_ptrs[0]->quantval[num];
}

void Image::setYRawImg(unsigned char* Yimg)
{
  rawImgY = Yimg;
}
void Image::setCbRawImg(unsigned char* Cbimg)
{
  rawImgCb = Cbimg;
}
void Image::setCrRawImg(unsigned char* Crimg)
{
  rawImgCr = Crimg;
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
