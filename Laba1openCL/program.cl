#define GETJSAMPLE(value)  ((int) (value))
#define DEQUANTIZE(coef,quantval)  (((float) (coef)) * (quantval))
#define DCTSIZE 8
#define DCTSIZE2 64
#define CENTERJSAMPLE  128
#define LOCAL_WINDOW_SIZE 16
#define MAXJSAMPLE  255
#define RANGE_MASK  (MAXJSAMPLE * 4 + 3)



kernel void sum(const float global *a, const float global *b, float global *c)
{
  uint x = get_global_id(0); // Номер текущего потока
  c[x] = a[x] + b[x];
}





void my_stride_jpeg_fdct_float (float * data, local float* sample_data, const unsigned int stride)
{
  float tmp0, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7;
  float tmp10, tmp11, tmp12, tmp13;
  float z1, z2, z3, z4, z5, z11, z13;
  float *dataptr;
  local float* elemptr;
  int ctr;
  
  /* Pass 1: process rows. */
  
  dataptr = data;
  for (ctr = 0; ctr < DCTSIZE; ctr++) {
    //elemptr = sample_data[ctr] + start_col;
    elemptr = &sample_data[ctr*stride];
    
    /* Load data into workspace */
    tmp0 = (float) (GETJSAMPLE(elemptr[0]) + GETJSAMPLE(elemptr[7]));
    tmp7 = (float) (GETJSAMPLE(elemptr[0]) - GETJSAMPLE(elemptr[7]));
    tmp1 = (float) (GETJSAMPLE(elemptr[1]) + GETJSAMPLE(elemptr[6]));
    tmp6 = (float) (GETJSAMPLE(elemptr[1]) - GETJSAMPLE(elemptr[6]));
    tmp2 = (float) (GETJSAMPLE(elemptr[2]) + GETJSAMPLE(elemptr[5]));
    tmp5 = (float) (GETJSAMPLE(elemptr[2]) - GETJSAMPLE(elemptr[5]));
    tmp3 = (float) (GETJSAMPLE(elemptr[3]) + GETJSAMPLE(elemptr[4]));
    tmp4 = (float) (GETJSAMPLE(elemptr[3]) - GETJSAMPLE(elemptr[4]));
    
    /* Even part */
    
    tmp10 = tmp0 + tmp3;  /* phase 2 */
    tmp13 = tmp0 - tmp3;
    tmp11 = tmp1 + tmp2;
    tmp12 = tmp1 - tmp2;
    
    /* Apply unsigned->signed conversion */
    dataptr[0] = tmp10 + tmp11 - 8 * CENTERJSAMPLE; /* phase 3 */
    dataptr[4] = tmp10 - tmp11;
    
    z1 = (tmp12 + tmp13) * ((float) 0.707106781); /* c4 */
    dataptr[2] = tmp13 + z1;  /* phase 5 */
    dataptr[6] = tmp13 - z1;
    
    /* Odd part */
    
    tmp10 = tmp4 + tmp5;  /* phase 2 */
    tmp11 = tmp5 + tmp6;
    tmp12 = tmp6 + tmp7;
    
    /* The rotator is modified from fig 4-8 to avoid extra negations. */
    z5 = (tmp10 - tmp12) * ((float) 0.382683433); /* c6 */
    z2 = ((float) 0.541196100) * tmp10 + z5; /* c2-c6 */
    z4 = ((float) 1.306562965) * tmp12 + z5; /* c2+c6 */
    z3 = tmp11 * ((float) 0.707106781); /* c4 */
    
    z11 = tmp7 + z3;    /* phase 5 */
    z13 = tmp7 - z3;
    
    dataptr[5] = z13 + z2;  /* phase 6 */
    dataptr[3] = z13 - z2;
    dataptr[1] = z11 + z4;
    dataptr[7] = z11 - z4;
    
    dataptr += DCTSIZE;    /* advance pointer to next row */
  }
  
  /* Pass 2: process columns. */
  
  dataptr = data;
  for (ctr = DCTSIZE-1; ctr >= 0; ctr--) {
    tmp0 = dataptr[DCTSIZE*0] + dataptr[DCTSIZE*7];
    tmp7 = dataptr[DCTSIZE*0] - dataptr[DCTSIZE*7];
    tmp1 = dataptr[DCTSIZE*1] + dataptr[DCTSIZE*6];
    tmp6 = dataptr[DCTSIZE*1] - dataptr[DCTSIZE*6];
    tmp2 = dataptr[DCTSIZE*2] + dataptr[DCTSIZE*5];
    tmp5 = dataptr[DCTSIZE*2] - dataptr[DCTSIZE*5];
    tmp3 = dataptr[DCTSIZE*3] + dataptr[DCTSIZE*4];
    tmp4 = dataptr[DCTSIZE*3] - dataptr[DCTSIZE*4];
    
    /* Even part */
    
    tmp10 = tmp0 + tmp3;  /* phase 2 */
    tmp13 = tmp0 - tmp3;
    tmp11 = tmp1 + tmp2;
    tmp12 = tmp1 - tmp2;
    
    dataptr[DCTSIZE*0] = tmp10 + tmp11; /* phase 3 */
    dataptr[DCTSIZE*4] = tmp10 - tmp11;
    
    z1 = (tmp12 + tmp13) * ((float) 0.707106781); /* c4 */
    dataptr[DCTSIZE*2] = tmp13 + z1; /* phase 5 */
    dataptr[DCTSIZE*6] = tmp13 - z1;
    
    /* Odd part */
    
    tmp10 = tmp4 + tmp5;  /* phase 2 */
    tmp11 = tmp5 + tmp6;
    tmp12 = tmp6 + tmp7;
    
    /* The rotator is modified from fig 4-8 to avoid extra negations. */
    z5 = (tmp10 - tmp12) * ((float) 0.382683433); /* c6 */
    z2 = ((float) 0.541196100) * tmp10 + z5; /* c2-c6 */
    z4 = ((float) 1.306562965) * tmp12 + z5; /* c2+c6 */
    z3 = tmp11 * ((float) 0.707106781); /* c4 */
    
    z11 = tmp7 + z3;    /* phase 5 */
    z13 = tmp7 - z3;
    
    dataptr[DCTSIZE*5] = z13 + z2; /* phase 6 */
    dataptr[DCTSIZE*3] = z13 - z2;
    dataptr[DCTSIZE*1] = z11 + z4;
    dataptr[DCTSIZE*7] = z11 - z4;
    
    dataptr++;      /* advance pointer to next column */
  }
}


void my_jpeg_idct_float (unsigned char* range_limit, float *quantptr, short *coef_block, float *output_buf)
{
  float tmp0, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7;
  float tmp10, tmp11, tmp12, tmp13;
  float z5, z10, z11, z12, z13;
  short *inptr;
  //FLOAT_MULT_TYPE * quantptr;
  float * wsptr;
  float* outptr;
  //JSAMPLE *range_limit = cinfo->sample_range_limit;
  int ctr;
  float workspace[DCTSIZE2]; /* buffers data between passes */
  
  /* Pass 1: process columns from input, store into work array. */
  
  inptr = coef_block;
  //quantptr = (FLOAT_MULT_TYPE *) compptr->dct_table;
  wsptr = workspace;
  for (ctr = DCTSIZE; ctr > 0; ctr--) {
    /* Due to quantization, we will usually find that many of the input
     * coefficients are zero, especially the AC terms.  We can exploit this
     * by short-circuiting the IDCT calculation for any column in which all
     * the AC terms are zero.  In that case each output is equal to the
     * DC coefficient (with scale factor as needed).
     * With typical images and quantization tables, half or more of the
     * column DCT calculations can be simplified this way.
     */
    
    if (inptr[DCTSIZE*1] == 0 && inptr[DCTSIZE*2] == 0 &&
        inptr[DCTSIZE*3] == 0 && inptr[DCTSIZE*4] == 0 &&
        inptr[DCTSIZE*5] == 0 && inptr[DCTSIZE*6] == 0 &&
        inptr[DCTSIZE*7] == 0) {
      /* AC terms all zero */
      float dcval = DEQUANTIZE(inptr[DCTSIZE*0], quantptr[DCTSIZE*0]);
      
      wsptr[DCTSIZE*0] = dcval;
      wsptr[DCTSIZE*1] = dcval;
      wsptr[DCTSIZE*2] = dcval;
      wsptr[DCTSIZE*3] = dcval;
      wsptr[DCTSIZE*4] = dcval;
      wsptr[DCTSIZE*5] = dcval;
      wsptr[DCTSIZE*6] = dcval;
      wsptr[DCTSIZE*7] = dcval;
      
      inptr++;      /* advance pointers to next column */
      quantptr++;
      wsptr++;
      continue;
    }
    
    /* Even part */
    
    tmp0 = DEQUANTIZE(inptr[DCTSIZE*0], quantptr[DCTSIZE*0]);
    tmp1 = DEQUANTIZE(inptr[DCTSIZE*2], quantptr[DCTSIZE*2]);
    tmp2 = DEQUANTIZE(inptr[DCTSIZE*4], quantptr[DCTSIZE*4]);
    tmp3 = DEQUANTIZE(inptr[DCTSIZE*6], quantptr[DCTSIZE*6]);
    
    tmp10 = tmp0 + tmp2;  /* phase 3 */
    tmp11 = tmp0 - tmp2;
    
    tmp13 = tmp1 + tmp3;  /* phases 5-3 */
    tmp12 = (tmp1 - tmp3) * ((float) 1.414213562) - tmp13; /* 2*c4 */
    
    tmp0 = tmp10 + tmp13;  /* phase 2 */
    tmp3 = tmp10 - tmp13;
    tmp1 = tmp11 + tmp12;
    tmp2 = tmp11 - tmp12;
    
    /* Odd part */
    
    tmp4 = DEQUANTIZE(inptr[DCTSIZE*1], quantptr[DCTSIZE*1]);
    tmp5 = DEQUANTIZE(inptr[DCTSIZE*3], quantptr[DCTSIZE*3]);
    tmp6 = DEQUANTIZE(inptr[DCTSIZE*5], quantptr[DCTSIZE*5]);
    tmp7 = DEQUANTIZE(inptr[DCTSIZE*7], quantptr[DCTSIZE*7]);
    
    z13 = tmp6 + tmp5;    /* phase 6 */
    z10 = tmp6 - tmp5;
    z11 = tmp4 + tmp7;
    z12 = tmp4 - tmp7;
    
    tmp7 = z11 + z13;    /* phase 5 */
    tmp11 = (z11 - z13) * ((float) 1.414213562); /* 2*c4 */
    
    z5 = (z10 + z12) * ((float) 1.847759065); /* 2*c2 */
    tmp10 = z5 - z12 * ((float) 1.082392200); /* 2*(c2-c6) */
    tmp12 = z5 - z10 * ((float) 2.613125930); /* 2*(c2+c6) */
    
    tmp6 = tmp12 - tmp7;  /* phase 2 */
    tmp5 = tmp11 - tmp6;
    tmp4 = tmp10 - tmp5;
    
    wsptr[DCTSIZE*0] = tmp0 + tmp7;
    wsptr[DCTSIZE*7] = tmp0 - tmp7;
    wsptr[DCTSIZE*1] = tmp1 + tmp6;
    wsptr[DCTSIZE*6] = tmp1 - tmp6;
    wsptr[DCTSIZE*2] = tmp2 + tmp5;
    wsptr[DCTSIZE*5] = tmp2 - tmp5;
    wsptr[DCTSIZE*3] = tmp3 + tmp4;
    wsptr[DCTSIZE*4] = tmp3 - tmp4;
    
    inptr++;      /* advance pointers to next column */
    quantptr++;
    wsptr++;
  }
  
  /* Pass 2: process rows from work array, store into output array. */
  
  wsptr = workspace;
  for (ctr = 0; ctr < DCTSIZE; ctr++) {
    outptr = &output_buf[ctr*DCTSIZE];
    /* Rows of zeroes can be exploited in the same way as we did with columns.
     * However, the column calculation has created many nonzero AC terms, so
     * the simplification applies less often (typically 5% to 10% of the time).
     * And testing floats for zero is relatively expensive, so we don't bother.
     */
    
    /* Even part */
    
    /* Apply signed->unsigned and prepare float->int conversion */
    z5 = wsptr[0] + ((float) CENTERJSAMPLE + (float) 0.5);
    tmp10 = z5 + wsptr[4];
    tmp11 = z5 - wsptr[4];
    
    tmp13 = wsptr[2] + wsptr[6];
    tmp12 = (wsptr[2] - wsptr[6]) * ((float) 1.414213562) - tmp13;
    
    tmp0 = tmp10 + tmp13;
    tmp3 = tmp10 - tmp13;
    tmp1 = tmp11 + tmp12;
    tmp2 = tmp11 - tmp12;
    
    /* Odd part */
    
    z13 = wsptr[5] + wsptr[3];
    z10 = wsptr[5] - wsptr[3];
    z11 = wsptr[1] + wsptr[7];
    z12 = wsptr[1] - wsptr[7];
    
    tmp7 = z11 + z13;
    tmp11 = (z11 - z13) * ((float) 1.414213562);
    
    z5 = (z10 + z12) * ((float) 1.847759065); /* 2*c2 */
    tmp10 = z5 - z12 * ((float) 1.082392200); /* 2*(c2-c6) */
    tmp12 = z5 - z10 * ((float) 2.613125930); /* 2*(c2+c6) */
    
    tmp6 = tmp12 - tmp7;
    tmp5 = tmp11 - tmp6;
    tmp4 = tmp10 - tmp5;
    
    /* Final output stage: float->int conversion and range-limit */
    
    outptr[0] = range_limit[((int) (tmp0 + tmp7)) & RANGE_MASK];
    outptr[7] = range_limit[((int) (tmp0 - tmp7)) & RANGE_MASK];
    outptr[1] = range_limit[((int) (tmp1 + tmp6)) & RANGE_MASK];
    outptr[6] = range_limit[((int) (tmp1 - tmp6)) & RANGE_MASK];
    outptr[2] = range_limit[((int) (tmp2 + tmp5)) & RANGE_MASK];
    outptr[5] = range_limit[((int) (tmp2 - tmp5)) & RANGE_MASK];
    outptr[3] = range_limit[((int) (tmp3 + tmp4)) & RANGE_MASK];
    outptr[4] = range_limit[((int) (tmp3 - tmp4)) & RANGE_MASK];
    
    wsptr += DCTSIZE;    /* advance pointer to next row */
  }
}

void strideFloatDCT(local float *sample_data, const unsigned int stride, float *data,  float constant *scaleMatQuant,  float constant *scaleMatDequant)
{
  
  //jpeg_fdct_float(data, sample_data, 0);
  my_stride_jpeg_fdct_float(data, sample_data, stride);
  //static const double aanscalefactor[DCTSIZE] = {    1.0, 1.387039845, 1.306562965, 1.175875602,    1.0, 0.785694958, 0.541196100, 0.275899379  };
  short shortData[DCTSIZE2];
  for (int i = 0; i < 8; ++i) {
    for (int j = 0; j < 8; ++j) {
      //shortData[i * 8 + j] = data[i * 8 + j] *(1/ (8 * aanscalefactor[i] * aanscalefactor[j] * cinfo.quant_tbl_ptrs[0]->quantval[i*DCTSIZE+j]/* * quantMatrix */));
      //shortData[i * 8 + j] = data[i * 8 + j] *(1/ (8 * aanscalefactor[i] * aanscalefactor[j] ));
      //shortData[i * 8 + j] = data[i * 8 + j] *scaleMatQuant[i * 8 + j]* cinfo.quant_tbl_ptrs[0]->quantval[i*DCTSIZE+j];
      shortData[i * 8 + j] = round(data[i * 8 + j] * scaleMatQuant[i * 8 + j]);
    }
  }
  float dctTable[DCTSIZE2];
  for (int i = 0; i<DCTSIZE; i++){
    for (int j = 0; j<DCTSIZE; j++){
      //dctTable[i*DCTSIZE+j] = aanscalefactor[i]*aanscalefactor[j]*0.125 * cinfo.quant_tbl_ptrs[0]->quantval[i*DCTSIZE+j];
      //dctTable[i*DCTSIZE+j] = aanscalefactor[i]*aanscalefactor[j]*0.125;
      //dctTable[i*DCTSIZE+j] = scaleMatDequant[i*DCTSIZE+j]/cinfo.quant_tbl_ptrs[0]->quantval[i*DCTSIZE+j];
      dctTable[i*DCTSIZE+j] = scaleMatDequant[i*DCTSIZE+j];
    }
  }
  unsigned char range_limit[1024];
  for (int i = 512; i<768; i++){
    range_limit[i] = i-512;
  }
  for (int i = 0; i<512; i++){
    range_limit[i] = 0;
  }
  for (int i = 768; i<1024; i++){
    range_limit[i] = 255;
  }
  my_jpeg_idct_float(&range_limit[512], dctTable, shortData, data);
}


kernel void reapplicationJPEG(const uchar global *rawImg, const uint rows, const uint cols, short global *resDCT,  float constant *scaleMatQuant,  float constant *dctTable, const uint offset, const uint stride)
{
  
  local union{
    float mat_1D[LOCAL_WINDOW_SIZE*LOCAL_WINDOW_SIZE];
    float mat1_2D[LOCAL_WINDOW_SIZE][LOCAL_WINDOW_SIZE];
    float4 mat4_2D[LOCAL_WINDOW_SIZE][LOCAL_WINDOW_SIZE/4];
  } mat;
  
  int x2 = get_local_id(0)%4;
  int y2 = get_local_id(0)/4;
  
  mat.mat4_2D[y2][x2] = convert_float4(*(global uchar4*)&rawImg[(get_group_id(1)*LOCAL_WINDOW_SIZE+y2)*stride+get_group_id(0)*LOCAL_WINDOW_SIZE+4*x2+offset]);
  //mat.mat4_2D[y2][x2] = *(global float4*)&rawImg[(get_group_id(1)*LOCAL_WINDOW_SIZE+y2)*stride+get_group_id(0)*LOCAL_WINDOW_SIZE+4*x2+offset];
  
  barrier(CLK_LOCAL_MEM_FENCE);
  
  int x = get_local_id(0)%8;
  int y = get_local_id(0)/8;
  
  float res[DCTSIZE2];
  strideFloatDCT(&mat.mat_1D[y*LOCAL_WINDOW_SIZE+x], LOCAL_WINDOW_SIZE, res, scaleMatQuant, dctTable);
  
  barrier(CLK_LOCAL_MEM_FENCE);
  mat.mat4_2D[y2][x2] = 0.0f;
  
  barrier(CLK_LOCAL_MEM_FENCE);
  for (int i=0; i<DCTSIZE; i++)
  {
    for (int j=0; j<DCTSIZE; j++)
    {
      mat.mat1_2D[y+i][x+j] += res[i*DCTSIZE+j];
      barrier(CLK_LOCAL_MEM_FENCE);
    }
  }
  
  *(global short4*)&resDCT[(get_group_id(1)*LOCAL_WINDOW_SIZE+y2)*stride+get_group_id(0)*LOCAL_WINDOW_SIZE+4*x2+offset] += convert_short4_sat_rte(mat.mat4_2D[y2][x2]);
  
}

kernel void justFloatDCT(const uchar global *rawImg, const uint rows, const uint cols, short global *resDCT, float constant *scaleMatQuant, float constant *dctTable, const uint offset, const uint stride)
{
  local uchar mat[DCTSIZE2];
  for (int i = 0; i<DCTSIZE2; i++)
    mat[i] = rawImg[i];
  float res[DCTSIZE2];
  strideFloatDCT(mat, DCTSIZE, res, scaleMatQuant, dctTable);
  for (int i = 0; i<DCTSIZE; i++)
    for (int j = 0; j<DCTSIZE; j++)
      resDCT[i*DCTSIZE+j] = res[i*DCTSIZE+j];
}


kernel void divMat(short global *resDCT, uchar global *resImg, const uint rows, const uint cols)
{
  int x = get_global_id(0);
  int y = get_global_id(1);
  
  
  short8 MatUL[DCTSIZE] = {
    (short8)(1, 2, 3, 4, 5, 6, 7, 8),
    (short8)(2, 4, 6, 8, 10, 12, 14, 16),
    (short8)(3, 6, 9, 12, 15, 18, 21, 24),
    (short8)(4, 8, 12, 16, 20, 24, 28, 32),
    (short8)(5, 10, 15, 20, 25, 30, 35, 40),
    (short8)(6, 12, 18, 24, 30, 36, 42, 48),
    (short8)(7, 14, 21, 28, 35, 42, 49, 56),
    (short8)(8, 16, 24, 32, 40, 48, 56, 64)
  };
  
  short8 MatDL[DCTSIZE] = {
    (short8)(8, 16, 24, 32, 40, 48, 56, 64),
    (short8)(7, 14, 21, 28, 35, 42, 49, 56),
    (short8)(6, 12, 18, 24, 30, 36, 42, 48),
    (short8)(5, 10, 15, 20, 25, 30, 35, 40),
    (short8)(4, 8, 12, 16, 20, 24, 28, 32),
    (short8)(3, 6, 9, 12, 15, 18, 21, 24),
    (short8)(2, 4, 6, 8, 10, 12, 14, 16),
    (short8)(1, 2, 3, 4, 5, 6, 7, 8)
  };
  
  short8 MatUR[DCTSIZE] = {
    (short8)(8, 7, 6, 5, 4, 3, 2, 1),
    (short8)(16, 14, 12, 10, 8, 6, 4, 2),
    (short8)(24, 21, 18, 15, 12, 9, 6, 3),
    (short8)(32, 28, 24, 20, 16, 12, 8, 4),
    (short8)(40, 35, 30, 25, 20, 15, 10, 5),
    (short8)(48, 42, 36, 30, 24, 18, 12, 6),
    (short8)(56, 49, 42, 35, 28, 21, 14, 7),
    (short8)(64, 56, 48, 40, 32, 24, 16, 8)
  };
  
  short8 MatDR[DCTSIZE] = {
    (short8)(64, 56, 48, 40, 32, 24, 16, 8),
    (short8)(56, 49, 42, 35, 28, 21, 14, 7),
    (short8)(48, 42, 36, 30, 24, 18, 12, 6),
    (short8)(40, 35, 30, 25, 20, 15, 10, 5),
    (short8)(32, 28, 24, 20, 16, 12, 8, 4),
    (short8)(24, 21, 18, 15, 12, 9, 6, 3),
    (short8)(16, 14, 12, 10, 8, 6, 4, 2),
    (short8)(8, 7, 6, 5, 4, 3, 2, 1)
  };
  
  short8 tmp = *(global short8*)&resDCT[y*cols+x*DCTSIZE];
  
  if (y>=DCTSIZE && (x*DCTSIZE) >= DCTSIZE && y<(rows-DCTSIZE) && (x*DCTSIZE) <(cols-DCTSIZE)) {
    tmp /= (short8)(64, 64, 64, 64, 64, 64, 64, 64);
  } else if (y<DCTSIZE && (x*DCTSIZE)<DCTSIZE) //UL
  {
    tmp /= MatUL[y%DCTSIZE];
  } else if (y>=(rows-DCTSIZE) && (x*DCTSIZE)<DCTSIZE) //DL
  {
    tmp /= MatDL[y%DCTSIZE];
  } else if (y>=(rows-DCTSIZE) && (x*DCTSIZE)>=(cols-DCTSIZE)) //DR
  {
    tmp /= MatDR[y%DCTSIZE];
  } else if (y<DCTSIZE && (x*DCTSIZE)>=(cols-DCTSIZE)) //UR
  {
    tmp /= MatUR [y%DCTSIZE];
  } else if (y<DCTSIZE && (x*DCTSIZE)>=DCTSIZE && (x*DCTSIZE) <(cols-DCTSIZE)) //up side
  {
    tmp /= (short8)(((y%DCTSIZE)+1)*DCTSIZE);
  } else if (y>=(rows-DCTSIZE) && (x*DCTSIZE)>=DCTSIZE && (x*DCTSIZE) <(cols-DCTSIZE)) //down sise
  {
    tmp /= (short8)((DCTSIZE-(y%DCTSIZE))*DCTSIZE);
  } else if ((x*DCTSIZE)>=(cols-DCTSIZE) && y>=(DCTSIZE) && y<(rows-DCTSIZE)) //right side
  {
    tmp /= MatUR[DCTSIZE-1];
  } else if ((x*DCTSIZE)<DCTSIZE && y>=(DCTSIZE) && y<(rows-DCTSIZE)) // left side
  {
    tmp /= MatUL[DCTSIZE-1];
  }
  
  
  //resDCT[y*cols+x] = convert_uchar_sat(tmp/64);
  *(global short8*)&resDCT[y*cols+x*DCTSIZE] = tmp;
}

