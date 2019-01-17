#define GETJSAMPLE(value)  (value)
#define DEQUANTIZE(coef,quantval)  ((coef) * (quantval))
#define DCTSIZE 8
#define DCTSIZE2 64
#define CENTERJSAMPLE  128
#define LOCAL_WINDOW_SIZE 16
#define MAXJSAMPLE  255
#define RANGE_MASK  (MAXJSAMPLE * 4 + 3)
#define RANGE_LIMIT(x) (x)



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
  
  
  
  
  /* Pass 1: process rows. */
  
  
#pragma unroll
  for (int ctr = 0; ctr < DCTSIZE; ctr++) {
    //elemptr = sample_data[ctr] + start_col;
    float *dataptr = &data[ctr*DCTSIZE];
    local float* elemptr = &sample_data[ctr*stride];
    
    /* Load data into workspace */
    tmp0 = (GETJSAMPLE(elemptr[0]) + GETJSAMPLE(elemptr[7]));
    tmp7 = (GETJSAMPLE(elemptr[0]) - GETJSAMPLE(elemptr[7]));
    tmp1 = (GETJSAMPLE(elemptr[1]) + GETJSAMPLE(elemptr[6]));
    tmp6 = (GETJSAMPLE(elemptr[1]) - GETJSAMPLE(elemptr[6]));
    tmp2 = (GETJSAMPLE(elemptr[2]) + GETJSAMPLE(elemptr[5]));
    tmp5 = (GETJSAMPLE(elemptr[2]) - GETJSAMPLE(elemptr[5]));
    tmp3 = (GETJSAMPLE(elemptr[3]) + GETJSAMPLE(elemptr[4]));
    tmp4 = (GETJSAMPLE(elemptr[3]) - GETJSAMPLE(elemptr[4]));
    
    /* Even part */
    
    tmp10 = tmp0 + tmp3;  /* phase 2 */
    tmp13 = tmp0 - tmp3;
    tmp11 = tmp1 + tmp2;
    tmp12 = tmp1 - tmp2;
    
    /* Apply unsigned->signed conversion */
    dataptr[0] = tmp10 + tmp11 - 8 * CENTERJSAMPLE; /* phase 3 */
    dataptr[4] = tmp10 - tmp11;
    
    z1 = (tmp12 + tmp13) * 0.707106781f; /* c4 */
    dataptr[2] = tmp13 + z1;  /* phase 5 */
    dataptr[6] = tmp13 - z1;
    
    /* Odd part */
    
    tmp10 = tmp4 + tmp5;  /* phase 2 */
    tmp11 = tmp5 + tmp6;
    tmp12 = tmp6 + tmp7;
    
    /* The rotator is modified from fig 4-8 to avoid extra negations. */
    z5 = (tmp10 - tmp12) * 0.382683433f; /* c6 */
    z2 = 0.541196100f * tmp10 + z5; /* c2-c6 */
    z4 = 1.306562965f * tmp12 + z5; /* c2+c6 */
    z3 = tmp11 * 0.707106781f; /* c4 */
    
    z11 = tmp7 + z3;    /* phase 5 */
    z13 = tmp7 - z3;
    
    dataptr[5] = z13 + z2;  /* phase 6 */
    dataptr[3] = z13 - z2;
    dataptr[1] = z11 + z4;
    dataptr[7] = z11 - z4;
    
  }
  
  /* Pass 2: process columns. */
  
#pragma unroll
  for (int ctr = 0; ctr < DCTSIZE; ctr++) {
    float *dataptr = &data[ctr];
    
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
    
    z1 = (tmp12 + tmp13) * 0.707106781f; /* c4 */
    dataptr[DCTSIZE*2] = tmp13 + z1; /* phase 5 */
    dataptr[DCTSIZE*6] = tmp13 - z1;
    
    /* Odd part */
    
    tmp10 = tmp4 + tmp5;  /* phase 2 */
    tmp11 = tmp5 + tmp6;
    tmp12 = tmp6 + tmp7;
    
    /* The rotator is modified from fig 4-8 to avoid extra negations. */
    z5 = (tmp10 - tmp12) * 0.382683433f; /* c6 */
    z2 = 0.541196100f * tmp10 + z5; /* c2-c6 */
    z4 = 1.306562965f * tmp12 + z5; /* c2+c6 */
    z3 = tmp11 * 0.707106781f; /* c4 */
    
    z11 = tmp7 + z3;    /* phase 5 */
    z13 = tmp7 - z3;
    
    dataptr[DCTSIZE*5] = z13 + z2; /* phase 6 */
    dataptr[DCTSIZE*3] = z13 - z2;
    dataptr[DCTSIZE*1] = z11 + z4;
    dataptr[DCTSIZE*7] = z11 - z4;
    
  }
}


void my_jpeg_idct_float (const float constant *quantptr_, float *coef_block)
{
  float tmp0, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7;
  float tmp10, tmp11, tmp12, tmp13;
  float z5, z10, z11, z12, z13;
  //FLOAT_MULT_TYPE * quantptr;
  //JSAMPLE *range_limit = cinfo->sample_range_limit;
  
  /* Pass 1: process columns from input, store into work array. */
  
  //quantptr = (FLOAT_MULT_TYPE *) compptr->dct_table;
  
#pragma unroll
  for (int ctr = 0; ctr < DCTSIZE; ctr++) {
    float *wsptr = &coef_block[ctr];
    const float constant *quantptr = quantptr_ + ctr;
    /* Even part */
    
    tmp0 = DEQUANTIZE(wsptr[DCTSIZE*0], quantptr[DCTSIZE*0]);
    tmp1 = DEQUANTIZE(wsptr[DCTSIZE*2], quantptr[DCTSIZE*2]);
    tmp2 = DEQUANTIZE(wsptr[DCTSIZE*4], quantptr[DCTSIZE*4]);
    tmp3 = DEQUANTIZE(wsptr[DCTSIZE*6], quantptr[DCTSIZE*6]);
    
    tmp10 = tmp0 + tmp2;  /* phase 3 */
    tmp11 = tmp0 - tmp2;
    
    tmp13 = tmp1 + tmp3;  /* phases 5-3 */
    tmp12 = (tmp1 - tmp3) * 1.414213562f - tmp13; /* 2*c4 */
    
    tmp0 = tmp10 + tmp13;  /* phase 2 */
    tmp3 = tmp10 - tmp13;
    tmp1 = tmp11 + tmp12;
    tmp2 = tmp11 - tmp12;
    
    /* Odd part */
    
    tmp4 = DEQUANTIZE(wsptr[DCTSIZE*1], quantptr[DCTSIZE*1]);
    tmp5 = DEQUANTIZE(wsptr[DCTSIZE*3], quantptr[DCTSIZE*3]);
    tmp6 = DEQUANTIZE(wsptr[DCTSIZE*5], quantptr[DCTSIZE*5]);
    tmp7 = DEQUANTIZE(wsptr[DCTSIZE*7], quantptr[DCTSIZE*7]);
    
    z13 = tmp6 + tmp5;    /* phase 6 */
    z10 = tmp6 - tmp5;
    z11 = tmp4 + tmp7;
    z12 = tmp4 - tmp7;
    
    tmp7 = z11 + z13;    /* phase 5 */
    tmp11 = (z11 - z13) * 1.414213562f; /* 2*c4 */
    
    z5 = (z10 + z12) * 1.847759065f; /* 2*c2 */
    tmp10 = z5 - z12 * 1.082392200f; /* 2*(c2-c6) */
    tmp12 = z5 - z10 * 2.613125930f; /* 2*(c2+c6) */
    
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
    
  }
  
  /* Pass 2: process rows from work array, store into output array. */
  
  
#pragma unroll
  for (int ctr = 0; ctr < DCTSIZE; ctr++) {
    /* Rows of zeroes can be exploited in the same way as we did with columns.
     * However, the column calculation has created many nonzero AC terms, so
     * the simplification applies less often (typically 5% to 10% of the time).
     * And testing floats for zero is relatively expensive, so we don't bother.
     */
    
    /* Even part */
    
    /* Apply signed->unsigned and prepare float->int conversion */
    float *wsptr = &coef_block[ctr*DCTSIZE];
    z5 = wsptr[0] + (CENTERJSAMPLE + 0.5f);
    tmp10 = z5 + wsptr[4];
    tmp11 = z5 - wsptr[4];
    
    tmp13 = wsptr[2] + wsptr[6];
    tmp12 = (wsptr[2] - wsptr[6]) * 1.414213562f - tmp13;
    
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
    tmp11 = (z11 - z13) * 1.414213562f;
    
    z5 = (z10 + z12) * 1.847759065f; /* 2*c2 */
    tmp10 = z5 - z12 * 1.082392200f; /* 2*(c2-c6) */
    tmp12 = z5 - z10 * 2.613125930f; /* 2*(c2+c6) */
    
    tmp6 = tmp12 - tmp7;
    tmp5 = tmp11 - tmp6;
    tmp4 = tmp10 - tmp5;
    
    /* Final output stage: float->int conversion and range-limit */
    
    wsptr[0] = RANGE_LIMIT(tmp0 + tmp7);
    wsptr[7] = RANGE_LIMIT(tmp0 - tmp7);
    wsptr[1] = RANGE_LIMIT(tmp1 + tmp6);
    wsptr[6] = RANGE_LIMIT(tmp1 - tmp6);
    wsptr[2] = RANGE_LIMIT(tmp2 + tmp5);
    wsptr[5] = RANGE_LIMIT(tmp2 - tmp5);
    wsptr[3] = RANGE_LIMIT(tmp3 + tmp4);
    wsptr[4] = RANGE_LIMIT(tmp3 - tmp4);
  }
}

void strideFloatDCT(local float *sample_data, const unsigned int stride, float *data,  const float constant *scaleMatQuant)
{
  
  my_stride_jpeg_fdct_float(data, sample_data, stride);
  
#pragma unroll
  for (int i = 0; i < 8; ++i) {
#pragma unroll
    for (int j = 0; j < 8; ++j) {
      data[i * 8 + j] = rint(data[i * 8 + j] * scaleMatQuant[i * 8 + j]);
    }
  }
  my_jpeg_idct_float(&scaleMatQuant[DCTSIZE2], data);
}


kernel __attribute__((reqd_work_group_size(64, 1, 1))) void reapplicationJPEG(const uchar global *rawImg, const uint rows, const uint cols, short global *resDCT,  constant const float  *scaleMatQuant, const uint offset, const uint stride)
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
  strideFloatDCT(&mat.mat_1D[y*LOCAL_WINDOW_SIZE+x], LOCAL_WINDOW_SIZE, res, scaleMatQuant);
  
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

kernel void dctRigthEdge(const uchar global *rawImg, const uint rows, const uint cols, short global *resDCT,  constant const float  *scaleMatQuant, const uint offset, const uint stride)
{
  local union{
    float mat_1D[LOCAL_WINDOW_SIZE*DCTSIZE];
    float mat_2D[LOCAL_WINDOW_SIZE][DCTSIZE];
    float8 mat8[LOCAL_WINDOW_SIZE];
  } mat;
  int x = get_local_id(0);
  
  mat.mat8[x] = convert_float8(*(global uchar8*)&rawImg[(get_group_id(0)*LOCAL_WINDOW_SIZE+x)*stride+cols-DCTSIZE+offset]);
  barrier(CLK_LOCAL_MEM_FENCE);
  mat.mat8[DCTSIZE+x] = convert_float8(*(global uchar8*)&rawImg[(get_group_id(0)*LOCAL_WINDOW_SIZE+DCTSIZE+x)*stride+cols-DCTSIZE+offset]);
  barrier(CLK_LOCAL_MEM_FENCE);
  
  float res[DCTSIZE2];
  strideFloatDCT(&mat.mat_1D[x*DCTSIZE], DCTSIZE, res, scaleMatQuant);
  barrier(CLK_LOCAL_MEM_FENCE);
  
  mat.mat8[x] = 0.0f;
  barrier(CLK_LOCAL_MEM_FENCE);
  mat.mat8[DCTSIZE+x] = 0.0f;
  barrier(CLK_LOCAL_MEM_FENCE);
  
  for (int i=0; i<DCTSIZE; i++)
  {
    for (int j=0; j<DCTSIZE; j++)
    {
      mat.mat_2D[x+i][j] += res[i*DCTSIZE+j];
      barrier(CLK_LOCAL_MEM_FENCE);
    }
  }
  
  *(global short8*)&resDCT[(get_group_id(0)*LOCAL_WINDOW_SIZE+x)*stride+cols-DCTSIZE+offset] += convert_short8_sat_rte(mat.mat8[x]);
  barrier(CLK_LOCAL_MEM_FENCE);
  *(global short8*)&resDCT[(get_group_id(0)*LOCAL_WINDOW_SIZE+DCTSIZE+x)*stride+cols-DCTSIZE+offset] +=  convert_short8_sat_rte(mat.mat8[DCTSIZE+x]);
}

kernel void dctDownEdge1(const uchar global *rawImg, const uint rows, const uint cols, short global *resDCT,  constant const float  *scaleMatQuant, const uint offset, const uint stride)
{
  local union{
    float mat_1D[DCTSIZE*LOCAL_WINDOW_SIZE];
    float mat_2D[DCTSIZE][LOCAL_WINDOW_SIZE];
    float8 mat8[DCTSIZE][LOCAL_WINDOW_SIZE/DCTSIZE];
  } mat;
  int x = get_local_id(0);
  
  mat.mat8[x][0] = convert_float8(*(global uchar8*)&rawImg[(rows-DCTSIZE+x)*stride+(get_group_id(0)*LOCAL_WINDOW_SIZE)+offset]);
  barrier(CLK_LOCAL_MEM_FENCE);
  mat.mat8[x][1] = convert_float8(*(global uchar8*)&rawImg[(rows-DCTSIZE+x)*stride+(get_group_id(0)*LOCAL_WINDOW_SIZE+DCTSIZE+offset)]);
  barrier(CLK_LOCAL_MEM_FENCE);
  
  float res[DCTSIZE2];
  strideFloatDCT(&mat.mat_1D[x], DCTSIZE, res, scaleMatQuant);
  barrier(CLK_LOCAL_MEM_FENCE);
  
  mat.mat8[x][0] = 0.0f;
  barrier(CLK_LOCAL_MEM_FENCE);
  mat.mat8[x][1] = 0.0f;
  barrier(CLK_LOCAL_MEM_FENCE);
  
  for (int i=0; i<DCTSIZE; i++)
  {
    for (int j=0; j<DCTSIZE; j++)
    {
      mat.mat_2D[i][x+j] += res[i*DCTSIZE+j];
      barrier(CLK_LOCAL_MEM_FENCE);
    }
  }
  
  *(global short8*)&resDCT[(rows-DCTSIZE+x)*stride+(get_group_id(0)*LOCAL_WINDOW_SIZE+offset)] += convert_short8_sat_rte(mat.mat8[x][0]);
  barrier(CLK_LOCAL_MEM_FENCE);
  *(global short8*)&resDCT[(rows-DCTSIZE+x)*stride+(get_group_id(0)*LOCAL_WINDOW_SIZE+DCTSIZE+offset)] +=  convert_short8_sat_rte(mat.mat8[x][1]);
}

kernel void dctDownEdge(const uchar global *rawImg, const uint rows, const uint cols, short global *resDCT,  constant const float  *scaleMatQuant, const uint offset, const uint stride)
{
  local union{
    float mat_1D[DCTSIZE*LOCAL_WINDOW_SIZE];
    float mat_2D[DCTSIZE][LOCAL_WINDOW_SIZE];
    float16 mat16[DCTSIZE];
  } mat;
  int x = get_local_id(0);
  
  mat.mat16[x] = convert_float16(*(global uchar16*)&rawImg[(rows-DCTSIZE+x)*stride+(get_group_id(0)*LOCAL_WINDOW_SIZE)+offset]);
  barrier(CLK_LOCAL_MEM_FENCE);
  
  float res[DCTSIZE2];
  strideFloatDCT(&mat.mat_1D[x], LOCAL_WINDOW_SIZE, res, scaleMatQuant);
  barrier(CLK_LOCAL_MEM_FENCE);
  
  mat.mat16[x] = 0.0f;
  barrier(CLK_LOCAL_MEM_FENCE);
  
  for (int i=0; i<DCTSIZE; i++)
  {
    for (int j=0; j<DCTSIZE; j++)
    {
      mat.mat_2D[i][x+j] += res[i*DCTSIZE+j];
      barrier(CLK_LOCAL_MEM_FENCE);
    }
  }
  
  *(global short16*)&resDCT[(rows-DCTSIZE+x)*stride+(get_group_id(0)*LOCAL_WINDOW_SIZE+offset)] += convert_short16_sat_rte(mat.mat16[x]);
}

kernel void dctCorner(const uchar global *rawImg, const uint rows, const uint cols, short global *resDCT,  constant const float  *scaleMatQuant, const uint stride)
{
  local union{
    float mat_1D[DCTSIZE2];
    float mat_2D[DCTSIZE][DCTSIZE];
    float8 mat8[DCTSIZE];
  } mat;
  int x = get_local_id(0);
  
  mat.mat8[x] = convert_float8(*(global uchar8*)&rawImg[(rows-DCTSIZE+x)*stride+cols-DCTSIZE]);
  barrier(CLK_LOCAL_MEM_FENCE);
  
  float res[DCTSIZE2];
  if (x == 0)
  {
    strideFloatDCT(mat.mat_1D, DCTSIZE, res, scaleMatQuant);
  }
  barrier(CLK_LOCAL_MEM_FENCE);
  
  mat.mat8[x] = 0.0f;
  barrier(CLK_LOCAL_MEM_FENCE);
  
  if (x == 0) {
    for (int i = 0; i<DCTSIZE2; i++){
      mat.mat_1D[i] = res[i];
    }
  }
  
  *(global short8*)&resDCT[(rows-DCTSIZE+x)*stride+cols-DCTSIZE] += convert_short8_sat_rte(mat.mat8[x]);
}


kernel void justFloatDCT(const uchar global *rawImg, const uint rows, const uint cols, short global *resDCT, float constant *scaleMatQuant, float constant *dctTable, const uint offset, const uint stride)
{
  local uchar mat[DCTSIZE2];
  for (int i = 0; i<DCTSIZE2; i++)
    mat[i] = rawImg[i];
  float res[DCTSIZE2];
  //strideFloatDCT(mat, DCTSIZE, res, scaleMatQuant, dctTable);
  for (int i = 0; i<DCTSIZE; i++)
    for (int j = 0; j<DCTSIZE; j++)
      resDCT[i*DCTSIZE+j] = res[i*DCTSIZE+j];
}



__constant const float8 MatUL[DCTSIZE] = {
  (float8)(1/1.0f,  1/2.0f,  1/3.0f,  1/4.0f,  1/5.0f,  1/6.0f,  1/7.0f,  1/8.0f),
  (float8)(1/2.0f,  1/4.0f,  1/6.0f,  1/8.0f,  1/10.0f,  1/12.0f,  1/14.0f,  1/16.0f),
  (float8)(1/3.0f,  1/6.0f,  1/9.0f,  1/12.0f,  1/15.0f,  1/18.0f,  1/21.0f,  1/24.0f),
  (float8)(1/4.0f,  1/8.0f,  1/12.0f,  1/16.0f,  1/20.0f,  1/24.0f,  1/28.0f,  1/32.0f),
  (float8)(1/5.0f,  1/10.0f,  1/15.0f,  1/20.0f,  1/25.0f,  1/30.0f,  1/35.0f,  1/40.0f),
  (float8)(1/6.0f,  1/12.0f,  1/18.0f,  1/24.0f,  1/30.0f,  1/36.0f,  1/42.0f,  1/48.0f),
  (float8)(1/7.0f,  1/14.0f,  1/21.0f,  1/28.0f,  1/35.0f,  1/42.0f,  1/49.0f,  1/56.0f),
  (float8)(1/8.0f,  1/16.0f,  1/24.0f,  1/32.0f,  1/40.0f,  1/48.0f,  1/56.0f,  1/64.0f)
};

__constant const float8 MatDL[DCTSIZE] = {
  (float8)(1/8.0f,  1/16.0f,  1/24.0f,  1/32.0f,  1/40.0f,  1/48.0f,  1/56.0f,  1/64.0f),
  (float8)(1/7.0f,  1/14.0f,  1/21.0f,  1/28.0f,  1/35.0f,  1/42.0f,  1/49.0f,  1/56.0f),
  (float8)(1/6.0f,  1/12.0f,  1/18.0f,  1/24.0f,  1/30.0f,  1/36.0f,  1/42.0f,  1/48.0f),
  (float8)(1/5.0f,  1/10.0f,  1/15.0f,  1/20.0f,  1/25.0f,  1/30.0f,  1/35.0f,  1/40.0f),
  (float8)(1/4.0f,  1/8.0f,  1/12.0f,  1/16.0f,  1/20.0f,  1/24.0f,  1/28.0f,  1/32.0f),
  (float8)(1/3.0f,  1/6.0f,  1/9.0f,  1/12.0f,  1/15.0f,  1/18.0f,  1/21.0f,  1/24.0f),
  (float8)(1/2.0f,  1/4.0f,  1/6.0f,  1/8.0f,  1/10.0f,  1/12.0f,  1/14.0f,  1/16.0f),
  (float8)(1/1.0f,  1/2.0f,  1/3.0f,  1/4.0f,  1/5.0f,  1/6.0f,  1/7.0f,  1/8.0f)
};

__constant const float8 MatUR[DCTSIZE] = {
  (float8)(1/8.0f,  1/7.0f,  1/6.0f,  1/5.0f,  1/4.0f,  1/3.0f,  1/2.0f,  1/1.0f),
  (float8)(1/16.0f,  1/14.0f,  1/12.0f,  1/10.0f,  1/8.0f,  1/6.0f,  1/4.0f,  1/2.0f),
  (float8)(1/24.0f,  1/21.0f,  1/18.0f,  1/15.0f,  1/12.0f,  1/9.0f,  1/6.0f,  1/3.0f),
  (float8)(1/32.0f,  1/28.0f,  1/24.0f,  1/20.0f,  1/16.0f,  1/12.0f,  1/8.0f,  1/4.0f),
  (float8)(1/40.0f,  1/35.0f,  1/30.0f,  1/25.0f,  1/20.0f,  1/15.0f,  1/10.0f,  1/5.0f),
  (float8)(1/48.0f,  1/42.0f,  1/36.0f,  1/30.0f,  1/24.0f,  1/18.0f,  1/12.0f,  1/6.0f),
  (float8)(1/56.0f,  1/49.0f,  1/42.0f,  1/35.0f,  1/28.0f,  1/21.0f,  1/14.0f,  1/7.0f),
  (float8)(1/64.0f,  1/56.0f,  1/48.0f,  1/40.0f,  1/32.0f,  1/24.0f,  1/16.0f,  1/8.0f)
};

__constant const float8 MatDR[DCTSIZE] = {
  (float8)(1/64.0f,  1/56.0f,  1/48.0f,  1/40.0f,  1/32.0f,  1/24.0f,  1/16.0f,  1/8.0f),
  (float8)(1/56.0f,  1/49.0f,  1/42.0f,  1/35.0f,  1/28.0f,  1/21.0f,  1/14.0f,  1/7.0f),
  (float8)(1/48.0f,  1/42.0f,  1/36.0f,  1/30.0f,  1/24.0f,  1/18.0f,  1/12.0f,  1/6.0f),
  (float8)(1/40.0f,  1/35.0f,  1/30.0f,  1/25.0f,  1/20.0f,  1/15.0f,  1/10.0f,  1/5.0f),
  (float8)(1/32.0f,  1/28.0f,  1/24.0f,  1/20.0f,  1/16.0f,  1/12.0f,  1/8.0f,  1/4.0f),
  (float8)(1/24.0f,  1/21.0f,  1/18.0f,  1/15.0f,  1/12.0f,  1/9.0f,  1/6.0f,  1/3.0f),
  (float8)(1/16.0f,  1/14.0f,  1/12.0f,  1/10.0f,  1/8.0f,  1/6.0f,  1/4.0f,  1/2.0f),
  (float8)(1/8.0f,  1/7.0f,  1/6.0f,  1/5.0f,  1/4.0f,  1/3.0f,  1/2.0f,  1/1.0f)
};

kernel void divMat(short global *resDCT, uchar global *resImg, const uint rows, const uint cols)
{
  int x = get_global_id(0);
  int y = get_global_id(1);
  
  
  float8 tmp = convert_float8(*(global short8*)&resDCT[y*cols+x*DCTSIZE]);
  
  if (y>=DCTSIZE && (x*DCTSIZE) >= DCTSIZE && y<(rows-DCTSIZE) && (x*DCTSIZE) <(cols-DCTSIZE)) {
    tmp *= (float8)(1/64.0f,  1/64.0f,  1/64.0f,  1/64.0f,  1/64.0f,  1/64.0f,  1/64.0f,  1/64.0f);
  } else if (y<DCTSIZE && (x*DCTSIZE)<DCTSIZE) //UL
  {
    tmp *= MatUL[y%DCTSIZE];
  } else if (y>=(rows-DCTSIZE) && (x*DCTSIZE)<DCTSIZE) //DL
  {
    tmp *= MatDL[y%DCTSIZE];
  } else if (y>=(rows-DCTSIZE) && (x*DCTSIZE)>=(cols-DCTSIZE)) //DR
  {
    tmp *= MatDR[y%DCTSIZE];
  } else if (y<DCTSIZE && (x*DCTSIZE)>=(cols-DCTSIZE)) //UR
  {
    tmp *= MatUR [y%DCTSIZE];
  } else if (y<DCTSIZE && (x*DCTSIZE)>=DCTSIZE && (x*DCTSIZE) <(cols-DCTSIZE)) //up side
  {
    tmp *= (float8)(1.0f/(((y%DCTSIZE)+1)*DCTSIZE));
  } else if (y>=(rows-DCTSIZE) && (x*DCTSIZE)>=DCTSIZE && (x*DCTSIZE) <(cols-DCTSIZE)) //down sise
  {
    tmp *= (float8)(1.0f/((DCTSIZE-(y%DCTSIZE))*DCTSIZE));
  } else if ((x*DCTSIZE)>=(cols-DCTSIZE) && y>=(DCTSIZE) && y<(rows-DCTSIZE)) //right side
  {
    tmp *= MatUR[DCTSIZE-1];
  } else if ((x*DCTSIZE)<DCTSIZE && y>=(DCTSIZE) && y<(rows-DCTSIZE)) // left side
  {
    tmp *= MatUL[DCTSIZE-1];
  }
  
  
  //resDCT[y*cols+x] = convert_uchar_sat(tmp/64);
  *(global uchar8*)&resImg[y*cols+x*DCTSIZE] = convert_uchar8_sat_rte(tmp);
}

