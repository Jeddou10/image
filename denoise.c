#include <stdlib.h>
#include <stdint.h>
#include <pthread.h>
#include <immintrin.h>

// data to be given to thread's methods
struct ThreadData
{
  uint8_t *Q;
  float *Ql;
  float *Qw;
  uint8_t *result;
  size_t width;
  size_t height;
}
;

// (producer )method to calculate Ql
void *
calculate_Ql (void *arg)
{
  struct ThreadData *data = (struct ThreadData *) arg;
  size_t width = data->width;
  size_t height = data->height;
  uint8_t *Q = data->Q;
  float *Ql = data->Ql;

  float Ml[9] = { 0, 1, 0, 1, -4, 1, 0, 1, 0 };
  for (size_t y = 0; y < height; y++)
    {
      for (size_t x = 0; x < width; x++)
	{
	  float res = 0.0;
	  for (size_t j = 0; j <= 2; j++)
	    {
	      for (size_t i = 0; i <= 2; i++)
		{
		  if (x + i  < width + 1 && x + i  >= 1 && y + j  >= 1
		      && y + j < height + 1)
		    {
		      res +=
			Ml[(i) * 3 + ( j )] * ((float)( Q[x + i - 1 + (y + j - 1) * width]));
		    }
		}
	    }
	  Ql[y * width + x] =( res>0 ? (res/ (float) (1020)) :( (-res) / (float) (1020))
	  );
	}
    }
  return NULL;
}

//(producer ) method to calculate Qw
void * calculate_Qw (void *arg)
{
  struct ThreadData *data = ((struct ThreadData *) arg);
  size_t width = data->width;
  size_t height = data->height;
  uint8_t *Q = data->Q;
  float *Qw = data->Qw;

  float Mw[9] = { 1, 2, 1, 2, 4, 2, 1, 2, 1 };
  for (size_t y = 0; y < height; y++)
    {
      for (size_t x = 0; x < width; x++)
	{
	  float res = 0.0;
	  for (size_t j = 0; j <= 2; j++)
	    {
	      for (size_t i = 0; i <= 2; i++)
		{
		  if (x + i < width + 1 && x + i >= 1 && y + j >= 1
		      && y + j < height + 1)
		    {
		      res +=
			Mw[( i) * 3 + ( j)] * ((float)(Q[x + i -1 + (y + j - 1) * width]));
		    }
		}
	    }
	  Qw[y * width + x] =  (res/ (float) 16);
	}
    }
  return NULL;
}
//(consumer) method to calculate the result 
void * calculate_result (void *arg)
{
  struct ThreadData *data = ((struct ThreadData *) arg);
  size_t width = data->width;
  size_t height = data->height;
  uint8_t *Q = data->Q;
  float  *Ql = data->Ql;
  float  *Qw = data->Qw;
  uint8_t *result = data->result;
  
  for (size_t index = 0; index < width * height; index++)
    {
      // waiting for calcculating Ql and Qw
      while (Ql[index] == -1 || Qw[index] == -1)
	{
	}

      float Qwxy = (float)Qw[index];
      float Qlxy1020 =(float)Ql[index] ;
      float Qprime = Qlxy1020 * ((float)Q[index] )+ (1 - Qlxy1020) * Qwxy;
      result[index] = (uint8_t) Qprime;
    }
  return NULL;
}

// tmp2 : width*height floats 
// tmp1 : with*height uint8_ts
void denoise_threading (const uint8_t * img, size_t width,
		   size_t height, float a, float b,
		   float c, uint8_t * tmp1, uint8_t * tmp2, uint8_t * result)
{
	
  size_t pixels_num = width * height;
  float  Qw[pixels_num];
  for(size_t i = 0  ; i < pixels_num ; i++){
	((float*)tmp2)[i] = -1 ;
	Qw[i] = -1 ;
  }
  pthread_t threads[3];
  float sum = a + b + c;
  for (size_t i = 0; i < width * height * 3; i += 3)
    {
      tmp1[i / 3] =(uint8_t)
	((a * img[i] + b * img[i + 1] + c * img[i + 2]) / sum);
    }
  struct ThreadData data = (struct ThreadData)
  {
    .Q =  tmp1,
    .Ql = (float*)tmp2,
    .Qw = Qw,
    .result = result,
    .width = width,
    .height = height
  };

  pthread_create (&threads[0], NULL, calculate_Ql, &data);
  pthread_create (&threads[1], NULL, calculate_Qw, &data);
  pthread_create (&threads[2], NULL, calculate_result, &data);
  pthread_join (threads[0], NULL);
  pthread_join (threads[1], NULL);
  pthread_join (threads[2], NULL);
}

void
denoise_optimised (const uint8_t * img, size_t width,
	  size_t height, float a, float b,
	  float c, uint8_t * tmp1, uint8_t * tmp2, uint8_t * result)
{

  float sum = a + b + c;
  for (size_t i = 0; i < width * height * 3; i += 3)
    {
       tmp1[i / 3] =(uint8_t)
	((a * img[i] + b * img[i + 1] + c * img[i + 2]) / sum);
    }

  float Ml[9] = { 0, 1, 0, 1, -4, 1, 0, 1, 0 };
  float Mw[9] = { 1, 2, 1, 2, 4, 2, 1, 2, 1 };
  size_t adress = 0;

  for (size_t y = 0; y < height; y++)
    {
      for (size_t x = 0; x < width; x++)
	{
	  float res1 = 0.0;
	  float res2 = 0.0;
	  float d = 0.0;
	  for (size_t j = 0; j <= 2; j++)
	    {

	      for (int i = 0; i <= 2; i++)
		{

		  if (x + i < width + 1 &&  x + i>= 1 &&  y + j >= 1
		      && y + j < height + 1)
		    {
		      d = tmp1[x + i - 1 + (y + j - 1) * width];
		      adress = 3*i + j;

		      res1 += Ml[adress] * d;
		      res2 += Mw[adress] * d;

		    }


		}
	    }
		 
	  d = tmp1[x + y * width];
	  res1 = res1 > 0 ? (res1 / (float) 1020 ): ((-res1) / (float) 1020) ;
	  if(x == 100 ) {
		//printf("x : %d -- y : %d -- Ql : %f\n" ,x,y, res1);
	  }
	  res2 /= (float) 16;

	  float Qprime = res1 * d + (1 - res1) * res2;

	  result[y * width + x ] = (uint8_t) (Qprime);

    }}

}

void
denoise_naive (const uint8_t * img, size_t width,
	  size_t height, float a, float b,
	  float c, uint8_t * tmp1, uint8_t * tmp2, uint8_t * result)
{
	
  float sum = a + b + c;
  float Ml[9] = { 0, 1, 0, 1, -4, 1, 0, 1, 0 };
  float Mw[9] = { 1, 2, 1, 2, 4, 2, 1, 2, 1 };
  size_t k = 0;
  //go through the pixels 
  for (size_t  y = 0; y < height; y++)
    {
      for (size_t x = 0; x < width; x++)
	{
	  k = 3 * (y * width + x);

	  //calculate the grayscale of the pixel
	  float D = (a * img[k] + b * img[k + 1] + c * img[k + 2]) / sum;

	  // calculate Ql for the pixel (x,y)
	  float Ql = 0.0;
	  for (size_t j = 0; j <= 2; j++)
	    {
	      for (size_t i = 0; i <= 2; i++)
		{
		  if (x + i < width + 1 && (x + i) >= 1 && (y + j) >= 1
		      && y + j < height + 1)
		    {
		      k = 3 * (x + i - 1 + (y + j - 1) * width);
		      Ql +=
			Ml[(i) * 3 +
			   (j)] * ((a * img[k] + b * img[k + 1] +
					c * img[k + 2]) / sum);
		    }

		}
	    }

	  //calculate Qw for the pixel (X,y)

	  float Qw = 0.0;
	  for (size_t j = 0; j <= 2; j++)
	    {
	      for (size_t i = 0; i <= 2; i++)
		{
		  if (x + i < width + 1 && x + i >= 1 && (y + j) >= 1
		      && y + j < height + 1)
		    {
		      k = 3 * (x + i - 1 + (y + j - 1) * width);
		      Qw +=
			Mw[(i) * 3 +
			   (j)] * ((a * img[k] + b * img[k + 1] +
					c * img[k + 2]) / sum);
		    }

		}
	    }

	  // calculate Q' for pixel (x,y)
	  Ql = Ql>0 ? (Ql / (float) 1020) : (-(Ql / (float) 1020));
	  Qw = (Qw / (float) 16);
	  uint8_t Qprime =(uint8_t) (Ql * D + (1 - Ql) * Qw);
	  result[y * width + x] = Qprime;


	}
    }
}




void
denoise_simd (const uint8_t * img, size_t width,
	      size_t height, float a, float b,
	      float c, uint8_t * tmp1, uint8_t * tmp2, uint8_t * result)
{
	
  // counter used to fill tmp1 
  size_t k = 0;

  // fill the first line of tmp1 all with 0
 /* for (k; k <=  width; k++){
    ((float *) tmp1)[k] = 0;
  }*/
  // fill tmp1 with d (the first column of tmp1 is all with 0)
  float sum = a + b + c;
  for (size_t i = 0; i < (width * height * 3); i += 3)
    {
      // fill the first column only  with 0
      if (k % (width + 1) == 0)
	{
	  ((float *) tmp1)[k] = 0;
	  k++;
	}
      // calculate d and fill tmp1 with it 
      ((float *) tmp1)[k] =(float)(
	((a * img[i] + b * img[i + 1] + c * img[i + 2]) / sum));
      k++;
    }

  // add a line at the end of tmp1 all with 0
  for (size_t i = k; i <= (k + width); i++)
    ((float *) tmp1)[i] = 0;

  // fill Ml and Mw and add a 0 after each 3 elements
  float Ml[12] = { 0, 1, 0, 0, 1, -4, 1, 0, 0, 1, 0, 0 };
  float Mw[12] = { 1, 2, 1, 0, 2, 4, 2, 0, 1, 2, 1, 0 };

  // load each 4 elements of Ml and Mw together using simd
  __m128 Mli1 = _mm_loadu_ps (&Ml[0]);
  __m128 Mli2 = _mm_loadu_ps (&Ml[4]);
  __m128 Mli3 = _mm_loadu_ps (&Ml[8]);
  __m128 Mwi1 = _mm_loadu_ps (&Mw[0]);
  __m128 Mwi2 = _mm_loadu_ps (&Mw[4]);
  __m128 Mwi3 = _mm_loadu_ps (&Mw[8]);

  // prepare some __m128 variables to use them later
  __m128 di1 = _mm_setzero_ps ();
  __m128 di2 = _mm_setzero_ps ();
  __m128 di3 = _mm_setzero_ps ();
  __m128 mul1 = _mm_setzero_ps ();
  __m128 mul2 = _mm_setzero_ps ();
  __m128 mul3 = _mm_setzero_ps ();
  __m128 add = _mm_setzero_ps ();

  // iterate over all pixels
  for (size_t y = 0; y < height; y++)
    {

      for (size_t x = 0; x < width; x++)
	{

	  // load each 3 elements from tmp1 in di1,di2 and di3 (the 4 th element will be multiplied by 0 in the next steps and ignored)
	  di1 =
	    _mm_loadu_ps (&((float *) tmp1)
			  [1 + x - 1 + (1 + y - 1) * (width + 1)]);
	  di2 =
	    _mm_loadu_ps (&((float *) tmp1)
			  [1 + x - 1 + (1 + y) * (width + 1)]);
	  di3 =
	    _mm_loadu_ps (&((float *) tmp1)
			  [1 + x - 1 + (1 + y + 1) * (width + 1)]);

	  // scalar multiplication of Mli and di
	  mul1 = _mm_mul_ps (Mli1, di1);
	  mul2 = _mm_mul_ps (Mli2, di2);
	  mul3 = _mm_mul_ps (Mli3, di3);
	  // addition of the obatined scalar products
	  add = _mm_add_ps (mul1, mul2);
	  __m128 result1 = _mm_add_ps (add, mul3);

	  // scalar multiplication of Mwi and di
	  mul1 = _mm_mul_ps (Mwi1, di1);
	  mul2 = _mm_mul_ps (Mwi2, di2);
	  mul3 = _mm_mul_ps (Mwi3, di3);
	  // addition of the obatined scalar products
	  add = _mm_add_ps (mul1, mul2);
	  __m128 result2 = _mm_add_ps (add, mul3);

	  // calculate Ql(x,y)
	  float *res1 = (float *) &result1;
	  float Ql = res1[0] + res1[1] + res1[2];
      
	  // calculate Qw(x,y)
	  float *res2 = (float *) &result2;
	  float Qw = res2[0] + res2[1] + res2[2];
	  Qw /= (float) 16;

	  // extraction of Q(x,y) from tmp1
	  float d = ((float *) tmp1)[1 + x + (1 + y) * (width + 1)];

	  // calculate Q'(x,y)
	  Ql = Ql > 0 ? (Ql / (float) 1020) : ((-Ql) / (float) 1020);
	  if(x == 100 ) {
		//printf("x : %d -- y : %d -- Ql_simd : %f\n" ,x,y, Ql);
	  }
	  float Qprime = Ql * d + (1 - Ql) * Qw;

	  // fill result with the calculated Q' (x, y) 
	  result[y * width  + x]  = (uint8_t) (Qprime);
    }}
}