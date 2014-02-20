
/* [sinusoidality~] is a spectral analysis/resynthsis exern for Pure Data */
/* copyright David Medine */
/* released under the GPL */

#include "m_pd.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>


static t_class *sinusoidality_tilde_class;

typedef struct _sinusoidality_tilde{
  t_object x_obj;
  t_sample f;
  t_int count;
  t_int technique, smoothing;
  t_float sampleRate;
  t_int	bufferPosition, circBuffReadOut, circBuffStart, delayCount;
  long blockSize;
  t_word *x_vec;
  int arraypts;
  t_symbol *x_arrayname;
  t_float *sigInBuff, *noiseInBuff;
  t_float *inBuffer, *inWindowed, *inNoiseBuffer, *inNoiseWindowed;
  t_float *analWindow, *synthWindow; 
  t_float *spectra, *noiseSpectra, *sinesSpectra;
  t_float *outbuf, *noiseOutbuf, *sinesOutbuf;
  t_float *inShift, *inNoiseShift;
  t_float *FFTbin;
  t_float *buffer;
  t_float *realOut, *realOutNoise, *realOutSines;
  t_float *imagOut, *imagOutNoise, *imagOutSines, *rsqrt;
  t_float *nonoverlappedOutBuff, *nonoverlappedSinesOutBuff, *nonoverlappedNoiseOutBuff;
  double *noiseBinSmooth, *smoothNum;
  double *phaseSinReal, *phaseSinImag;
  double *phaseNoiseReal, *phaseNoiseImag;
  double *sinFullSpecReal, *sinFullSpecImag;
  double *noiseTimeSmooth;
  double *oldNoiseMagSpec, *newNoiseMagSpec;
  double *noiseFullSpecReal, *noiseFullSpecImag;
  double *pows, *mags, *noiseMags, *pows_minus_one, *pows_minus_two;
  double *sinCoefs, *scaledSinCoefs;
  double *sinMagSpec;
  double sfm, sfmNum, sfmDenom, zcScalar;
  long fftHalfSize, fftSize, log2n, zcCount;
  long inputTime, outputTime;
  long n,buffer_limit, overlap;
  t_float pi, twoPi;
  long circBuffLength;//size of the circular buffer(in samples)
  long N;//index along the circular buffer
  t_float *real, *noiseReal;
  t_float *imag, *noiseImag;
  t_float thresh, z;
  long processframes, framesLeft;
  int FFTtick, DSPtick;
  double startPoint, windowIndex;

  t_float ts; /*timestretch*/
  t_float ps; //pitchshift

}t_sinusoidality_tilde;

/*-----------------------------functiondeclaratins----------------------------*/
static void sinusoidality_tilde_fft(t_sinusoidality_tilde *x);
static void init_window(t_sinusoidality_tilde *x);
static void reinit(t_sinusoidality_tilde *x);
//------------------------------------------------------------------------//
static void sinusoidality_tilde_smoothing(t_sinusoidality_tilde *x)
{
  //this doesn't work at all!!!!!!!!!!


  long i, j, k;
  t_float a = 0.4; //this is the coefficient used in noise smoothing
  t_float one_minus_a = (1.0-a);
  int b = 10;//this guy is for spectral smoothing accross bins
  int invWindowSize = (1/(x->fftHalfSize + 1));

  //noise smoothing
  for(i=0; i<=x->fftHalfSize; i++)
    {
      //first smooth accross time
      x->noiseTimeSmooth[i] = ((a*x->oldNoiseMagSpec[i]) + ((one_minus_a) * x->newNoiseMagSpec[i]));
      
      //copy the new to the old
      x->oldNoiseMagSpec[i] = x->newNoiseMagSpec[i];
      
      //then smooth accross bins
      //loop to find average over b bins
      for(j = -b/2; j < b/2; j++)
	{
	  k = j + i;
	  if (k < 0)
	    k+=(x->fftHalfSize + 1);
	  if (k>x->fftHalfSize + 1)
	    k -= x->fftHalfSize+1;
	  x->smoothNum[i] += x->noiseTimeSmooth[i + k];
	}
  
  
      //find the average
      x->newNoiseMagSpec[i] = x->smoothNum[i]/b;
    }
}
/****************these methods are the sinusoidality d**********************/
static void sinusoidality_tilde_sigmund(t_sinusoidality_tilde *x)
{
}
static void sinusoidality_tilde_powpersist(t_sinusoidality_tilde *x)
{
 long i, j, k;
  t_float a = 0.4; //this is the coefficient used in noise smoothing
  t_float one_minus_a = (1.0-a);
  int b = 10;//this guy is for spectral smoothing accross bins
  int invWindowSize = (1/(x->fftHalfSize + 1));

  double maxpow;
  double sfm, sfmNum, sfmDenom, roots[x->fftHalfSize+1];

  //power spectrum defines sinusoidality, for more about technique see Apel[20]
  //find the sines
  for(i=0; i<=x->fftHalfSize; i++)
    {
      x->mags[i] = q8_sqrt((x->real[i] * x->real[i]) + (x->imag[i] * x->imag[i]));//find the magnitudes
      x->pows[i] = x->mags[i] * x->mags[i];//find the powers
      //find the noise magnitudes while we're at it, we will need them
      x->noiseMags[i] = q8_sqrt((x->noiseReal[i] * x->noiseReal[i]) + (x->noiseImag[i] * x->noiseImag[i]));
    }


  //first find the max value in each frame
  maxpow=x->pows[2];
  for(i=3; i<(x->fftHalfSize-1); i++)
    {
      if(maxpow < x->pows[i]/(.6*x->pows[i-2]*x->pows[i+2]))
	 maxpow =  x->pows[i]/(.6*x->pows[i-2]*x->pows[i+2]);
    }
  
  //then calculate the coefficients and copy the new into the old


  for(i=2; i<x->fftHalfSize-1; i++)	
      x->sinCoefs[i] = x->pows[i]/(.6*x->pows[i-2]*x->pows[i+2])/maxpow;


  for(i=0; i<=x->fftHalfSize; i++)
    {
      
      x->sinMagSpec[i] = x->mags[i] * x->sinCoefs[i];
      
      x->newNoiseMagSpec[i] = x->mags[i] - x->sinMagSpec[i];
      if(x->sinCoefs[i]>x->thresh)  x->newNoiseMagSpec[i] = 0;
      else x->sinMagSpec[i] = 0;
      
    }


}

static void sinusoidality_tilde_power(t_sinusoidality_tilde *x)
{

  
  long i, j, k;
  t_float a = 0.4; //this is the coefficient used in noise smoothing
  t_float one_minus_a = (1.0-a);
  int b = 10;//this guy is for spectral smoothing accross bins
  int invWindowSize = (1/(x->fftHalfSize + 1));

  double maxpow;
  double sfm, sfmNum, sfmDenom, roots[x->fftHalfSize+1];

  //power spectrum defines sinusoidality, for more about technique see Apel[20]
  //find the sines
  for(i=0; i<=x->fftHalfSize; i++)
    {
      x->mags[i] = q8_sqrt((x->real[i] * x->real[i]) + (x->imag[i] * x->imag[i]));//find the magnitudes
      x->pows[i] = x->mags[i] * x->mags[i];//find the powers
      //find the noise magnitudes while we're at it, we will need them
      x->noiseMags[i] = q8_sqrt((x->noiseReal[i] * x->noiseReal[i]) + (x->noiseImag[i] * x->noiseImag[i]));
    }
    
  //find the coefficients
  //assuming sinusoidal peaks are greater in power than noise peaks the normalized power spectrum
  //indicates sinusoidality

  //first find the max value in each frame
  maxpow=x->pows[0];
  for(i=1; i<=(x->fftHalfSize+1); i++)
    {
      if(maxpow < x->pows[i])
	maxpow = x->pows[i];
    }
  
  //then calculate the coefficients
  for(i=0; i<=x->fftHalfSize; i++)	
    x->sinCoefs[i] = x->pows[i]/maxpow;
   
  //now find the spectral flatness measure

  //fuck this shit!!!!!!!!!!!!!!!!!!!!:

  //sfm is the ratio of the geometric and arithmetic means of the magnitude spectrum fo this frame
 /*  sfmNum = 1.0; */
/*   sfmDenom = 0.0; */
/*   for(i = 0; i<=x->fftHalfSize; i++) */
/*     { */
/*       if(x->pows[i]>0.0)//if it is zero, we don't want it????not sure about this */
/* 	{ */
/* 	  roots[i] *= pow(x->pows[i]*1000000, invWindowSize); */
/* 	  sfmDenom += x->pows[i]*1000000;//scale these guys so that they are not so close to zero */
/* 	} */
/*     } */
/*   for(i = 0; i<=x->fftHalfSize; i++)sfmNum*=roots[i]; */

/*   x->sfmNum = sfmNum; */
/*   x->sfmDenom = (sfmDenom * invWindowSize)+.000000000000001; */
/*   x->sfm = sfm = sfmNum/sfmDenom; */
  
/*   //weight the coeficients with the sfm */
/*   for(i=0; i<=x->fftHalfSize; i++) */
/*     x->scaledSinCoefs[i] = x->sinCoefs[i] * sfm; */
  

//scale the coefs with zero crossing percentage for each frame

  for(i=0; i<=x->fftHalfSize; i++)
  x->sinCoefs[i] = pow(x->sinCoefs[i], x->zcScalar); 

  //now that we have the sinusoidal coefs., find the sin and noise spectra above threshhold
  if(x->x_arrayname)
    {
      for(i=0; i<=x->fftHalfSize; i++)
	{
	  x->thresh = x->x_vec[i].w_float;
	  x->sinMagSpec[i] = x->mags[i] * x->sinCoefs[i];
	  
	  x->newNoiseMagSpec[i] = x->mags[i] - x->sinMagSpec[i];
	  if(x->sinCoefs[i]>x->thresh)  x->newNoiseMagSpec[i] = 0;
	  else x->sinMagSpec[i] = 0;
	  
	}
    }
  

  else
    {
      for(i=0; i<=x->fftHalfSize; i++)
	{
	  
	  x->sinMagSpec[i] = x->mags[i] * x->sinCoefs[i];
	  
	  x->newNoiseMagSpec[i] = x->mags[i] - x->sinMagSpec[i];
	  if(x->sinCoefs[i]>x->thresh)  x->newNoiseMagSpec[i] = 0;
	  else x->sinMagSpec[i] = 0;
	  
	}
    }
  
  //here begins the routine to rescale the energy and recalculate the noise element
}


static void sinusoidality_tilde_fftGo(t_sinusoidality_tilde *x)
{
  
  long i, j, k;
  

  /***************************************************************/
  //buffer stuff in in preparation for fft:

  //pushout last block
  //set to zero for safety
  for(i=0; i<x->blockSize; i++)
    {
      x->sigInBuff[i] = 0;
      x->noiseInBuff[i] = 0;
    }
  
  //then shift over by one block
  for(i=0; i<(x->fftSize-x->blockSize); i++)
    {
      x->sigInBuff[i]=x->sigInBuff[x->blockSize+i];
      x->noiseInBuff[i]=x->noiseInBuff[x->blockSize+i];
    }

  //push in new block
  for(i=0; i<x->blockSize; i++)
    {
      x->sigInBuff[x->fftSize-x->blockSize+i]=x->inBuffer[i];
      x->noiseInBuff[x->fftSize-x->blockSize+i]=x->inNoiseBuffer[i];
    }
  for(i=0; i<x->fftSize; i++)
    {
      x->inWindowed[i] = x->sigInBuff[i] * x->analWindow[i];
      x->inNoiseWindowed[i] = x->noiseInBuff[i] * x->analWindow[i];
    }

  //find the zero crossings:
  int tmp, tmp_minus_one;
  
  x->zcCount = 0;
  tmp_minus_one = -1;
  if(x->sigInBuff[0]>= 0.0)
    tmp_minus_one = 1; 

  for(i=1; i<x->fftSize; i++)
    {  
      if(x->sigInBuff[i]>=0.0)
	tmp = 1;
      else tmp = -1;
      if(tmp!=tmp_minus_one)x->zcCount++;
      tmp_minus_one = tmp;
    }
  x->zcScalar = (double)x->zcCount/(double)x->fftSize;

  //do the FFT
  mayer_realfft(x->fftSize, x->inWindowed);
  mayer_realfft(x->fftSize, x->inNoiseWindowed);
  
  for(i=0; i<=x->fftHalfSize; i++)//  nyquist
    {
      x->real[i] = x->inWindowed[i];
      x->noiseReal[i] = x->inNoiseWindowed[i];
    }
  
  
  //imag routine
  x->imag[0]= x->noiseImag[0]= 0;  // 0 DC
  
  for(i=(x->fftSize-1), j=1; i>x->fftHalfSize; i--, j++)
    {
      x->imag[j] = x->inWindowed[i];
      x->noiseImag[j] = x->inNoiseWindowed[i];
    }
  x->imag[x->fftHalfSize]=x->noiseImag[x->fftHalfSize]=0;


  //if(x->technique == 1)
    sinusoidality_tilde_power(x);
  

  //first, normalize the spectrum to extract phase only info
  for(i=0; i<=x->fftHalfSize; i++)
    {
      x->phaseSinReal[i] = x->real[i]/(x->mags[i] +
				       .00000000000000001);
      x->phaseSinImag[i] = x->imag[i]/(x->mags[i] +
				       .00000000000000001);
      x->phaseNoiseReal[i] = x->noiseReal[i]/(x->noiseMags[i] +
					      .000000000000000000001);
      x->phaseNoiseImag[i] = x->noiseImag[i]/(x->noiseMags[i] +
					      .000000000000000000001);
      
      
      //then find the sinusoidal spectrum
      x->sinFullSpecReal[i] = x->phaseSinReal[i] * x->sinMagSpec[i];
      x->sinFullSpecImag[i] = x->phaseSinImag[i] * x->sinMagSpec[i];
    }
  if(x->smoothing)  
    sinusoidality_tilde_smoothing(x);
    
  //use the smoothed noise magnitude spectrum to find the amplitudes of the noise phase spectrum
  
  for(i=0; i<=x->fftHalfSize; i++)
    {
      x->noiseFullSpecReal[i] = x->phaseNoiseReal[i] * x->newNoiseMagSpec[i];
      x->noiseFullSpecImag[i] = x->phaseNoiseImag[i] * x->newNoiseMagSpec[i];
      
      
      //add the new noise with the scaled sines
      x->realOut[i] = (x->sinFullSpecReal[i] +  x->noiseFullSpecReal[i]);
      x->imagOut[i] = (x->sinFullSpecImag[i] +  x->noiseFullSpecImag[i]);
      
    }



  /*-----------------------------------------------*/
  //ifft and output prep methods:
  
  //copy real[0-512] into spectra [0-512]
  for(i=0; i<=x->fftHalfSize; i++)  // +1 to include Nyquist
    {
      x->spectra[i] = x->realOut[i];
      x->noiseSpectra[i] = x->noiseFullSpecReal[i];
      x->sinesSpectra[i] = x->sinFullSpecReal[i];
    }
  
  //copyimag[511-1] into spectra[513 - 1023]
  for(j = x->fftHalfSize -1, i = x->fftHalfSize + 1; i < x->fftSize; j--, i++)
    {
      x->spectra[i] = x->imagOut[j];
      x->noiseSpectra[i] = x->noiseFullSpecImag[j];
      x->sinesSpectra[i] = x->sinFullSpecImag[j];
    }
  
  mayer_realifft(x->fftSize, x->spectra);
  mayer_realifft(x->fftSize, x->noiseSpectra);
  mayer_realifft(x->fftSize, x->sinesSpectra);
  
  //WB routine for overlap add:
  //first, multiply by the synth window
  for(i=0; i<x->fftSize; i++)
    {
      x->spectra[i] *= x->synthWindow[i];
      x->noiseSpectra[i] *= x->synthWindow[i];
      x->sinesSpectra[i] *= x->synthWindow[i];
    }
  
  //  first shift output in nonoverlapped buffer
  for(i=0; i<((x->overlap-1)*x->fftSize); i++)
    {
      x->nonoverlappedOutBuff[i] = x->nonoverlappedOutBuff[x->fftSize+i];
      x->nonoverlappedNoiseOutBuff[i] = x->nonoverlappedNoiseOutBuff[x->fftSize+i];
      x->nonoverlappedSinesOutBuff[i] = x->nonoverlappedSinesOutBuff[x->fftSize+i];
    }
  
  //  then, write a new window in
  for(i=0; i<x->fftSize; i++)
    {
      x->nonoverlappedOutBuff[((x->overlap-1)*x->fftSize)+i] = x->spectra[i];
      x->nonoverlappedNoiseOutBuff[((x->overlap-1)*x->fftSize)+i] = x->noiseSpectra[i];
      x->nonoverlappedSinesOutBuff[((x->overlap-1)*x->fftSize)+i] = x->sinesSpectra[i];
    }
  
  // init this chunk of the final output so it can be summed in the for() below
  for(i=0; i<(x->blockSize); i++)
    {
      x->outbuf[i] = 0.0;
      x->noiseOutbuf[i] = 0.0;
      x->sinesOutbuf[i] = 0.0;
    }
  
  
  // do the overlap/add
  for(i=0; i<x->overlap; i++)
    for(j=0; j<(x->blockSize); j++)
      {
	x->outbuf[j] += x->nonoverlappedOutBuff[(i*x->fftSize)+((x->overlap-i-1)*x->blockSize)+j/*overlap_chunk*/];
	x->noiseOutbuf[j] += x->nonoverlappedNoiseOutBuff[(i*x->fftSize)+((x->overlap-i-1)*x->blockSize)+j/*overlap_chunk*/];
	x->sinesOutbuf[j] += x->nonoverlappedSinesOutBuff[(i*x->fftSize)+((x->overlap-i-1)*x->blockSize)+j/*overlap_chunk*/];
      }
}



/*-------------------------------process samples------------------------------*/

t_int *sinusoidality_tilde_perform (t_int *w)
{
  t_sinusoidality_tilde *x = (t_sinusoidality_tilde *)(w[1]);
  t_sample  *in = (t_sample *)(w[2]);
  t_sample *noiseIn = (t_sample *)(w[3]);
  t_sample *out1 = (t_sample *)(w[4]);
  t_sample *out2 = (t_sample *)(w[5]);
  t_sample *out3 = (t_sample *)(w[6]);
  int n = (int)(w[7]);
   
  long i;


  long frames;
  long framesLeft, processframes;

  //this is code copied from Tom Erbe:

  framesLeft = n;	
  while ( framesLeft > 0 )
    {
     
      if (framesLeft + x->bufferPosition < x->blockSize)
	processframes = framesLeft;
      else
	processframes = x->blockSize - x->bufferPosition;
     
      memcpy(x->inBuffer+(x->bufferPosition), in, processframes * sizeof(t_float));
      memcpy(x->inNoiseBuffer+(x->bufferPosition), noiseIn, processframes * sizeof(t_float));
      
      for (i=0; i<processframes; i++)
	{
	  *out1++ = x->outbuf[i+x->bufferPosition];
	  *out2++ = x->noiseOutbuf[i+x->bufferPosition];
	  *out3++ = x->sinesOutbuf[i+x->bufferPosition];
	}
      // increment in and out pointers	
    /*   out1 += processframes; */
/*       out2 += processframes; */
/*       out3 += processframes; */
      in += processframes;
      noiseIn += processframes;

      // increment bufferPostion, if the bufferPosition hits the blockSize (1/4 of FFT size)
      // perform another FFT.
      x->bufferPosition += processframes;
      if (x->bufferPosition >= x->blockSize)
	{
	  x->bufferPosition = 0;
	  sinusoidality_tilde_fftGo(x);
	}
      // decrement framesLeft by the number of frames (samples) processed
      framesLeft -= processframes;
    }

  return(w+8);
}




/*-----------------------------dsp-----------------------------------*/
void sinusoidality_tilde_dsp(t_sinusoidality_tilde *x, t_signal **sp)
{
  x->sampleRate = sp[0]->s_sr;
  dsp_add(sinusoidality_tilde_perform, 7, x, sp[0]->s_vec, sp[1]->s_vec, sp[2]->s_vec, sp[3]->s_vec, sp[4]->s_vec, sp[0]->s_n);
  // FFTbinTable(x);
}

/*--------------------thresh-----------------------------*/
void sinusoidality_tilde_thresh(t_sinusoidality_tilde *x, t_floatarg g)
{
  if(g<0.0) x->thresh = 0.0;
  else if(g>=1.0)x->thresh = 1.0;
  else x->thresh = g;
}

/*--------------------scaleSFM-----------------------------*/
void sinusoidality_tilde_scaleSFM(t_sinusoidality_tilde *x, t_floatarg g)
{
  x->z = g * 100;
  
}

/*----------------------table for thresh--------------*/
static void sinusoidality_tilde_table(t_sinusoidality_tilde *x, t_symbol *s, int argc, t_atom *argv)
{
  t_garray *a;
  int vecsize = argc-1;
  int i;

  x->x_arrayname = atom_getsymbol(argv);

  if(!(a = (t_garray *)pd_findbyclass(x->x_arrayname, garray_class)))
    pd_error(x, "%s: no such array", s->s_name);
 else if(!garray_getfloatwords(a, &x->arraypts, &x->x_vec))
    	pd_error(x, "%s: bad template for table", x->x_arrayname->s_name);
 /* else if(vecsize!=x->fftSize) */
 /*   pd_error(x, "%s: array must have %d points", x->x_arrayname->s_name, x->fftSize);   */
else
    {
      for(i=1; i<argc; i++)
	x->x_vec[i-1].w_float = atom_getfloat(argv+i);
      
      garray_redraw(a);
    }
}

/*--------------------post-----------------------------*/
void sinusoidality_tilde_post(t_sinusoidality_tilde *x)
{

  post("%f", x->zcScalar);
/*   int i; */
/*   for(i=0; i<=200; i++) */
/*     post("%d: %f", i, x->pows[i]); */

}


void sinusoidality_tilde_smooth(t_sinusoidality_tilde *x, t_floatarg g)
{
  if(g)
    x->smoothing = 1;
  else
    x->smoothing = 0;
}

/*----------------------------sinusoidality_new--------------------------------*/
void *sinusoidality_tilde_new(t_symbol *s, int argc, t_atom *argv)
{
  int i;
  t_sinusoidality_tilde *x = (t_sinusoidality_tilde *)pd_new(sinusoidality_tilde_class);

  x->technique = 1;
  
  while (argc > 0)
    {
      t_symbol *firstarg = atom_getsymbolarg(0, argc, argv);
      if (!strcmp(firstarg->s_name, "-amp"))
        {
	  t_symbol*secondarg = atom_getsymbolarg(1, argc, argv);
	  if(!strcmp(secondarg->s_name, "power"))
	    {
	      x->technique = 1;
  	      post("sinusoidality measure based on power");
	    }
	  else if(!strcmp(secondarg->s_name, "pow-persist"))
	    {
	      x->technique = 2;
	      post("sinusoidality measure based on power peristence");
	    }
	  else if(!strcmp(secondarg->s_name, "sigmund"))
	    {
	      x->technique = 3;
	      post("sinusoidality measure based on  sigmund");
	    }


	  else
	    {
	      pd_error(x, "sinusoidality: %s: unknown flag or argument missing",
		       firstarg->s_name);
	      
	    }
	  argc-=2, argv+=2;
	}
      else if(!strcmp(firstarg->s_name, "-phase"))
	{
	  t_symbol*secondarg = atom_getsymbolarg(1, argc, argv);
	  if(!strcmp(secondarg->s_name, "charpentier"))
	    {
	      x->technique = 4;
	      post("sinusoidality measure based on phase, charpentier");
	    }
	  argc-=2, argv+=2;
	}
      else
	{
	  pd_error(x, "sinusoidality: %s: unknown flag or argument missing",
		   firstarg->s_name);
	  argc--, argv++;
	}
    }

    
  inlet_new(&x->x_obj, &x->x_obj.ob_pd, &s_signal, &s_signal);
  inlet_new(&x->x_obj, &x->x_obj.ob_pd, gensym("float"), gensym("thresh"));

  outlet_new(&x->x_obj, &s_signal);
  outlet_new(&x->x_obj, &s_signal);
  outlet_new(&x->x_obj, &s_signal);
  x->log2n = (long) log2f ( (float) x->fftSize);
  x->bufferPosition = x->circBuffReadOut = x->circBuffStart = 0;
  x->inputTime = x->outputTime = 0;
  x->fftSize = 1024;
  x->fftHalfSize = x->fftSize >> 1;
  x->blockSize = x->fftSize >> 2;
  x->pi = 4.0f * atanf(1.0f);
  x->twoPi = 8.0f * atanf(1.0f);
  x->count = 0;
  x->ts = x->ps = 1.0;   
  x->startPoint = x->windowIndex = x->delayCount = x->outputTime = 0;
  x->FFTtick = 0;
  x->circBuffLength = x->fftSize;//for now, circBuf is 60 seconds long @ 48k samplerate
  x->N = 0;//intialize circBuf index to zero
  x->processframes = x->framesLeft = 0;
  x->DSPtick = 0;   
  x->n = 64;
  x->overlap = x->fftSize/x->blockSize;
  x->buffer_limit = x->fftSize/x->n/x->overlap;
  x->thresh = 10;
  x->z = 10000;
  x->sfm = 0;
  x->smoothing= 0;
  x->zcScalar = 0.0;
  x->zcCount = 0;


	
    //set up memory
    //
    x->sigInBuff = (t_float *)t_getbytes(0);
    x->sigInBuff = (t_float *) t_resizebytes(x->sigInBuff, 0, sizeof(t_float) * x->fftSize);
    //
    x->noiseInBuff = (t_float *)t_getbytes(0);
    x->noiseInBuff = (t_float *) t_resizebytes(x->noiseInBuff, 0, sizeof(t_float) * x->fftSize);	
    //
    x->inBuffer = (t_float *)t_getbytes(0);
    x->inBuffer = (t_float *) t_resizebytes(x->inBuffer, 0, sizeof(t_float) * x->fftSize);
    //
    x->inNoiseBuffer = (t_float *)t_getbytes(0);
    x->inNoiseBuffer = (t_float *) t_resizebytes(x->inNoiseBuffer, 0, sizeof(t_float) * x->fftSize);    
    //
    x->inWindowed = (t_float *)t_getbytes(0);
    x->inWindowed = (t_float *) t_resizebytes(x->inWindowed, 0, sizeof(t_float) * x->fftSize);
    //
    x->inNoiseWindowed = (t_float *)t_getbytes(0);
    x->inNoiseWindowed = (t_float *) t_resizebytes(x->inNoiseWindowed, 0, sizeof(t_float) * x->fftSize);
    //
    x->inShift = (t_float *)t_getbytes(0);
    x->inShift = (t_float *) t_resizebytes(x->inShift, 0, sizeof(t_float) * x->fftSize);
    //
    x->inNoiseShift = (t_float *)t_getbytes(0);
    x->inNoiseShift = (t_float *) t_resizebytes(x->inNoiseShift, 0, sizeof(t_float) * x->fftSize);
    //
    x->analWindow = (t_float *)t_getbytes(0);
    x->analWindow = (t_float *) t_resizebytes(x->analWindow, 0, sizeof(t_float) * x->fftSize);
    //
    x->synthWindow = (t_float *)t_getbytes(0);
    x->synthWindow = (t_float *) t_resizebytes(x->synthWindow, 0, sizeof(t_float) * x->fftSize);
    //	
    x->real = (t_float *)t_getbytes(0);
    x->real = (t_float *) t_resizebytes(x->real, 0, sizeof(t_float) * (x->fftHalfSize +1));
    //
    x->imag = (t_float *)t_getbytes(0);
    x->imag = (t_float *) t_resizebytes(x->imag, 0, sizeof(t_float) * (x->fftHalfSize+1));
    //	
    x->noiseReal = (t_float *)t_getbytes(0);
    x->noiseReal = (t_float *) t_resizebytes(x->noiseReal, 0, sizeof(t_float) * (x->fftHalfSize +1));
    //
    x->noiseImag = (t_float *)t_getbytes(0);
    x->noiseImag = (t_float *) t_resizebytes(x->noiseImag, 0, sizeof(t_float) * (x->fftHalfSize+1));
    //
    x->spectra = (t_float *)t_getbytes(0);
    x->spectra = (t_float *) t_resizebytes(x->spectra, 0, sizeof(t_float) * x->fftSize);
    //
    x->noiseSpectra = (t_float *)t_getbytes(0);
    x->noiseSpectra = (t_float *) t_resizebytes(x->noiseSpectra, 0, sizeof(t_float) * x->fftSize);
    //
   x->sinesSpectra = (t_float *)t_getbytes(0);
    x->sinesSpectra = (t_float *) t_resizebytes(x->sinesSpectra, 0, sizeof(t_float) * x->fftSize);
    //    
    x->outbuf = (t_float *)t_getbytes(0);
    x->outbuf = (t_float *) t_resizebytes(x->outbuf, 0, sizeof(t_float) * x->fftSize);
    //
    x->noiseOutbuf = (t_float *)t_getbytes(0);
    x->noiseOutbuf = (t_float *) t_resizebytes(x->noiseOutbuf, 0, sizeof(t_float) * x->fftSize);
   //
    x->sinesOutbuf = (t_float *)t_getbytes(0);
    x->sinesOutbuf = (t_float *) t_resizebytes(x->sinesOutbuf, 0, sizeof(t_float) * x->fftSize);
    //
    x->FFTbin = (t_float *)t_getbytes(0);
    x->FFTbin = (t_float *) t_resizebytes(x->FFTbin, 0, sizeof(t_float) * x->fftSize);
    //
    x->buffer = (t_float *)t_getbytes(0);
    x->buffer = (t_float *) t_resizebytes(x->buffer, 0, sizeof(t_float) * x->fftSize);
    //
    x->realOut = (t_float *)t_getbytes(0);
    x->realOut = (t_float *)t_resizebytes(x->realOut, 0, sizeof(t_float) * x->fftHalfSize+1);
    //
    x->imagOut = (t_float *)t_getbytes(0);
    x->imagOut = (t_float *)t_resizebytes(x->imagOut, 0, sizeof(t_float) * x->fftHalfSize+1);
    //
    x->realOutNoise = (t_float *)t_getbytes(0);
    x->realOutNoise = (t_float *)t_resizebytes(x->realOutNoise, 0, sizeof(t_float) * x->fftHalfSize+1);
    //
    x->imagOutNoise = (t_float *)t_getbytes(0);
    x->imagOutNoise = (t_float *)t_resizebytes(x->imagOutNoise, 0, sizeof(t_float) * x->fftHalfSize+1);
    //
    x->rsqrt = (t_float *)t_getbytes(0);
    x->rsqrt = (t_float *)t_resizebytes(x->rsqrt, 0, sizeof(t_float) * x->fftHalfSize+1);
    //
    x->noiseReal = (t_float *)t_getbytes(0);
    x->noiseReal = (t_float *)t_resizebytes(x->noiseReal, 0, sizeof(t_float) * x->fftHalfSize+1);
    //
    x->noiseImag = (t_float *)t_getbytes(0);
    x->noiseImag = (t_float *)t_resizebytes(x->noiseImag, 0, sizeof(t_float) * x->fftHalfSize+1);
    //
    x->phaseSinReal = (double *)t_getbytes(0);
    x->phaseSinReal = (double *)t_resizebytes(x->phaseSinReal, 0, sizeof(double) * x->fftHalfSize+1);
    //
    x->phaseSinImag = (double *)t_getbytes(0);
    x->phaseSinImag = (double *)t_resizebytes(x->phaseSinImag, 0, sizeof(double) * x->fftHalfSize+1);
    //
    x->phaseNoiseReal = (double *)t_getbytes(0);
    x->phaseNoiseReal = (double *)t_resizebytes(x->phaseNoiseReal, 0, sizeof(double) * x->fftHalfSize+1);
    //
    x->phaseNoiseImag = (double *)t_getbytes(0);
    x->phaseNoiseImag = (double *)t_resizebytes(x->phaseNoiseImag, 0, sizeof(double) * x->fftHalfSize+1);
    //
    x->sinFullSpecReal = (double *)t_getbytes(0);
    x->sinFullSpecReal = (double *)t_resizebytes(x->sinFullSpecReal, 0, sizeof(double) * x->fftHalfSize+1);
    //
    x->sinFullSpecImag = (double *)t_getbytes(0);
    x->sinFullSpecImag = (double *)t_resizebytes(x->sinFullSpecImag, 0, sizeof(double) * x->fftHalfSize+1);
    //
    x->noiseTimeSmooth = (double *)t_getbytes(0);
    x->noiseTimeSmooth = (double *)t_resizebytes(x->noiseTimeSmooth, 0, sizeof(double) * x->fftHalfSize+1);
    //
    x->noiseBinSmooth = (double *)t_getbytes(0);
    x->noiseBinSmooth = (double *)t_resizebytes(x->noiseBinSmooth, 0, sizeof(double) * x->fftHalfSize+1);
    //
    x->smoothNum = (double *)t_getbytes(0);
    x->smoothNum = (double *)t_resizebytes(x->smoothNum, 0, sizeof(double) * x->fftHalfSize+1);
    //
    x->oldNoiseMagSpec = (double *)t_getbytes(0);
    x->oldNoiseMagSpec = (double *)t_resizebytes(x->oldNoiseMagSpec, 0, sizeof(double) * x->fftHalfSize+1);
    //
    x->newNoiseMagSpec = (double *)t_getbytes(0);
    x->newNoiseMagSpec = (double *)t_resizebytes(x->newNoiseMagSpec, 0, sizeof(double) * x->fftHalfSize+1);
    //
    x->noiseFullSpecReal = (double *)t_getbytes(0);
    x->noiseFullSpecReal = (double *)t_resizebytes(x->noiseFullSpecReal, 0, sizeof(double) * x->fftHalfSize+1);
    //
    x->noiseFullSpecImag = (double *)t_getbytes(0);
    x->noiseFullSpecImag = (double *)t_resizebytes(x->noiseFullSpecImag, 0, sizeof(double) * x->fftHalfSize+1);
    //
    x->pows = (double *)t_getbytes(0);
    x->pows = (double *)t_resizebytes(x->pows, 0, sizeof(double) * x->fftHalfSize+1);
    //
    x->mags = (double *)t_getbytes(0);
    x->mags = (double *)t_resizebytes(x->mags, 0, sizeof(double) * x->fftHalfSize+1);
    //
    x->noiseMags = (double *)t_getbytes(0);
    x->noiseMags = (double *)t_resizebytes(x->noiseMags, 0, sizeof(double) * x->fftHalfSize+1);
    //
    x->sinCoefs = (double *)t_getbytes(0);
    x->sinCoefs = (double *)t_resizebytes(x->sinCoefs, 0, sizeof(double) * x->fftHalfSize+1);
    //
    x->scaledSinCoefs = (double *)t_getbytes(0);
    x->scaledSinCoefs = (double *)t_resizebytes(x->scaledSinCoefs, 0, sizeof(double) * x->fftHalfSize+1);
    //
    x->sinMagSpec = (double *)t_getbytes(0);
    x->sinMagSpec = (double *)t_resizebytes(x->sinMagSpec, 0, sizeof(double) * x->fftHalfSize+1);
    //
    x->nonoverlappedOutBuff = (t_float *)t_getbytes(0);
    x->nonoverlappedOutBuff = (t_float *)t_resizebytes(x->nonoverlappedOutBuff, 0, sizeof(t_float) * x->fftSize * x->overlap);
    //
    x->nonoverlappedNoiseOutBuff = (t_float *)t_getbytes(0);
    x->nonoverlappedNoiseOutBuff = (t_float *)t_resizebytes(x->nonoverlappedNoiseOutBuff, 0, sizeof(t_float) * x->fftSize * x->overlap);
    //
    x->nonoverlappedSinesOutBuff = (t_float *)t_getbytes(0);
    x->nonoverlappedSinesOutBuff = (t_float *)t_resizebytes(x->nonoverlappedSinesOutBuff, 0, sizeof(t_float) * x->fftSize * x->overlap);
   
    if(x->technique = 2)
    x->pows_minus_one = (double *)t_getbytes(0);
    x->pows_minus_one = (double *)t_resizebytes(x->pows_minus_one, 0, sizeof(double) * x->fftHalfSize+1);
    //
    x->pows_minus_two = (double *)t_getbytes(0);
    x->pows_minus_two = (double *)t_resizebytes(x->pows_minus_two, 0, sizeof(double) * x->fftHalfSize+1);


  
    //call windowing inits
    init_window(x);
    
  
	 return (void *)x;
}


	/***********init window*********/
static void init_window(t_sinusoidality_tilde *x)
{
  int i;
  for ( i = 0; i < x->fftSize; i++ )
    {
    //hann
     x->analWindow[i] = x->synthWindow[i] = 0.5f * (1 - cosf(x->twoPi*i/(x->fftSize-1)));
     x->synthWindow[i] *= 2.0/3.0/x->fftSize;
    }
     //hamm
  //x->analWindow[i] = x->synthWindow[i] = (float) (.54f - (.46f * cosf(x->twoPi * i / (x->fftSize - 1)) ) );
    
 
}

/*---------------------------reinit------------------------------------*/


/****************************destructor*********************************/
static void sinusoidality_tilde_free(t_sinusoidality_tilde *x)
{
  if(x->sigInBuff != 0) t_freebytes(x->sigInBuff, x->fftSize * sizeof(t_float));
  if(x->noiseInBuff != 0) t_freebytes(x->noiseInBuff, x->fftSize * sizeof(t_float));
  if(x->inBuffer != 0) t_freebytes(x->inBuffer, x->fftSize * sizeof(t_float));
  if(x->inNoiseBuffer != 0) t_freebytes(x->inBuffer, x->fftSize * sizeof(t_float));
  if(x->inWindowed != 0) t_freebytes(x->inWindowed, x->fftSize * sizeof(t_float));
  if(x->inNoiseWindowed != 0) t_freebytes(x->inNoiseWindowed, x->fftSize * sizeof(t_float));
  if(x->inShift != 0) t_freebytes(x->inShift, x->fftSize * sizeof(t_float));
  if(x->inNoiseShift != 0) t_freebytes(x->inNoiseShift, x->fftSize * sizeof(t_float));
  if(x->analWindow != 0) t_freebytes(x->analWindow, x->fftSize * sizeof(t_float));
  if(x->synthWindow != 0) t_freebytes(x->synthWindow, x->fftSize * sizeof(t_float));
  if(x->outbuf != 0) t_freebytes(x->outbuf, x->fftSize * sizeof(t_float));
  if(x->noiseOutbuf != 0) t_freebytes(x->noiseOutbuf, x->fftSize * sizeof(t_float));
  if(x->sinesOutbuf != 0) t_freebytes(x->sinesOutbuf, x->fftSize * sizeof(t_float));
  if(x->spectra != 0) t_freebytes(x->spectra, x->fftSize * sizeof(t_float));
  if(x->noiseSpectra != 0) t_freebytes(x->noiseSpectra, x->fftSize * sizeof(t_float));
  if(x->sinesSpectra != 0) t_freebytes(x->sinesSpectra, x->fftSize * sizeof(t_float));
  if(x->real != 0) t_freebytes(x->real, x->fftSize * sizeof(t_float));
  if(x->imag != 0) t_freebytes(x->imag, x->fftSize * sizeof(t_float));
  if(x->noiseReal != 0) t_freebytes(x->noiseReal, x->fftSize * sizeof(t_float));
  if(x->noiseImag != 0) t_freebytes(x->noiseImag, x->fftSize * sizeof(t_float));
  if(x->FFTbin != 0) t_freebytes(x->FFTbin, x->fftSize * sizeof(t_float));
  if(x->buffer != 0) t_freebytes(x->buffer, x->fftSize * sizeof(t_float));
  if(x->realOut != 0) t_freebytes(x->realOut, sizeof(t_float) * x->fftHalfSize+1);
  if(x->imagOut != 0) t_freebytes(x->imagOut, sizeof(t_float) * x->fftHalfSize+1);
  if(x->realOutNoise != 0) t_freebytes(x->realOutNoise, sizeof(t_float) * x->fftHalfSize+1);
  if(x->imagOutNoise != 0) t_freebytes(x->imagOutNoise, sizeof(t_float) * x->fftHalfSize+1);
  if(x->rsqrt != 0) t_freebytes(x->rsqrt, sizeof(t_float) * x->fftHalfSize+1);
  if(x->noiseReal != 0) t_freebytes(x->noiseReal, sizeof(t_float) * x->fftHalfSize+1);
  if(x->noiseImag != 0) t_freebytes(x->noiseImag, sizeof(t_float) * x->fftHalfSize+1);
  if(x->phaseSinReal != 0) t_freebytes(x->phaseSinReal, sizeof(double) * x->fftHalfSize+1);
  if(x->phaseSinImag != 0) t_freebytes(x->phaseSinImag, sizeof(double) * x->fftHalfSize+1);
  if(x->phaseNoiseReal != 0) t_freebytes(x->phaseNoiseReal, sizeof(double) * x->fftHalfSize+1);
  if(x->phaseNoiseImag != 0) t_freebytes(x->phaseNoiseImag, sizeof(double) * x->fftHalfSize+1);
  if(x->sinFullSpecReal != 0) t_freebytes(x->sinFullSpecReal, sizeof(double) * x->fftHalfSize+1);
  if(x->sinFullSpecImag != 0) t_freebytes(x->sinFullSpecImag, sizeof(double) * x->fftHalfSize+1); 
  if(x->noiseTimeSmooth != 0) t_freebytes(x->noiseTimeSmooth, sizeof(double) * x->fftHalfSize+1); 
  if(x->noiseBinSmooth != 0) t_freebytes(x->noiseBinSmooth, sizeof(double) * x->fftHalfSize+1); 
  if(x->smoothNum != 0) t_freebytes(x->smoothNum, sizeof(double) * x->fftHalfSize+1); 
  if(x->oldNoiseMagSpec != 0) t_freebytes(x->oldNoiseMagSpec, sizeof(double) * x->fftHalfSize+1); 
  if(x->newNoiseMagSpec != 0) t_freebytes(x->newNoiseMagSpec, sizeof(double) * x->fftHalfSize+1); 
  if(x->noiseFullSpecReal != 0) t_freebytes(x->noiseFullSpecReal, sizeof(double) * x->fftHalfSize+1); 
  if(x->noiseFullSpecImag != 0) t_freebytes(x->noiseFullSpecImag, sizeof(double) * x->fftHalfSize+1); 
  if(x->pows != 0) t_freebytes(x->pows, sizeof(double) * x->fftHalfSize+1); 
  if(x->pows_minus_one != 0) t_freebytes(x->pows_minus_one, sizeof(double) * x->fftHalfSize+1); 
  if(x->pows_minus_two != 0) t_freebytes(x->pows_minus_two, sizeof(double) * x->fftHalfSize+1); 
  if(x->mags != 0) t_freebytes(x->mags, sizeof(double) * x->fftHalfSize+1); 
  if(x->noiseMags != 0) t_freebytes(x->noiseMags, sizeof(double) * x->fftHalfSize+1); 
  if(x->sinCoefs != 0) t_freebytes(x->sinCoefs, sizeof(double) * x->fftHalfSize+1); 
  if(x->scaledSinCoefs != 0) t_freebytes(x->scaledSinCoefs, sizeof(double) * x->fftHalfSize+1); 
  if(x->sinMagSpec != 0) t_freebytes(x->sinMagSpec, sizeof(double) * x->fftHalfSize+1); 
  if(x->nonoverlappedOutBuff != 0) t_freebytes(x->nonoverlappedOutBuff, sizeof(t_float) * x->fftSize * x->overlap);
  if(x->nonoverlappedNoiseOutBuff != 0) t_freebytes(x->nonoverlappedNoiseOutBuff, sizeof(t_float) * x->fftSize * x->overlap);
  if(x->nonoverlappedSinesOutBuff != 0) t_freebytes(x->nonoverlappedSinesOutBuff, sizeof(t_float) * x->fftSize * x->overlap);


}

/**************************set up routine*********************************/
void sinusoidality_tilde_setup(void)
{
  sinusoidality_tilde_class = class_new(gensym("sinusoidality~"),
					(t_newmethod)sinusoidality_tilde_new,
					(t_method)sinusoidality_tilde_free,
					sizeof(t_sinusoidality_tilde),
					0,
					A_GIMME,
					0);


  CLASS_MAINSIGNALIN(sinusoidality_tilde_class, t_sinusoidality_tilde, f);

  class_addmethod(sinusoidality_tilde_class, 
		  (t_method)sinusoidality_tilde_dsp, 
		  gensym("dsp"), 
		  0);

  class_addmethod(sinusoidality_tilde_class, 
		  (t_method)sinusoidality_tilde_post, 
		  gensym("post"), 
		  A_DEFFLOAT, 
		  0);

  class_addmethod(sinusoidality_tilde_class, 
		  (t_method)sinusoidality_tilde_thresh, 
		  gensym("thresh"), 
		  A_DEFFLOAT, 
		  0);

 class_addmethod(sinusoidality_tilde_class, 
		  (t_method)sinusoidality_tilde_table, 
		  gensym("table"), 
		  A_GIMME, 
		  0);

  /* class_addmethod(sinusoidality_tilde_class,  */
/* 		  (t_method)sinusoidality_tilde_scaleSFM,  */
/* 		  gensym("scaleSFM"),  */
/* 		  A_DEFFLOAT,  */
/* 		  0); */

 class_addmethod(sinusoidality_tilde_class, 
		  (t_method)sinusoidality_tilde_smooth, 
		  gensym("smooth"), 
		  A_DEFFLOAT, 
		  0);
}
  ////////////////////////////////////////////////////////
