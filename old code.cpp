

#include <fstream>
#include <iostream>

#include <cmath>
#include <ctime>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

using namespace std;

#define SWAP(a,b) tempr=(a);(a)=(b);(b)=tempr
double pi=4*atan((double)1);

#define		WINDOW(j,a,b) (1.0-fabs((((j)-1)-(a))*(b)))		//Bartlett
//#define		WINDOW(j,a,b) 1.0							//Square
//#define		WINDOW(j,a,b) (1.0-SQR((((j)-1)-(a))*(b)))	//Welch

void generate_signal(float data[], float fftdata[], unsigned long captured_samples, unsigned long sample_rate, double frequency )	{

	//example of a sin signal
	double amplitude = 2;
	for(unsigned long i=0; i<(captured_samples); i++)	{
		double temp=(double)(2*pi*frequency*((float)i/sample_rate));		
		data[i]=amplitude*sin(temp);	
		fftdata[i]=data[i];
		
	}

} 

// like abs() but for doubles. necessary because
// abs() only accepts and returns ints
double flabs(double f)	{
	return f * (f<0? -1: 1);
}

// returns the maximum value at a given bit-depth
int depthMax(int bits)	{
	return pow(2.0, bits-1);
}

// returns maximum value in an array of doubles
double getMax(float data[], long captured_samples)	{
double c, m;
int i;
for (i = 0, m = 0; i < captured_samples; ++i) {
	c = flabs(data[i]);
	if (c > m)
		m = c;
	}
return m;
}

void quantisation(float data[], long captured_samples,int dataMax, int depthMax  )	{

	for(int i=0; i<captured_samples; i++)	{			
			data[i]= floor(depthMax*(data[i]/dataMax));	
	}

}

void four1(float data[], unsigned long nn, int isign)
/*******************************************************************************
Replaces data[1..2*nn] by its discrete Fourier transform, if isign is input as
1; or replaces data[1..2*nn] by nn times its inverse discrete Fourier transform,
if isign is input as -1.  data is a complex array of length nn or, equivalently,
a real array of length 2*nn.  nn MUST be an integer power of 2 (this is not
checked for!).
*******************************************************************************/
{
	unsigned long n,mmax,m,j,istep,i;
	double wtemp,wr,wpr,wpi,wi,theta;
	float tempr,tempi;

	n=nn << 1;
	j=1;
	for (i=1;i<n;i+=2) {	/* This is the bit-reversal section of the routine. */
		if (j > i) {
			SWAP(data[j],data[i]);	/* Exchange the two complex numbers. */
			SWAP(data[j+1],data[i+1]);
		}
		m=n >> 1;
		while (m >= 2 && j > m) {
			j -= m;
			m >>= 1;
		}
		j += m;
	}
	mmax=2;
	while (n > mmax) {	/* Outer loop executed log2 nn times. */
		istep=mmax << 1;
		theta=isign*(6.28318530717959/mmax);	/* Initialize the trigonometric recurrence. */
		wtemp=sin(0.5*theta);
		wpr = -2.0*wtemp*wtemp;
		wpi=sin(theta);
		wr=1.0;
		wi=0.0;
		for (m=1;m<mmax;m+=2) {		/* Here are the two nested inner loops. */
			for (i=m;i<=n;i+=istep) {
				j=i+mmax;				/* This is the Danielson-Lanczos formula. */
				tempr=wr*data[j]-wi*data[j+1];
				tempi=wr*data[j+1]+wi*data[j];
				data[j]=data[i]-tempr;
				data[j+1]=data[i+1]-tempi;
				data[i] += tempr;
				data[i+1] += tempi;
			}
			wr=(wtemp=wr)*wpr-wi*wpi+wr;		/* Trigonometric recurrence. */
			wi=wi*wpr+wtemp*wpi+wi;
		}
		mmax=istep;
	}
}

void output_data (float data[],float fftdata[], unsigned long captured_samples, unsigned long sample_rate)	{
	
	ofstream outputfile("output.dat");
	unsigned long nn = captured_samples/2;

	outputfile << data[0] << " " <<  0.0 << " " << fftdata[0] << " " << 0.0 << " " << fabs(data[0]) << endl;
	for (unsigned long i = 1; i < nn; i++)	{
		outputfile << data[i]  << "	" << " " << (float)i/(float)captured_samples*(float)sample_rate << " ";
		outputfile <<  fftdata[2*i] << " " << fftdata[2*i + 1 ] << " ";
		outputfile << sqrt( fftdata[2*i]* fftdata[2*i] +  fftdata[2*i+1]* fftdata[2*i+1]) << endl;
	}
	outputfile << data[nn] << " " << 0.5*sample_rate << " " << fftdata[1] << " " << 0.0 << " " << fabs(fftdata[1]) << endl;

	outputfile.close();


	ofstream fftoutput("indata.dat");
	for (int i=0; i<captured_samples; i++)	{
		fftoutput<< data[i] << endl;
	}
	fftoutput.close();
	
}

void realft(float data[], unsigned long n, int isign)
{
	void four1(float data[], unsigned long nn, int isign);
	unsigned long i,i1,i2,i3,i4,np3;
	float c1=0.5,c2,h1r,h1i,h2r,h2i;
	double wr,wi,wpr,wpi,wtemp,theta;

	theta=3.141592653589793/(double) (n>>1);
	if (isign == 1) {
		c2 = -0.5;
		four1(data,n>>1,1);
	} else {
		c2=0.5;
		theta = -theta;
	}
	wtemp=sin(0.5*theta);
	wpr = -2.0*wtemp*wtemp;
	wpi=sin(theta);
	wr=1.0+wpr;
	wi=wpi;
	np3=n+3;
	for (i=2;i<=(n>>2);i++) {
		i4=1+(i3=np3-(i2=1+(i1=i+i-1)));
		h1r=c1*(data[i1]+data[i3]);
		h1i=c1*(data[i2]-data[i4]);
		h2r = -c2*(data[i2]+data[i4]);
		h2i=c2*(data[i1]-data[i3]);
		data[i1]=h1r+wr*h2r-wi*h2i;
		data[i2]=h1i+wr*h2i+wi*h2r;
		data[i3]=h1r-wr*h2r+wi*h2i;
		data[i4] = -h1i+wr*h2i+wi*h2r;
		wr=(wtemp=wr)*wpr-wi*wpi+wr;
		wi=wi*wpr+wtemp*wpi+wi;
	}
	if (isign == 1) {
		data[1] = (h1r=data[1])+data[2];
		data[2] = h1r-data[2];
	} else {
		data[1]=c1*((h1r=data[1])+data[2]);
		data[2]=c1*(h1r-data[2]);
		four1(data,n>>1,-1);
	}
}

int main() {
	
	
	// number of bits data is sampled
	int bits = 4;

	//sample rate of the signal 
	unsigned long sample_rate = 200; // Msps

	//number of samples you want to send for processing in the fft (any)
	//for example 100 samples
	//unsigned long captured_samples = 1 << 7; // this is 2**16
        unsigned long captured_samples = 1 << 8;
        
        cout <<"Captured Samples        " <<     captured_samples << endl;

        unsigned long Ksegments = 128;

        //frequency of the signal (has to be smaller than sample_rate/2)
        //for example 46
        double frequency = 4.0*sample_rate/captured_samples; // MHz
        frequency = 6.5;
        //frequency = (int)(frequency*captured_samples*sample_rate)*(double)sample_rate/captured_samples;
	
        generate_signal(data, fftdata, captured_samples, sample_rate, frequency);
        
        /*
	float *data    = (float*)malloc((captured_samples)*sizeof(float));
	float *fftdata = (float*)malloc((captured_samples)*sizeof(float));
	*/

	float *data    = (float*)malloc((Ksegments)*sizeof(float));
	float *fftdata = (float*)malloc((Ksegments)*sizeof(float));
        
        
	
//	quantisation(data, captured_samples, getMax(data,captured_samples), depthMax(bits) );
//	quantisation(fftdata, captured_samples, getMax(data,captured_samples), depthMax(bits) );

	

	//output_data (data, fftdata, captured_samples);
	
	// window data
	data[0] = 0.0;
	for(int k = 1 ; k < captured_samples/2; k++) {
		data[k] *= 2.0*(float)k/captured_samples;
		data[captured_samples - k] *= 2.0*(float)k/captured_samples;
	}


	realft(fftdata-1, captured_samples, 1);
        
        

	output_data (data, fftdata, captured_samples, sample_rate);
	
	//apply the FFT to the signal
	

	free(data);
	free(fftdata);
	return 0;
}