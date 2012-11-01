
#include <fstream>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <ctime>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

using namespace std;

#define SWAP(a,b) tempr=(a);(a)=(b);(b)=tempr
double pi = 4 * atan((double) 1);


void tone(float *data, int ndata, float *f, int nf, float delta) {
/* tone generation function */
    static float Tlast = 0.0;
    for (int k = 0; k < ndata; k++) {
        data[k] = 0;
        for (int i = 0; i < nf; i++) {
            data[k] += sin(2 * pi * f[i]*((k * delta)+ Tlast));
            
        }
        Tlast += ndata*delta;
      //  cout << "Tlast[" << k << "]  :" << Tlast << "   ndata: " << ndata << "   delta: " << delta << endl;
    }
    cout << "Tlast :" << Tlast << endl;
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
//this is a test
    cout.precision(6);
    cout.setf(ios::fixed | ios::showpoint);
    
    /* Simulation initialisation variables */

    unsigned long sample_rate = 256; //Msps 
    unsigned long captured_samples = 1 << 10; //1024 samples

    /* Simulation data vectors */

    float *data = (float*) malloc((captured_samples) * sizeof (float));
    float *data2 = (float*) malloc((captured_samples) * sizeof (float));   
    
    /* Simulation frequencies vector */
    
    float *f = (float*) malloc((sample_rate/2) * sizeof (float));

    
    //f[0] = 1.5;
    cout << "sample_rate :" << sample_rate << endl;
    cout << "captured_samples :" << captured_samples << endl;
    float delta = 0.0; 
    delta =  1/(float)sample_rate;
    cout <<"Delta:   "<< delta << endl;
    cout << "sample_rate/captured_samples :" << ((float)sample_rate/(float)captured_samples) << endl;

    
    /* Generate Signals */
    int nf = 6;
    f[0] = 20; 
    f[1] =  40;
    f[2] =  60;
    f[3] =  80;
    f[4] =  100;
    f[5] =  120;    
    tone(data, captured_samples, f, nf, delta);
    tone(data2, captured_samples, f, nf, delta);
    
    
    /* FFT data vector */
    
    float *fftdata = (float*)malloc((captured_samples)*sizeof(float));
    float *fftdata2 = (float*)malloc((captured_samples)*sizeof(float));
    
    for(unsigned long i=0; i<(captured_samples); i++)	{
	fftdata[i]=data[i];
        fftdata2[i]=data2[i];
    }
       
    realft(fftdata-1, captured_samples, 1);
    
    
        
    /* Output simulation results to data file */

    ofstream output("output.dat");
    for (int i = 0; i < captured_samples; i++) {
        output << data[i] << endl;
    }
    for (int i = 0; i < captured_samples; i++) {
        output << data2[i] << endl;
    }
    output.close();
    
    ofstream fftoutput("fftdata.dat");
    for (int i=0; i < captured_samples; i++)	{
            fftoutput << 0.5*i*(float)sample_rate/(float)captured_samples << "     ";
            fftoutput << fftdata[i] << endl;
    }
    fftoutput.close();

    /* End simulation */

    return 0;
}
