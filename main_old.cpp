
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
#define TAB '\t'

double pi = 4.0 * atan((double) 1);

/* like abs() but for doubles. necessary because */
/* abs() only accepts and returns ints */
double flabs(float f)	{
return f * (f<0? -1: 1);
}

/* returns the maximum value at a given bit-depth */
int depthMax(int bits)	{
return (int)pow(2.0, bits-1);
}

/* returns maximum value in an array of doubles */
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

/* www.strchr.com/standard_deviation_in_one_pass */
/* - this version takes two */
double std_dev1(double a[], int n) {
    if(n == 0)
        return 0.0;
    double sum = 0;
    for(int i = 0; i < n; ++i)
       sum += a[i];
    double mean = sum / n;
    double sq_diff_sum = 0;
    for(int i = 0; i < n; ++i) {
       double diff = a[i] - mean;
       sq_diff_sum += diff * diff;
    }
    double variance = sq_diff_sum / n;
    return sqrt(variance);
}

void quantisation(float data[], long captured_samples,int dataMax, int depthMax )	{

for(int i=0; i<captured_samples; i++)	{	
data[i]= floor(depthMax*(data[i]/dataMax));	
}

}

void tone(float *data, int ndata, float *f, int nf, float delta) {
/* tone generation function */
    static double Tlast = 0.0;
    for (int k = 0; k < ndata; k++) {
        data[k] = 0;
        for (int i = 0; i < nf; i++) {
            data[k] += sin(2 * M_PI * f[i]*(k * delta + Tlast));
        }
    }
    Tlast += ndata*delta;
}

void four1(float data[], unsigned long nn, int isign)
/*******************************************************************************
Replaces data[1..2*nn] by its discrete Fourier transform, if isign is input as
1; or replaces data[1..2*nn] by nn times its inverse discrete Fourier transform,
if isign is input as -1. data is a complex array of length nn or, equivalently,
a real array of length 2*nn. nn MUST be an integer power of 2 (this is not
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
for (m=1;m<mmax;m+=2) {	/* Here are the two nested inner loops. */
for (i=m;i<=n;i+=istep) {
j=i+mmax;	/* This is the Danielson-Lanczos formula. */
tempr=wr*data[j]-wi*data[j+1];
tempi=wr*data[j+1]+wi*data[j];
data[j]=data[i]-tempr;
data[j+1]=data[i+1]-tempi;
data[i] += tempr;
data[i+1] += tempi;
}
wr=(wtemp=wr)*wpr-wi*wpi+wr;	/* Trigonometric recurrence. */
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
    cout.precision(6);
    cout.setf(ios::fixed | ios::showpoint);
    
    /* Simulation initialisation variables */
    unsigned long sample_rate = 150; //Msps
    unsigned long captured_samples = 1 << 10; //1024 samples
    
    // number of bits data is sampled
    int bits = 14;
    ofstream output("output.dat");
    const int tsegments = 1000;
    
    /* Noise parameters */
    const float Po = 1E-3; // Watts
    const float Vpp = 4; // Volts
    const float Vmax = 2; // Volts
    const float Rt = 50; // Ohms
    float Pfs = (Vmax*Vmax)/(2*Rt); // full scale Power
    float Noise_dBm = 0.0;
    float ADC_squared = 0.0; //to replace
    
    Noise_dBm = -260.0*log10(2.0)+10*log10(Pfs/Po);
    cout << "Noise_dBm :" << Noise_dBm << endl;
    Noise_dBm +=log10(ADC_squared); //to replace
            
    /* Signal Parameters */
    float *f = (float*) malloc((sample_rate/2) * sizeof (float));
    float delta = 1/(float)sample_rate;
    

    double N_dBm[750];
    for (int l=0;l<750;l++)
    {
        N_dBm[l] = 0.0;
    }
    
    double Tseg[tsegments];
 
    double freq;
    /* Simulation frequency vector */
    int nf = 1;
    f[0] = 0.0;
    freq = f[0];
    
    for (int j=0;j<3;j++)
    {
        f[0] = f[0] +0.1;          
        for (int t=0 ;t<tsegments ;t++ )
        {
            /* Simulation data vectors */
            float *data = (float*) malloc((captured_samples) * sizeof (float));

            /* Generate Signal */
            tone(data, captured_samples, f, nf, delta);

            /* Quantise Data*/
            quantisation(data, captured_samples, getMax(data,captured_samples), depthMax(bits) );

            /* FFT data vector */
            float *fftdata = (float*)malloc((captured_samples)*sizeof(float));
            for(unsigned long i=0; i<(captured_samples); i++)
            {
                fftdata[i]=data[i];
            }

            /*Window data*/
            fftdata[0] = 0.0;
            for(int k = 1 ; k < captured_samples/2; k++) 
            {
                fftdata[k] *= 2.0*(float)k/captured_samples;
                fftdata[captured_samples - k] *= 2.0*(float)k/captured_samples;
            }

            /* FFT data */
            realft(fftdata-1, captured_samples, 1);

            /* Normalise FFT data */
            for(int k = 0.0; k< captured_samples; k++)fftdata[k] *= sqrt(2.0)/captured_samples;

            /* Output simulation results to data file */

            for(int k = 0; k < nf; k++) 
            {
                int idx = (int)trunc(f[k]/sample_rate*captured_samples);
                Tseg[t] = sqrt(fftdata[2*idx+1]*fftdata[2*idx+1]+fftdata[2*idx]*fftdata[2*idx]);
            }
            
            /* delete data vectors */
            free(data);
            free(fftdata);

          }
            
          N_dBm[j] = std_dev1(Tseg, tsegments);
        
        /* */
    }
    
    ofstream dBm_out("dBm_out.dat");
    for (int l=0;l<750;l++)
    {
        freq = freq + 0.1;
        dBm_out << freq << TAB << N_dBm[l] << endl;
    }
    
    
    /* End simulation */
    output.close();
    return 0;
}
