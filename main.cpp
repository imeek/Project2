
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
}

int main() {

    cout.precision(6);
    cout.setf(ios::fixed | ios::showpoint);
    
    /* Simulation initialisation variables */

    float sample_rate = 256; //Msps 
    //unsigned long captured_samples = 1 << 10; //1024 samples
    float captured_samples = 1024; //1024 samples

    /* Simulation data vectors */

    float *data = (float*) malloc((captured_samples) * sizeof (float));
    
    /* Simulation frequencies vector */
    
    float *f = (float*) malloc((sample_rate/2) * sizeof (float));
    int nf = 6;
    f[0] = 20; //10*sample_rate/captured_samples; //MHz
    f[1] =  40;
    f[2] =  60;
    f[3] =  80;
    f[4] =  100;
    f[5] =  120;
    
    //f[0] = 1.5;
    cout << "sample_rate :" << sample_rate << endl;
    cout << "captured_samples :" << captured_samples << endl;
    float delta = 0.0; 
    delta =  1/sample_rate;
    cout <<"Delta:   "<< delta << endl;
    cout << "sample_rate/captured_samples :" << (sample_rate/captured_samples) << endl;
    cout << "f[0] " << f[0] << endl;
    
    /* Generate Signals */
    

    tone(data, captured_samples, f, nf, delta);
    

// this is a test


    /* Output simulation results to data file */

    ofstream output("output.dat");
    for (int i = 0; i < captured_samples; i++) {
        output << data[i] << endl;
    }
    output.close();

    /* End simulation */

    return 0;
}