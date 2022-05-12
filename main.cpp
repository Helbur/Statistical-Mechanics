#include<iostream>
#include<vector>
#include<cstdlib>
#include<math.h>
#include<complex>
#include<ctime>
#include<fstream>

using namespace std;

#define PI 3.141592653589793238462643
const int q=3;      // q spin states
const int L=16;      // linear system size
double T=0.25;   // temperature in units of J

const int N=L^2; // total number of spins
double pconnect;//1-exp(-1/T); // connection probability

const int NCLUSTERS=1; // # of cluster builds in one MC step
const int NESTEPS=10000; // # of equilibrium MC step
const int NMSTEPS=10000; // # of measurement MC step
const int NBINS=100; // # of measurement bins

vector<int> S(N); // the spin array
vector<int> M(q); // number of spins in different states
vector<complex<double> > W(q); // order parameter weights

//lattice handling;
enum dirs{RIGHT, LEFT, UP, DOWN};
int indx(int x, int y){return y*L+x;} // make an indx on every site
int xpos(int i){return i%L;}
int ypos(int i){return i/L;}

int Nbr(int i, int dir)
{
    int x=xpos(i);
    int y=ypos(i);
    switch(dir)
    {
        case RIGHT: return indx((x+1)%L,y);
        case LEFT: return indx((x-1+L)%L,y);
        case UP: return indx(x,(y+1)%L);
        case DOWN: return indx(x,(y-1+L)%L);
    }
}

void FlipandBuildFrom(int s)
{
    int oldstate(S[s]), newstate((S[s]+1)%q);
    S[s]=newstate; // flip spin
    M[oldstate]--; M[newstate]++; // update spin counts

    for(int dir=0; dir<4; dir++) // go thru neighbors
    {
        int j=Nbr(s, dir);
        if(S[j] == oldstate)
        if(rand()/(RAND_MAX+1.) < pconnect) {FlipandBuildFrom(j);}
    }
}

int main()
{
    cout << "main called successfully" << endl;
    ofstream file;
    file.open("correlator.txt");
    //initialize order parameter weights
    for(int s=0; s<q; s++)
        W[s]=complex<double>(cos(2*PI*s/q), sin(2*PI*s/q));

    for(int n=0; n<NBINS; n++)
    {
        for(int i=0; i<N; i++) S[i]=0; //initialize spins
        for(int s=1; s<q; s++) M[s]=0; // initialize counters
        M[0]=N;
        srand((unsigned) time(0)); // initialize RNG
        // equilibrate
        //for(int t=0; t<NESTEPS; t++)
        //    for(int c=0; c<NCLUSTERS; c++)
        //    {
        //        FlipandBuildFrom(rand()%N);
        //    }
        //measure
        T = 0.0001 + (4-0.0001)*double(n)/NBINS;
        pconnect = 1-exp(-1./T);
        complex<double> m(0.,0.);
        //complex<double> m0(0.,0.);
        //complex<double> mr(0.,0.);
        double m1=0, m2=0, m4=0; // measurement results

        for(int t=0; t<NMSTEPS; t++)
        {

            for(int c=0; c<NCLUSTERS; c++) FlipandBuildFrom(rand()%N);
            complex<double> tm(0.,0.);
            //complex<double> m0(0.,0.);
            //complex<double> mr(0.,0.);
            for(int s=0; s<q; s++){tm+=W[s]*double(M[s]);}
            tm/=N;
            double tm1=abs(tm);
            double tm2=tm1*tm1;
            //m0+=conj(W[S[0]]);
            //mr+=W[S[r]];
            m+=tm; m1+=tm1; m2+=tm2; m4+=tm2*tm2;
        }
        m/=NMSTEPS;//m0/=NMSTEPS; mr/=NMSTEPS;// m1/=NMSTEPS; m2/=NMSTEPS; m4/=NMSTEPS;
    file << std::real(m) << endl;
    //cout << std::real(m) << endl;
    }
    file.close();
}
