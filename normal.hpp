#include<stdio.h>
#include<fstream>
#include<iostream>
#include<sstream>
#include<stdlib.h>
#include <string.h>
#include<time.h>
#include<vector>
#include <algorithm>
#include<math.h>
#include<float.h>
#include <sys/stat.h>

using namespace std;

typedef vector<char> c1d;
typedef vector<c1d> c2d;
typedef vector<c2d> c3d;
typedef vector<c3d> c4d;
typedef vector<int> i1d;
typedef vector<i1d> i2d;
typedef vector<i2d> i3d;
typedef vector<i3d> i4d;
typedef vector<double>d1d;
typedef vector<d1d> d2d;
typedef vector<d2d> d3d;
typedef vector<d3d> d4d;


class LSHADE{
    public:
        void RUN(int ITER,int POP,int DIM,int A,int H,int pbest)
        {
            srand( time(NULL) );
            INI(ITER,POP,DIM,A,H);
            Particle_INI(DIM);
            RANK();
           
        }
    private:
        d2d Particle;
        d1d Objective_Value;
        d2d Archieve;
        d2d H_Table;
        d1d Objective_Rank_INDEX;
    
    
    private:
    void INI(int ITER,int POP,int DIM,int A,int H)
    {
        Particle.clear();
        Particle.swap((Particle));

        Objective_Value.clear();
        Objective_Value.swap(Objective_Value);

        Archieve.clear();
        Archieve.swap(Archieve);
        
      

        H_Table.clear();
        H_Table.swap(H_Table);

        Objective_Rank_INDEX.clear();
        Objective_Rank_INDEX.swap(Objective_Rank_INDEX);

        Particle.resize(POP,d1d(DIM));

        Objective_Value.resize(POP);

        Archieve.resize(A,d1d(DIM));

        H_Table.resize(H,d1d(2,0.5));

        Objective_Rank_INDEX.resize(POP);
        for(int i=0;i<POP;i++)
        {
            Objective_Rank_INDEX[i] = i;
        }
    }

    void Mutation()
    {
        d2d S_Table;
        for(int i=0;i<Particle.size();i++)
        {
            int r = rand() % (H_Table.size() - 0 + 1) + 0;
            double CR = Normal_Distribution(i,r);
            double F = Cauchy_Distribution(i,r);

        }
        

    }
    void Particle_INI(int DIM)
    {
        for(int i=0;i<Particle.size();i++)
        {
            ACKLEY(DIM,i);
        }
    }
    double Normal_Distribution(int i ,int r)//公式待確認
    {
        double z = (5 - -5) * rand() / (RAND_MAX + 1.0) + -5;
        double x = H_Table[r][i]+z*0.1;
        if (x > 1)
            x = 1 ;
        else if (x < 0)
            x = 0;


        return x;
    }
    double Cauchy_Distribution(int i ,int r)//公式待確認
    {
        double z = (5 - -5) * rand() / (RAND_MAX + 1.0) + -5;
        double x = H_Table[r][i]+z*0.1;
        if (x > 1)
            x = 1 ;
        else if (x < 0)
            x = 0;
        return x;       
    }
    void pairsort(d1d a, d1d &b) 
    { 
        int n =a.size();
        pair<double, double> pairt[n]; 
    
        // Storing the respective array 
        // elements in pairs. 
        for (int i = 0; i < n; i++)  
        { 
            pairt[i].first = a[i]; 
            pairt[i].second = b[i]; 
        } 
    
        // Sorting the pair array. 
        sort(pairt, pairt + n); 
        
        // Modifying original arrays 
        for (int i = 0; i < n; i++)  
        { 
            b[i] = pairt[i].second; 
        } 
    } 
    void RANK()
    {
      pairsort(Objective_Value,Objective_Rank_INDEX); 
    }
    void ACKLEY(int DIM,int index) //random initial in ACKLEY Function and using RADVIZ calculate 2 dimension coordinates
    {
        double max = 40.0;
        double min = -40.0;
   
        for(int i=0;i<DIM;i++)
        {
            double a = ((float(rand()) / float(RAND_MAX)) * (max - min)) + min;
            Particle[index][i] = a;

        }
        double F = ACKLEY_OBJECTIVE_VALUE(DIM,Particle[index]);
        Objective_Value[index] = F;
    }
    double ACKLEY_OBJECTIVE_VALUE(int DIM,d1d arr) //Calculate the objective value at ACKLEY Function 
    {
        double sum1= 0;
        double sum2 = 0;
        for(int i=0;i<DIM;i++)
        {

            sum1 += pow(arr[i],2);
            sum2 += cos(2*M_PI*arr[i]);

        }
        double F = -20*(exp((-0.2)*sqrt(sum1/DIM)))-exp(sum2/DIM)+20+exp(1);
        return F;
    }
};