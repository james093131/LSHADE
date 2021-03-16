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

        int Archieve_Coef;
        int H_Table_Coef;
    
    private:
    void INI(int ITER,int POP,int DIM,int A,int H)
    {
        Particle.clear();
        Particle.swap((Particle));

        Objective_Value.clear();
        Objective_Value.swap(Objective_Value);

        Archieve.clear();
        Archieve.swap(Archieve);
        
        Archieve_Coef = 0;
        H_Table_Coef = 0;

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

    void Mutation_Selection(int pbest,int A)
    {
        d2d S_Table;
        for(int i=0;i<Particle.size();i++)
        {
            int r = rand() % (H_Table.size() - 0 + 1) + 0;
            double CR = Normal_Distribution(i,r);
            double F = Cauchy_Distribution(i,r);
            d1d V(Particle[i].size());
            for(int j=0;j<Particle[i].size();j++)
            {
                double a = ((float(rand()) / float(RAND_MAX)) * (1.0 - 0.0)) + 0.0;
                if (a < CR)
                    V[j] = Mutation( pbest, CR, F,i,j);
                else 
                    V[j] = Particle[i][j];
            }
            double V_Objective_Value = Function_Evaluate(V.size(),V);
            if( V_Objective_Value < Objective_Value[i])
            {
                S_Table.push_back({CR,F,Objective_Value[i] - V_Objective_Value});
                
            
                

                Particle[i].assign(V.begin(), V.end());
                
                Objective_Value[i] = V_Objective_Value;

                if(Archieve_Coef == A)///這邊要改成random pop
                {
                    Archieve_Coef = 0;
                }
                
                Archieve[Archieve_Coef].assign(V.begin(),V.end());
                Archieve_Coef ++;
            }
        }
        

    }

    double Mutation(int pbest,double CR,double F,int Pop,int Dim)
    {
        int X_P_best = Current_To_Pbest(pbest);

        int r1 = rand() % ( (Objective_Value.size() -1) - 0 + 1) + 0;
        while( r1 == Pop)
        {
            r1 = rand() % ( (Objective_Value.size() -1) - 0 + 1) + 0;
        }
        
        int r2 = rand() % ( (Objective_Value.size() -1) - 0 + 1) + 0;
        while (r2 ==r1 || r2 == Pop)
        {
            r2 = rand() % ( (Objective_Value.size() -1) - 0 + 1) + 0;
        }
        
        double V = Particle[Pop][Dim] + F * (Particle[X_P_best][Dim] - Particle[Pop][Dim]) + F *(Particle[r1][Dim] - Particle[r2][Dim]); 
        return V;
    }
    int  Current_To_Pbest(int pbest)
    {
        int r = rand() % ( (pbest-1) - 0 + 1) + 0;
        return Objective_Rank_INDEX[r] ;
    }

    void Update_Htable(int H,d2d S_table)//next cccccode 
    {
        double S_CR_SUM = 0;
        double S_F_SUM = 0;
        double S_FIT_SUM = 0;
        for(int i=0;i<S_table.size();i++)
        {
            S_CR_SUM += S_table[i][0];
            S_F_SUM += S_table[i][1];
            S_FIT_SUM += S_table[i][2];
        }
        double s1 = 0;
        double s2 = 0;
        double f1 = 0;
        double f2 = 0;

        for(int i=0;i<S_table.size();i++)
        {
            s1 += (S_table[i][2]/S_FIT_SUM)*pow(S_table[i][0],2);
            s2 += (S_table[i][2]/S_FIT_SUM)* S_table[i][0];

            f1 += (S_table[i][2]/S_FIT_SUM)*pow(S_table[i][1],2);
            f2 += (S_table[i][2]/S_FIT_SUM)* S_table[i][1];
        }

        if(H_Table_Coef == H)
            H_Table_Coef = 0;

        H_Table[H_Table_Coef][0] = s1/s2;

        H_Table[H_Table_Coef][1] = f1/f2; 
        
        H_Table_Coef++;

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
    double Function_Evaluate(int DIM,d1d arr)
    {
        return ACKLEY_OBJECTIVE_VALUE( DIM, arr);
    }

};