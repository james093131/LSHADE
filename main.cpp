#include "normal.hpp"

int main(int argc,const char *argv[])
{
    int run = atoi(argv[1]);//just run
    int iteration = atoi(argv[2]);//just iteration
    int pop = atoi(argv[3]);
    int DIM = atoi(argv[4]);
    int A = atoi(argv[5]);
    int H = atoi(argv[6]);
    int pbest = atoi(argv[7]);

    if( argc > 1 )
    {
        LSHADE lshade;
        lshade.RUN( iteration,pop, DIM, A, H, pbest);
    }
      
    else 
    {
        LSHADE lshade;
        lshade.RUN(1000,30,50,30,10,4);

    }
   
    return 0;
}