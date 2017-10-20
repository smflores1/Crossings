//this simulation measures fk cluster crossings with independently wired top/bottom sides and free left/right sides
#include "stdlib.h"
#include "stdio.h"
#include "math.h"

#define INTEGER 17
#define SIDELENGTHSUM 2048
#define RMAX 16
#define SEED 170
#define RUNFILE "RunFile.txt"
#define OUTFILE "RecFkQ3XingIndep.txt"
#define Q 3
#define PRINTFREQ 32767
#define RUNSMAX 4194304

void    randinit(long seed);
void    find(long x, long y, long s, long t);

#define M 16383
#define NewRandomInteger (++nd,ra[nd&M]=ra[(nd-471)&M]^ra[(nd-1586)&M]^ra[(nd-6988)&M]^ra[(nd-9689)&M])

long    ra[M+1], nd;
long    dx[4] = {1, 0, -1,  0}, dy[4] = {0, 1, 0, -1}, prob, TopFlag, BottomFlag;
int     length, height, lat[SIDELENGTHSUM][SIDELENGTHSUM]; //true array size is lat[length][height]

int main(void)
{
	long  x, y, l, m, runs, VertXing;
    FILE  *fprunfile, *fpoutfile;
    
    randinit(SEED);
    
    length = floor(SIDELENGTHSUM*1.0/(RMAX+1)+SIDELENGTHSUM*(RMAX-1)*(INTEGER-1)*1.0/(32*(RMAX+1)));
    height = SIDELENGTHSUM-length;
	
    prob =(int)(sqrt(Q)/(1+sqrt(Q))*2147483648.0);
    
    for (x = 1; x < length-1; ++x)
        lat[x][0] = lat[x][height-1] = 15;
    
    for (y = 0; y < height; ++y)
        lat[0][y] = lat[length-1][y] = 15;
    
    for (x = 1; x < length-1; ++x)
        for (y = 1; y < height-1; ++y)
        {
            if (y == 1 || y == height-2) lat[x][y] = Q;
            else lat[x][y] = Q+(int)((NewRandomInteger/2147483648.0)*Q);
        }
    
    VertXing = 0;
    
    for (runs = 0; runs <= RUNSMAX; ++runs)
    {
        
        for (x = 1; x < length-1; ++x)
            for (y = 1; y < height-1; ++y)
                lat[x][y] -= Q;
        
        for (x = 1; x < length-1; ++x)
            for (y = 1; y < height-1; ++y)
                if ((m = lat[x][y]) < Q)
                {
                    TopFlag = BottomFlag = 0;
                    l = Q+(int)((NewRandomInteger/2147483648.0)*Q);
                    lat[x][y] = l;
                    find(x, y, l, m);
                    if (TopFlag == 1 && BottomFlag == 1) ++VertXing;
                }
        
        if ((runs & PRINTFREQ) == 0){//"bitwise and;" note that 2^p-1 & 2^p = 0.
            
            fprunfile = fopen(RUNFILE, "w");
            fprintf(fprunfile, "%15ld\n", runs);
            fclose(fprunfile);
            
            fpoutfile = fopen(OUTFILE, "w");
            fprintf(fpoutfile, "%f%15f\n", (length-2)*1.0/(height-2), VertXing*1.0/runs); //aspect ratio is inverse of normal, wired sides in paper are top and bottom
            fclose(fpoutfile);
            
        }
        
    }/* for runs */
    
}


void find(long x, long y, long s, long t)
{
    int j, SideFlag;
    long xp, yp;
	
	if (y == 1)
		BottomFlag = 1; 
	if (y == height-2)
		TopFlag = 1;
	
    for (j = 0; j < 4; ++j)
    {
        
        xp = x+dx[j];
        yp = y+dy[j];
        
        SideFlag = 0;
        
        if (y == 1 || y == height-2) if (j == 0 || j == 2) SideFlag = 1;
        if (lat[xp][yp] == t)
            if (NewRandomInteger < prob || SideFlag == 1)
            {
                lat[xp][yp] = s;
                find(xp, yp, s, t);
            }
    }
}

void randinit(long seed)
{
    double a, ee = -1+1/2147483648.0;
    long i;
    extern long nd, ra[M+1];
    
    a = seed/2147483648.0;
    for (nd = 0; nd <= M; nd++)
    {
        a *= 16807;
        a += ee*(long)(a);
        if (a >= 1) a += ee;
        ra[nd] = (long)(2147483648.0*a);
    }
    for (i = 1; i<10000001; i++)
        NewRandomInteger;
}


