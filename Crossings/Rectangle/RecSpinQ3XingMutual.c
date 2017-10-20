//this simulation measures spin cluster crossings with mutually wired top/bottom sides and fluctuating left/right sides
#include "stdlib.h"
#include "stdio.h"
#include "math.h"

#define INTEGER 17
#define SIDELENGTHSUM 2048
#define RMAX 16
#define SEED 170
#define RUNFILE "RunFile.txt"
#define OUTFILE "RecSpinQ3XingMutual.txt"
#define Q 3
#define PRINTFREQ 32767
#define RUNSMAX 4194304

void    randinit(long seed);
void    find(long x, long y, long s, long t);
void	walk(void);

#define M 16383
#define NewRandomInteger (++nd,ra[nd&M]=ra[(nd-471)&M]^ra[(nd-1586)&M]^ra[(nd-6988)&M]^ra[(nd-9689)&M])

long    ra[M+1], nd;
long    dx[4] = {1, 0, -1,  0}, dy[4] = {0, 1, 0, -1}, prob, VertXing;
int     length, height, lat[SIDELENGTHSUM][SIDELENGTHSUM]; //true array size is lat[length][height]

int main(void)
{
	long  x, y, l, m, runs;
    FILE  *fprunfile, *fpoutfile;
    
    randinit(SEED);
    
    length = floor(SIDELENGTHSUM*1.0/(RMAX+1)+SIDELENGTHSUM*(RMAX-1)*(INTEGER-1)*1.0/(32*(RMAX+1)));
    height = SIDELENGTHSUM-length;
	
    prob = (int)(sqrt(Q)/(1+sqrt(Q))*2147483648.0);
    
    for (x = 1; x < length-1; ++x)
        lat[x][0] = lat[x][height-1] = 15;
    
    for (y = 0; y < height; ++y)
        lat[0][y] = lat[length-1][y] = 15;
    
    for (x = 1; x < length-1; ++x)
        for (y = 1; y < height-1; ++y)
        {
            if (y == 1 || y == height-2) lat[x][y] = Q;
            else if (x == 1 || x == length-2) lat[x][y] = Q+1+(int)((NewRandomInteger/2147483648.0)*(Q-1));
            else lat[x][y] = Q+(int)((NewRandomInteger/2147483648.0)*Q);
        }
    
    VertXing = 0;
    
    for (runs = 0; runs <= RUNSMAX; ++runs)
    {
        
        for (x = 1; x < length-1; ++x)
            for (y = 1; y < height-1; ++y)
                lat[x][y] -= Q;
        
        y = 1;
        for (x = 1; x < length-1; ++x)
            if ((m = lat[x][y]) < Q)
            {
                l = Q;
                lat[x][y] = l;
                find(x, y, l, m);
            }
        
        y = height-2;
        for (x = 1; x < length-1; ++x)
            if ((m = lat[x][y]) < Q)
            {
                l = Q;
                lat[x][y] = l;
                find(x, y, l, m);
            }
        
        x = 1;
        for (y = 1; y < height-1; ++y)
            if ((m = lat[x][y]) < Q)
            {
                l = Q+1+(int)((NewRandomInteger/2147483648.0)*(Q-1));
                lat[x][y] = l;
                find(x, y, l, m);
            }
        
        x = length-2;
        for (y = 1; y < height-1; ++y)
            if ((m = lat[x][y]) < Q)
            {
                l = Q+1+(int)((NewRandomInteger/2147483648.0)*(Q-1));
                lat[x][y] = l;
                find(x, y, l, m);
            }
        
        for (x = 1; x < length-1; ++x)
            for (y = 1; y < height-1; ++y)
                if ((m = lat[x][y]) < Q)
                {
                    l = Q+(int)((NewRandomInteger/2147483648.0)*Q);
                    lat[x][y] = l;
                    find(x, y, l, m);
                }
        
        walk();
        
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
    int j, TopBottomFlag, LeftRightFlag;
    long xp, yp;
	
    for (j = 0; j < 4; ++j)
    {
        
        xp = x+dx[j];
        yp = y+dy[j];
        
        TopBottomFlag = LeftRightFlag = 0;
        
		if (lat[xp][yp] == t)
        {
            if (yp == 1 || yp == height-2) TopBottomFlag = 1;
            if ((xp == 1 || xp == length-2) && TopBottomFlag != 1) LeftRightFlag = 1;
			if ((TopBottomFlag == 1 && s == Q) || (LeftRightFlag == 1 && s != Q) || (LeftRightFlag == 0 && TopBottomFlag == 0))
				if (NewRandomInteger < prob)
				{
					lat[xp][yp] = s;
					find(xp, yp, s, t);
				}
        }
    }
}


void walk(void)
{
	
	long x, y, xp, yp, dir;
	
	dir = 41;
	x = 1;
	y = 1;
	
	do {	//spin cluster perimeter walk
        
		xp = x+dx[dir%4];
		yp = y+dy[dir%4];
        
		if (lat[xp][yp] == Q)
        {
            x = xp;
            y = yp;
            ++dir;
        }
        
		else --dir;
        
	} while ((x != 1 || y != height-2) && (x != length-2 || y != 1)); //while we haven't reached the top left or bottom right corner
    
    if (x == 1 && y == height-2) ++VertXing; //vertical spin cluster crossing
    
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


