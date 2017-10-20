//this simulation measures spin cluster crossings in a hexagon with mutually wired bottom/upper-right/upper-left sides and fluctuating top/bottom-left/bottom-right sides
#include "stdlib.h"
#include "stdio.h"
#include "math.h"

#define INTEGER 17 //between 1 and 33, inclusive
#define SIDELENGTHSUM 1024
#define RMAX 16
#define SEED 170
#define RUNFILE "RunFile.txt"
#define	OUTFILE "HexSpinQ4XingMutual.txt"
#define Q 4
#define PRINTFREQ 100//32767
#define RUNSMAX 4194304
#define PROB 0.5

void    randinit(long seed);
void    find(long x, long y, long s, long t);
void	walk(void);

#define M  16383
#define NewRandomInteger (++nd,ra[nd&M] = ra[(nd-471)&M]^ra[(nd-1586)&M]^ra[(nd-6988)&M]^ra[(nd-9689)&M])

long 	ra[M+1], nd;
long	dx[6] = {-1, -1, 0, 1, 1, 0}, dy[6] = {0, 1, 1, 0, -1, -1}, prob, NoXing, BottomUpperRightXing, BottomUpperLeftXing, UpperLeftUpperRightXing, AllXing;
long	L1, L2, L3, L4, L5, L6; //side lengths of the hexagon
long	Vertex1x, Vertex2x, Vertex3x, Vertex4x, Vertex5x, Vertex6x;
long    Vertex1y, Vertex2y, Vertex3y, Vertex4y, Vertex5y, Vertex6y;
long	length, height, lat[SIDELENGTHSUM+3][SIDELENGTHSUM+3]; //equivalently, lat[length][height]

int main(void)

{
    
	long 	x, y, l, m, runs;
    FILE	*fprunfile, *fpoutfile;
    
    randinit(SEED);
	
	L2 = floor(SIDELENGTHSUM*1.0/(RMAX+1)+SIDELENGTHSUM*(RMAX-1)*(INTEGER-1)*1.0/(32*(RMAX+1)));
	L1 = SIDELENGTHSUM-L2;

    //this causes the hexagon side length to alternate
	L3 = L1;
	L4 = L2;

    //this closes the hexagon
    if (L4 > L1+L2) {
        printf("error, hexagon does not close\n");
        exit(1);
    }
    else L5 = L1+L2-L4;
    
    if (L1 > L3+L4) {
        printf("error, hexagon does not close\n");
        exit(1);
    }
    else L6 = L3+L4-L1;
	
	Vertex1x = L1+1;
	Vertex1y = 1;
	Vertex2x = L1+L2+1;
	Vertex2y = 1;
	Vertex3x = L1+L2+1;
	Vertex3y = L3+1;
	Vertex4x = L1+L2-L4+1;
	Vertex4y = L3+L4+1;
	Vertex5x = 1;
	Vertex5y = L3+L4+1;
	Vertex6x = 1;
	Vertex6y = L1+1;
	
	length = L1+L2+3; //this is the length of the array for the rectangle that contains the hexagon
	height = L3+L4+3; //this is the height of the array for the rectangle that contains the hexagon
	
	prob = (long) (2147483648.0*PROB);
	
	for (x = 0; x < length; ++x)
		for (y = 0; y < height; ++y)
        {
			if (x+y <= Vertex1x) lat[x][y] = 15; //always blocked
            else if (x > 1 && x <= Vertex1x && x+y == Vertex1x+1) lat[x][y] = Q+1+(int)((NewRandomInteger/2147483648.0)*(Q-1)); //always between 1 and Q-1
			else if (Vertex1x <= x && y == 0) lat[x][y] = 15; //always blocked
            else if (Vertex1x < x && x < length-1 && y == 1) lat[x][y] = Q; //always 0
			else if (y <= Vertex3y && x == length-1) lat[x][y] = 15; //always blocked
            else if (y > 1 && y <= Vertex3y && x == length-2) lat[x][y] = Q+1+(int)((NewRandomInteger/2147483648.0)*(Q-1)); //always between 1 and Q-1
			else if (x+y >= length+L3) lat[x][y] = 15; //always blocked
            else if (Vertex4x <= x && x < length-2 && x+y == length+L3-1) lat[x][y] = Q; //alway 0
			else if (x <= Vertex4x && y == height-1) lat[x][y] = 15; //alwyas blocked
            else if (x > 0 && x < Vertex4x && y == height-2) lat[x][y] = Q+1+(int)((NewRandomInteger/2147483648.0)*(Q-1)); //always between 1 and Q-1
			else if (Vertex6y <= y && x == 0) lat[x][y] = 15; //always blocked
            else if (Vertex6y <= y && y < Vertex5y && x == 1) lat[x][y] = Q; //always 0
			else lat[x][y] = Q+(int)((NewRandomInteger/2147483648.0)*Q);
		}
    
	NoXing = AllXing = BottomUpperRightXing = BottomUpperLeftXing = UpperLeftUpperRightXing = 0;
	
	for (runs = 1; runs <= RUNSMAX; ++runs)
	{
		
		for (x = 0; x < length; ++x)
			for (y = 0; y < height; ++y)
				if (lat[x][y] != 15) lat[x][y] -= Q;
        
        for (x = Vertex1x+1; x < length-1; ++x)
            if ((m = lat[x][1]) < Q)
            {
                y = 1;
                l = Q;
                lat[x][y] = l;
                find(x, y, l, m);
            }
        
        for (x = Vertex4x; x < length-2; ++x)
            if ((m = lat[x][length+L3-1-x]) < Q)
            {
                y = length+L3-1-x;
                l = Q;
                lat[x][y] = l;
                find(x, y, l, m);
            }
        
        for (y = Vertex6y; y < Vertex5y; ++y)
            if ((m = lat[1][y]) < Q)
            {
                x = 1;
                l = Q;
                lat[x][y] = l;
                find(x, y, l, m);
            }
        
        for (x = 2; x <= Vertex1x; x++)
            if ((m = lat[x][Vertex1x+1-x]) < Q)
            {
                y = Vertex1x+1-x;
                l = Q+1+(int)((NewRandomInteger/2147483648.0)*(Q-1));
                lat[x][y] = l;
                find(x, y, l, m);
            }
        
        for (y = 2; y <= Vertex3y; ++y)
            if ((m = lat[length-2][y]) < Q)
            {
                x = length-2;
                l = Q+1+(int)((NewRandomInteger/2147483648.0)*(Q-1));
                lat[x][y] = l;
                find(x, y, l, m);
            }
        
        for (x = 1; x < Vertex4x; ++x)
            if ((m = lat[x][height-2]) < Q)
            {
                y = height-2;
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

        if ((runs % PRINTFREQ) == 0){//"bitwise and;" note that 2^p-1 & 2^p = 0.
			
			fprunfile = fopen(RUNFILE, "w");
			fprintf(fprunfile, "%15ld\n", runs);
			fclose(fprunfile);
            
            NoXing = runs-AllXing-BottomUpperRightXing-BottomUpperLeftXing-UpperLeftUpperRightXing;
			
			fpoutfile = fopen(OUTFILE, "w");
			fprintf(fpoutfile, "%f%15f\n", L2*1.0/L1, NoXing*1.0/runs);
			fprintf(fpoutfile, "%f%15f\n", L2*1.0/L1, AllXing*1.0/runs);
			fprintf(fpoutfile, "%f%15f\n", L2*1.0/L1, BottomUpperRightXing*1.0/runs);
			fprintf(fpoutfile, "%f%15f\n", L2*1.0/L1, BottomUpperLeftXing*1.0/runs);
			fprintf(fpoutfile, "%f%15f\n", L2*1.0/L1, UpperLeftUpperRightXing*1.0/runs);
			fclose(fpoutfile);
		}
		
	}// for runs
	
}

void find(long x, long y, long s, long t)
{
    int	j, BottomUpperLeftUpperRightFlag, TopLowerLeftLowerRightFlag;
	long xp, yp;
	
    for (j = 0; j < 6; ++j)
    {
        
        xp = x+dx[j];
		yp = y+dy[j];

        BottomUpperLeftUpperRightFlag = TopLowerLeftLowerRightFlag = 0;

		if (lat[xp][yp] == t)
        {
            if (Vertex1x < x && x < length-1 && y == 1) BottomUpperLeftUpperRightFlag = 1;
            if (Vertex4x <= x && x < length-2 && x+y == length+L3-1) BottomUpperLeftUpperRightFlag = 1;
            if (Vertex6y <= y && y < Vertex5y && x == 1) BottomUpperLeftUpperRightFlag = 1;
            
            if (x > 1 && x <= Vertex1x && x+y == Vertex1x+1) TopLowerLeftLowerRightFlag = 1;
            if (y > 1 && y <= Vertex3y && x == length-2) TopLowerLeftLowerRightFlag = 1;
            if (x > 0 && x < Vertex4x && y == height-2) TopLowerLeftLowerRightFlag = 1;
            
            if ((BottomUpperLeftUpperRightFlag == 1 && s == Q) || (TopLowerLeftLowerRightFlag == 1 && s != Q) || (BottomUpperLeftUpperRightFlag == 0 && TopLowerLeftLowerRightFlag == 0))
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
	
	long x, y, xp, yp, dir, End1x, End1y, End2x, End2y;
	
	dir = 40;
	x = Vertex1x+1;
	y = Vertex1y;
	
	do {	//first spin cluster perimeter walk
        
		xp = x+dx[dir%6];
		yp = y+dy[dir%6];
        
		if (lat[xp][yp] == Q)
        {
            x = xp;
            y = yp;
            --dir;
        }
		else ++dir;
        
	} while ((x != Vertex2x || y != Vertex2y) && (x != Vertex4x || y != Vertex4y) && (x != Vertex6x || y != Vertex6y));
    
    End1x = x;
	End1y = y;
	
	dir = 62;
	x = Vertex5x;
	y = Vertex5y-1;
    
    if (End1x != Vertex4x || End1y != Vertex4y) {
		do { //second walk
			
			xp = x+dx[dir%6];
			yp = y+dy[dir%6];
            
			if (lat[xp][yp] == Q)
            {
				x = xp;
				y = yp;
				--dir;
			}
			else ++dir; //deactivated
			
			
		} while ((x != Vertex2x || y != Vertex2y) && (x != Vertex4x || y != Vertex4y) && (x != Vertex6x || y != Vertex6y));
        
		End2x = x;
		End2y = y;
	}
    
    if (End1x == Vertex2x && End1y == Vertex2y && End2x == Vertex6x && End2y == Vertex6y) ++NoXing;
    else if (End1x == Vertex6x && End1y == Vertex6y && End2x == Vertex4x && End2y == Vertex4y) ++AllXing;
    else if (End1x == Vertex4x && End1y == Vertex4y) ++BottomUpperRightXing;
    else if (End1x == Vertex6x && End1y == Vertex6y && End2x == Vertex2x && End2y == Vertex2y) ++BottomUpperLeftXing;
    else if (End1x == Vertex2x && End1y == Vertex2y && End2x == Vertex4x && End2y == Vertex4y) ++UpperLeftUpperRightXing;
    else {printf("fail  %10ld%10ld%10ld%10ld\n", End1x, End1y, End2x, End2y); exit(1);}
    
}
	
	
void randinit(long seed)
{	
	double a, ee = -1 + 1/2147483648.0;
	long i;
	extern long nd, ra[M+1];
	
	a = seed/2147483648.0;
	for (nd = 0; nd <= M; nd++)
	{
		a *= 16807;
		a += ee * (long)(a);
		if (a >= 1) a += ee;
		ra[nd] = (long) (2147483648.0 * a);
	}
	nd = M;
	for(i = 0; i<100001; i++)
		NewRandomInteger;
}



