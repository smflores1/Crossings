//In this simulation, side one is the bottom side of the hexagon, side two is the upper left side, and side three is the lower right side.  This may be backwards from some notes with 
//sides two and three switched.  The vertices are numbered counterclockwise starting from the lower left vertex.
#include "stdlib.h"
#include "stdio.h"
#include "math.h"

void randinit(long seed);

#define integer 17
#define RUNSMAX 3276805//4194304//262144//8388608//4194304//32768//1048576//2147483647//131072//
#define PRINTFREQ2 32767//131071//32767
#define PROB 0.5
#define DIRMAX 6
#define S 1048575//32767
#define SEED 200
#define M  16383					
#define BLOCK 2147483647  //"blocked," we never determine if the site is occupied or vaccant
#define RUNFILE "RunFile.txt"
#define	OUTFILE "HexPercXing.txt"

#define NewRandomInteger (++nd,ra[nd&M] = ra[(nd-471)&M]^ra[(nd-1586)&M]^ra[(nd-6988)&M]^ra[(nd-9689)&M])								 

float	r;
long	L1, L2, L3, L4, L5, L6;
long	vertex1, vertex2, vertex3, vertex4, vertex5, vertex6;
long	length, height;
long	lat[4192][4192]; //true array size is lat[L1+L2+3][L3+L4+3]
long 	ra[M+1], nd;

int main(void)

{

	L2 = floor(120.471+56.4706*(integer-1));
	L1 = 2048-L2;
	L3 = L1; 
	L4 = L2; 
	L5 = L1+L2-L4;
	L6 = L3+L4-L1;
	
	vertex1 = L1+1;
	vertex2 = L1+L2+1;
	vertex3 = L3+1;
	vertex4 = L1+L2-L4+1;
	vertex5 = L3+L4+1;
	vertex6 = L1+1;
	length = L1+L2+2;
	height = L3+L4+2;
	
	long	dx[6] = {-1, -1, 0, 1, 1, 0}, dy[6] = {0, 1, 1, 0, -1, -1};
	long 	x, y, xp, yp, dir, prob; 
	long	hitflag11, hitflag12, hitflag13, hitflag15, hitflag23, sideflag2, marker, runs;
	long    nohitcounter, hitcounter12not3, hitcounter13not2, hitcounter23not1, hitcounter123nodual;
	long	hitcounter12and13, hitcounter12and23, hitcounter13and23, hitcounter123dual123;
	FILE	*fprunfile, *fpoutfile;
	
	prob = (long) (2147483648.0*PROB);
	randinit(SEED);
	
	
	
	if (L4 > L1+L2) {
		printf("error40\n");
		exit(1);
	}
	
	if (L1 > L3+L4) {
		printf("error41\n");
		exit(1);
	}
	
	for (x = 0; x <= length; ++x)
		for (y = 0; y <= height; ++y)		
		{
			lat[x][y] = 0;
			
			if (x+y == vertex1)	lat[x][y] = BLOCK; //always blocked
			
			if (x+y == vertex1+1)	lat[x][y] = BLOCK; //always vacant

			if ((vertex1 <= x) && (y == 0))	lat[x][y] = BLOCK; //always blocked

			if ((y <= vertex3) && (x == length))	lat[x][y] = BLOCK; //always blocked

			if ((0 < y) && (y <= vertex3) && (x == length-1))	lat[x][y] = BLOCK; //always vacant

			if (x+y == length+L3+1)	lat[x][y] = BLOCK; //always blocked

			if ((x <= vertex4) && (y == height))	lat[x][y] = BLOCK; //alwyas blocked

			if ((0 < x) && (x <= vertex4) && (y == height-1))	lat[x][y] = BLOCK; //always vacant

			if ((vertex6 <= y) && (y <= height) && (x == 0))	lat[x][y] = BLOCK; //always blocked
		}
	
	nohitcounter = hitcounter12not3 = hitcounter13not2 = hitcounter23not1 = hitcounter123nodual = 0;
	hitcounter12and13 = hitcounter12and23 = hitcounter13and23 = hitcounter123dual123 = 0;
	
	for (runs = 1; runs <= RUNSMAX; ++runs)
	{
		//printf("%10ld\n", runs);
		hitflag11 = hitflag12 = hitflag13 = hitflag15 = hitflag23 = sideflag2 = 0;
		
		for (x = vertex1+1; x < length; ++x)	lat[x][1] = 2*runs-1;
		
		for (y = vertex6; y < vertex5; ++y)	lat[1][y] = BLOCK;
		
		for (x = vertex4; x < length-1; ++x){
			y = length+L3-x;
			lat[x][y] = 2*runs-1;
		}
		
		dir = 60;
		x = vertex1+1;
		y = 1;
		marker = vertex6;
		
		do
		{
			xp = x+dx[dir%6];
			yp = y+dy[dir%6];
			if (lat[xp][yp] >= 2*runs) ++dir;
			else if (lat[xp][yp] == 2*runs-1 || NewRandomInteger < prob){
				lat[xp][yp] = 2*runs-1;
				x = xp;
				y = yp;
				--dir;
			}
			else {
				lat[xp][yp] = 2*runs;
				++dir;
			}
			
			if (x == 2) {
				if (y >= vertex6 && y > marker) marker = y;
				if (hitflag11 + hitflag13 == 5){
					--hitflag11;
					--hitflag13;
				}
				if (hitflag11 == 3 && hitflag13 <= 1) --hitflag11;
				if (hitflag11 <= 1 && hitflag13 == 3) --hitflag13;
				hitflag12 = 3;
				sideflag2 = 1;
			}
			
			if (y == 1) {
				if (hitflag12 + hitflag13 == 5){
					--hitflag12;
					--hitflag13;
				}
				if (hitflag12 == 3 && hitflag13 <= 1) --hitflag12;
				if (hitflag12 <= 1 && hitflag13 == 3) --hitflag13;
				hitflag11 = 3;
			}
			
			if (y == height-2) hitflag15 = 1;
			
			if (x+y == length+L3-1) {
				if (hitflag11 + hitflag12 == 5){
					--hitflag11;
					--hitflag12;
				}
				if (hitflag11 == 3 && hitflag12 <= 1) --hitflag11;
				if (hitflag11 <= 1 && hitflag12 == 3) --hitflag12;
				hitflag13 = 3;
			}
			
			if (x+y == length+L3) break;
		}
		while ((y != 1) || (x != vertex2));
		
		for (y = vertex6; y < vertex5; ++y){
			lat[1][y] = 2*runs-1;
		}
		
		dir = 62;
		x = 1;
		y = vertex5-1;
		
		if (hitflag15 == 0 && y > marker+sideflag2)
			do
			{
				xp = x+dx[dir%6];
				yp = y+dy[dir%6];
				if (lat[xp][yp] >= 2*runs) ++dir;
				else if (lat[xp][yp] == 2*runs-1 || NewRandomInteger < prob){
					lat[xp][yp] = 2*runs-1;
					x = xp;
					y = yp;
					--dir;
				}
				else {
					lat[xp][yp] = 2*runs;
					++dir;
				}
				if (x+y == length+L3-1) {
					hitflag23 = 1;
					break;
				}
			}
		while ((x != 1) || (y != marker+sideflag2));
		
		if (hitflag12 == 2 && hitflag13 == 3) ++hitcounter123nodual;
		
		else if (hitflag12 == 0 && hitflag13 == 0 && hitflag23 == 0) ++nohitcounter;
		
		else if (hitflag11 == 0 && hitflag12 == 3 && hitflag13 == 0 && hitflag23 == 0) ++hitcounter12not3;

		else if (hitflag11 == 3 && hitflag12 == 2 && hitflag13 == 0 && hitflag23 == 0) ++hitcounter12not3;

		else if (hitflag11 == 2 && hitflag12 == 3 && hitflag13 == 0 && hitflag23 == 0) printf("error1%10ld\n", runs);
		
		else if (hitflag11 == 0 && hitflag12 == 0 && hitflag13 == 3 && hitflag23 == 0) ++hitcounter13not2;

		else if (hitflag11 == 2 && hitflag12 == 0 && hitflag13 == 3 && hitflag23 == 0) ++hitcounter13not2;

		else if (hitflag11 == 3 && hitflag12 == 0 && hitflag13 == 2 && hitflag23 == 0) printf("error2%10ld\n", runs);
		
		else if (hitflag11 == 0 && hitflag12 == 0 && hitflag13 == 0 && hitflag23 == 1) ++hitcounter23not1;

		else if (hitflag11 == 3 && hitflag12 == 0 && hitflag13 == 0 && hitflag23 == 1) ++hitcounter23not1;
		
		else if (hitflag11 == 0 && hitflag12 == 3 && hitflag13 == 0 && hitflag23 == 1) ++hitcounter12and23;

		else if (hitflag11 == 3 && hitflag12 == 2 && hitflag13 == 0 && hitflag23 == 1) ++hitcounter12and23;

		else if (hitflag11 == 2 && hitflag12 == 3 && hitflag13 == 0 && hitflag23 == 1) printf("error3%10ld\n", runs);
		
		else if (hitflag11 == 0 && hitflag12 == 0 && hitflag13 == 3 && hitflag23 == 1) ++hitcounter13and23;

		else if (hitflag11 == 2 && hitflag12 == 0 && hitflag13 == 3 && hitflag23 == 1) ++hitcounter13and23;

		else if (hitflag11 == 3 && hitflag12 == 0 && hitflag13 == 2 && hitflag23 == 1) printf("error4%10ld\n", runs);
		
		else if (hitflag11 == 2 && hitflag12 == 1 && hitflag13 == 3 && hitflag23 == 0) ++hitcounter12and13;
		
		else if (hitflag11 == 2 && hitflag12 == 1 && hitflag13 == 3 && hitflag23 == 1) ++hitcounter123dual123;
		
		else printf("error5%10ld\n", runs);
		
		if ((runs & PRINTFREQ2) == 0){//"bitwise and;" note that 2^p-1 & 2^p = 0.
			
			fprunfile = fopen(RUNFILE, "w");
			fprintf(fprunfile, "%15ld\n", runs);
			fclose(fprunfile);
			
			fpoutfile = fopen(OUTFILE, "w");
			fprintf(fpoutfile, "%f%15f\n", r, nohitcounter*1.0/runs);
			fprintf(fpoutfile, "%f%15f\n", r, hitcounter12not3*1.0/runs);
			fprintf(fpoutfile, "%f%15f\n", r, hitcounter13not2*1.0/runs);
			fprintf(fpoutfile, "%f%15f\n", r, hitcounter23not1*1.0/runs);
			fprintf(fpoutfile, "%f%15f\n", r, hitcounter12and13*1.0/runs);
			fprintf(fpoutfile, "%f%15f\n", r, hitcounter12and23*1.0/runs);
			fprintf(fpoutfile, "%f%15f\n", r, hitcounter13and23*1.0/runs);
			fprintf(fpoutfile, "%f%15f\n", r, hitcounter123nodual*1.0/runs);
			fprintf(fpoutfile, "%f%15f\n", r, hitcounter123dual123*1.0/runs);
			fclose(fpoutfile);
		}
		
	}// for runs
	
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


