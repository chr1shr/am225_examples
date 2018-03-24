#include "pss_read.hh"

#include <cstdio>
#include <cstdlib>

/** Sets up the perfect squared square grid, readin
 * \param[in] n the side length of the perfect squared square.
 * \param[in] filename the name of the file with the constituent squares. */
pss_read::pss_read(int n_,const char* filename) : n(n_), r(n-1), rr(r*r),
	c(new char[rr]) {
	
	// Clear the grid point status array
	for(int i=0;i<rr;i++) c[i]=0;

	// Open the file and check for an error
	FILE *fp=fopen(filename,"r");
	if(fp==NULL) {
		fputs("Can't open input file\n",stderr);
		exit(1);
	}

	// Loop over the squares in the file
	int j=1,q,k,l,s,x,y;
	while((q=fscanf(fp,"%d %d %d",&x,&y,&s))==3) {

		// For each square, label the interior gridpoints in the square with a
		// unique ID. If the square has side length s, then it has (s-1)*(s-1)
		// interior gridpoints.
		for(l=y;l<y+s-1;l++) for(k=x;k<x+s-1;k++)
			c[k+l*r]=j;
		j++;

		// Record the positions of the squares. This information is not used in
		// this example, but would be useful in a complete Schur solution routine.
		xs.push_back(x);ys.push_back(y);
		ss.push_back(s);
	}

	// Print an error message if there was an error in reading the file
	if(q!=EOF) {
		fputs("Error reading file\n",stderr);
		exit(1);
	}
	fclose(fp);
}

/** The class destructor frees the dynamically allocated memory. */
pss_read::~pss_read() {
	delete [] c;
}

/** Prints an ASCII art image of the grid, with gridpoints in each component
 * square being labeled by different letters, and the glueing region being
 * labeled with dots. */ 
void pss_read::ascii_art() {
	int k,l;
	char *z=new char[n];z[n-1]=0;

	// Loop over each row of the grid. Since lines are printed from top to
	// bottom, but the y coordinate goes from bottom to top, the vertical index
	// is considered in reverse order.
	for(l=r-1;l>=0;l--) {

		// Assemble an ASCII art string for this row, and then print it 
		for(k=0;k<r;k++) z[k]=c[k+l*r]==0?'.':64+c[k+l*r];
		puts(z);
	}
	delete [] z;
}

/** Scans the grid and prints out the coordinates of each glueing point. */
void pss_read::analyze_glue() {
	int k,l,tot=0;char *cp=c;

	// Loop over the grid and find the zero entries corresponding to glueing
	// gridpoints
	for(l=0;l<r;l++) for(k=0;k<r;k++,cp++) if(*cp==0) {
		tot++;
		
		// For each glueing gridpoint, print the statuses of the four
		// orthogonal neighbors
		printf("(%d,%d):",k,l);
		if(l>0) printf(" D->%d",int(cp[-r]));
		if(k>0) printf(" L->%d",int(cp[-1]));
		if(k<r-1) printf(" R->%d",int(cp[1]));
		if(l<r-1) printf(" U->%d",int(cp[r]));
		putchar('\n');
	}
	printf("%d total glueing gridpoints\n",tot);
}
