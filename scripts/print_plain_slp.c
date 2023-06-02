/*

Repair -- an implementation of Larsson and Moffat's compression and
decompression algorithms.
Copyright (C) 2010-current_year Gonzalo Navarro

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

Author's contact: Gonzalo Navarro, Dept. of Computer Science, University of
Chile. Blanco Encalada 2120, Santiago, Chile. gnavarro@dcc.uchile.cl

Modified by Gogis to output the SLP in a plain format

*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

typedef struct
  { unsigned int left,right;
  } Tpair;

long u; // |text| and later current |C| with gaps


unsigned int alph; // size of terminal alphabet, or smallest non terminal symbol

Tpair *R; // rules

size_t n; // |R|
// int n; // |R|

char *ff;
FILE *f;

void print_rule(unsigned int i, FILE *f) {
    fprintf(f, "%u %u\n", R[i-alph].left, R[i-alph].right);
}

int main (int argc, char **argv)

   { char fname[1024]; char outname[1024];
     //char *text;
     size_t len;
     FILE *Tf,*Rf,*Cf;
     struct stat s;
     fputs("==== Command line:\n",stderr);
     for(int i=0;i<argc;i++)
       fprintf(stderr," %s",argv[i]);
     fputs("\n",stderr);     
     if (argc != 2)
  { fprintf (stderr,"Usage: %s <filename>\n"
        "Generate <filename>.plainslp from its .C and .R "
        "extensions.\n"
        "This is a version for prefix-free parsing\n",argv[0]);
    exit(1);
  }
  
  // read .R file, store data in alpha and R[]
     strcpy(fname,argv[1]);
     strcat(fname,".R");
     if (stat (fname,&s) != 0)
  { fprintf (stderr,"Error: cannot stat file %s\n",fname);
    exit(1);
  }
     len = s.st_size;
     Rf = fopen (fname,"r");
     if (Rf == NULL)
  { fprintf (stderr,"Error: cannot open file %s for reading\n",fname);
    exit(1);
  }
     if (fread(&alph,sizeof(int),1,Rf) != 1)
  { fprintf (stderr,"Error: cannot read file %s\n",fname);
    exit(1);
  }
  // note that in the original char-based repair the R file contains also
  // a map between the 0...alph-1 and the actual symbols in the input file
  // here alph is 256 and there is no such map (this is why the two .R .C
  // formats are not compatible).
  
     // n is the number of rules, sizeof(int) accounts for alpha
     n = (len-sizeof(int))/sizeof(Tpair);
     // allocate and read array of rules stored as pairs 
     R = (void*)malloc(n*sizeof(Tpair));
     if (fread(R,sizeof(Tpair),n,Rf) != n)
  { fprintf (stderr,"Error: cannot read file %s\n",fname);
    exit(1);
  }
     fclose(Rf);

   // open C file and get the nuber of symbols in it 
     strcpy(fname,argv[1]);
     strcat(fname,".C");
     if (stat (fname,&s) != 0)
  { fprintf (stderr,"Error: cannot stat file %s\n",fname);
    exit(1);
  }
     Cf = fopen (fname,"r");
     if (Cf == NULL)
  { fprintf (stderr,"Error: cannot open file %s for reading\n",fname);
    exit(1);
  }

  // open output file
     strcpy(outname,argv[1]);
     strcat(outname,".plainslp");
     Tf = fopen (outname,"w");
     if (Tf == NULL)
  { fprintf (stderr,"Error: cannot open file %s for writing\n",outname);
    exit(1);
  }

  for (int i = 0; i < n; i++) {
    print_rule(alph+i, Tf);
  }
  fclose(Tf);
 return 0;
}

