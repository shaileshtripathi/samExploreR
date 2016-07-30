/***************************************************************

   The Subread and Rsubread software packages are free
   software packages:
 
   you can redistribute it and/or modify it under the terms
   of the GNU General Public License as published by the 
   Free Software Foundation, either version 3 of the License,
   or (at your option) any later version.

   Subread is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty
   of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
   
   See the GNU General Public License for more details.

   Authors: Drs Yang Liao and Wei Shi

  ***************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <Rversion.h>

#ifndef WIN32
#define WIN32 _WIN32
#endif

#if (R_VERSION >= R_Version(2,3,0))
#define R_INTERFACE_PTRS 1
#define CSTACK_DEFNS 1
#include <Rinterface.h>
#endif 



int main_junction(int argc,char ** argv);
int main_align(int argc,char ** argv);
int main_buildindex(int argc,char ** argv);
int sam2bed(int argc,char *argv[]);
int propmapped(int argc,char *argv[]);
int readSummary(int argc,char *argv[]);
int main_snp_calling_test(int argc,char *argv[]);
int main_repeated_test(int argc, char *argv[]);
int main_qualityScores(int argc, char *argv[]);
int findCommonVariants(int argc, char *argv[]);

void R_mergeVCF(int * nargs, char ** argv)
{
	char * r_argv, ** c_argv;
	int i,n;

	n = *nargs;
	r_argv = (char *)calloc(15000, sizeof(char));
	strcpy(r_argv,*argv);

//	printf("N=%d; V=%s\n", n, r_argv);
	c_argv = (char **) calloc(n+1,sizeof(char *));
	for(i=0;i<1+n;i++) c_argv[i] = (char *)calloc(300,sizeof(char));
	strcpy(c_argv[0],"R_mergeVCF");
	strcpy(c_argv[1],strtok(r_argv,";"));
	for(i=2;i<n+1;i++) strcpy(c_argv[i],strtok(NULL,";"));

	findCommonVariants(n+1, c_argv);

	free(r_argv);
	for(i=0;i<n+1;i++) free(c_argv[i]);
	free(c_argv);
}

void R_buildindex_wrapper(int * nargs, char ** argv)
{

	char * r_argv, ** c_argv;
	int i,n;
	
	r_argv = (char *)calloc(1000, sizeof(char));
	strcpy(r_argv,*argv);
	
	n = *nargs;
	c_argv = (char **) calloc(n,sizeof(char *));
	for(i=0;i<n;i++) c_argv[i] = (char *)calloc(200,sizeof(char));
	strcpy(c_argv[0],strtok(r_argv,","));
	for(i=1;i<n;i++) strcpy(c_argv[i],strtok(NULL,","));

	main_buildindex(n,c_argv);

	for(i=0;i<n;i++) free(c_argv[i]);
	free(c_argv);
	free(r_argv);

}

void R_align_wrapper(int * nargs, char ** argv)
{
	uintptr_t old_cstack_limit = R_CStackLimit;
	R_CStackLimit =(uintptr_t)-1;
	
        char * r_argv, ** c_argv;
        int i,n;
    
        r_argv = (char *)calloc(1000, sizeof(char));
        strcpy(r_argv,*argv);
    	
        n = *nargs;
        c_argv = (char **) calloc(n,sizeof(char *));
        for(i=0;i<n;i++) c_argv[i] = (char *)calloc(200,sizeof(char));
        strcpy(c_argv[0],strtok(r_argv,","));
        for(i=1;i<n;i++) strcpy(c_argv[i],strtok(NULL,","));

        main_align(n,c_argv);

        for(i=0;i<n;i++) free(c_argv[i]);
        free(c_argv);
        free(r_argv);

	R_CStackLimit = old_cstack_limit;
}

void R_junction_wrapper(int * nargs, char ** argv)
{
	uintptr_t old_cstack_limit = R_CStackLimit;
	R_CStackLimit =(uintptr_t)-1;

        char * r_argv, ** c_argv;
        int i,n;
    
        r_argv = (char *)calloc(1000, sizeof(char));
        strcpy(r_argv,*argv);
    
        n = *nargs;
        c_argv = (char **) calloc(n,sizeof(char *));
        for(i=0;i<n;i++) c_argv[i] = (char *)calloc(200,sizeof(char));
        strcpy(c_argv[0],strtok(r_argv,","));
        for(i=1;i<n;i++) strcpy(c_argv[i],strtok(NULL,","));

        main_junction(n,c_argv);

        for(i=0;i<n;i++) free(c_argv[i]);
        free(c_argv);
        free(r_argv);

	R_CStackLimit = old_cstack_limit;
}


void R_sam2bed_wrapper(int * nargs, char ** argv)
{

        char * r_argv, ** c_argv;
        int i,n;

        r_argv = (char *)calloc(1000, sizeof(char));
        strcpy(r_argv,*argv);

        n = *nargs;
        c_argv = (char **) calloc(n,sizeof(char *));
        for(i=0;i<n;i++) c_argv[i] = (char *)calloc(200,sizeof(char));
        strcpy(c_argv[0],strtok(r_argv,","));
        for(i=1;i<n;i++) strcpy(c_argv[i],strtok(NULL,","));

        sam2bed(n,c_argv);

        for(i=0;i<n;i++) free(c_argv[i]);
        free(c_argv);
        free(r_argv);

}


void R_propmapped_wrapper(int * nargs, char ** argv)
{

        char * r_argv, ** c_argv;
        int i,n;

        r_argv = (char *)calloc(1000, sizeof(char));
        strcpy(r_argv,*argv);

        n = *nargs;
        c_argv = (char **) calloc(n,sizeof(char *));
        for(i=0;i<n;i++) c_argv[i] = (char *)calloc(200,sizeof(char));
        strcpy(c_argv[0],strtok(r_argv,","));
        for(i=1;i<n;i++) strcpy(c_argv[i],strtok(NULL,","));

        propmapped(n,c_argv);

        for(i=0;i<n;i++) free(c_argv[i]);
        free(c_argv);
        free(r_argv);

}

void R_readSummary_wrapper(int * nargs, char ** argv)
{
	uintptr_t old_cstack_limit = R_CStackLimit;
	R_CStackLimit =(uintptr_t)-1;
	//printf("RCL=%ld\n", R_CStackLimit);

        char * r_argv, ** c_argv;
        int i,n, arg_len;

	arg_len = strlen(*argv);
        r_argv = (char *)calloc(1+arg_len, sizeof(char));
        strcpy(r_argv,*argv);

        n = *nargs;
        c_argv = (char **) calloc(n,sizeof(char *));
	if(strstr(r_argv, ",,")==NULL )
	{
		for(i=0; i<n; i++)
		{
			char * current_arg = strtok(i?NULL:r_argv,",");
			if(current_arg == NULL)
				break;
			arg_len = strlen(current_arg);
			c_argv[i] = (char *)calloc(1+arg_len,sizeof(char));
			strcpy(c_argv[i], current_arg);
			//modified code
			//printf("\n%s\n", current_arg);
			// end modified code
		}

		n=i;

		readSummary(n,c_argv);

		for(i=0;i<n;i++) free(c_argv[i]);
	}
	else Rprintf("No input files are provided. \n");
        free(c_argv);
        free(r_argv);

	R_CStackLimit = old_cstack_limit;
}


void R_SNPcalling_wrapper(int * nargs, char ** argv)
{
	uintptr_t old_cstack_limit = R_CStackLimit;
	R_CStackLimit =(uintptr_t)-1;

        char * r_argv, ** c_argv;
        int i,n;

        r_argv = (char *)calloc(1000, sizeof(char));
        strcpy(r_argv,*argv);

        n = *nargs;
        c_argv = (char **) calloc(n,sizeof(char *));
        for(i=0;i<n;i++) c_argv[i] = (char *)calloc(200,sizeof(char));
        strcpy(c_argv[0],strtok(r_argv,","));
        for(i=1;i<n;i++) strcpy(c_argv[i],strtok(NULL,","));

        main_snp_calling_test(n,c_argv);

        for(i=0;i<n;i++) free(c_argv[i]);
        free(c_argv);
        free(r_argv);

	R_CStackLimit = old_cstack_limit;
}


void R_removeDupReads_wrapper(int * nargs, char ** argv)
{

        char * r_argv, ** c_argv;
        int i,n;

        r_argv = (char *)calloc(1000, sizeof(char));
        strcpy(r_argv,*argv);

        n = *nargs;
        c_argv = (char **) calloc(n,sizeof(char *));
        for(i=0;i<n;i++) c_argv[i] = (char *)calloc(200,sizeof(char));
        strcpy(c_argv[0],strtok(r_argv,","));
        for(i=1;i<n;i++) strcpy(c_argv[i],strtok(NULL,","));

        main_repeated_test(n,c_argv);

        for(i=0;i<n;i++) free(c_argv[i]);
        free(c_argv);
        free(r_argv);

}


void R_qualityScores_wrapper(int * nargs, char ** argv)
{

	char * r_argv, ** c_argv;
	int i,n;
	
	r_argv = (char *)calloc(1000, sizeof(char));
	strcpy(r_argv,*argv);
	
	n = *nargs;
	c_argv = (char **) calloc(n,sizeof(char *));
	for(i=0;i<n;i++) c_argv[i] = (char *)calloc(200,sizeof(char));
	strcpy(c_argv[0],strtok(r_argv,","));
	for(i=1;i<n;i++) strcpy(c_argv[i],strtok(NULL,","));

	main_qualityScores(n,c_argv);

	for(i=0;i<n;i++) free(c_argv[i]);
	free(c_argv);
	free(r_argv);

}


