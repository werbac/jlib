#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*
  gcc fread.c -shared -fPIC -o fread.so
  push!(Base.DL_LOAD_PATH, pwd())
  x=ccall((:xmean,"/mapr/mapr04p/analytics0001/analytic_users/jmdl/Hormel_Rev_8_967/fread"),Float64,(Float64,Float64),2.0,5.0)
  
  x=ccall((:myFunction,"/mapr/mapr04p/analytics0001/analytic_users/jmdl/Hormel_Rev_8_967/fread"),String)
  
  push!(Base.DL_LOAD_PATH, pwd())
  x=ccall((:myFunction,"fread"), Cstring, (Cstring,), "tst.csv")
*/



/*const char* xfread()  / * http://www.fundza.com/c4serious/fileIO_reading_all/  */
/*   http://docs.julialang.org/en/stable/manual/calling-c-and-fortran-code/
 * Julia Code *
  x=ccall((:xfread,"fread"), Cstring, (Cstring,), "tst.csv")
  unsafe_string(x) 

  fread(fname::String)=unsafe_string(ccall((:xfread,"fread"), Cstring, (Cstring,), fname))
*/

const char* xfread(char* fname)
{
    FILE    *infile; /* declare a file pointer */
    char    *buffer;
    long    numbytes;
    infile = fopen(fname, "r");   /* open an existing file for reading */
    if(infile == NULL)  /* quit if the file does not exist */
        return "NULL"; //return 1;
    fseek(infile, 0L, SEEK_END); /* Get the number of bytes */
    numbytes = ftell(infile);
    printf("The file called orig.csv s size %ld\n", numbytes);

    fseek(infile, 0L, SEEK_SET);	/* reset the file position indicator to * the beginning of the file */
    buffer = (char*)calloc(numbytes, sizeof(char));	 /* grab sufficient memory for the  buffer to hold the text */
    if(buffer == NULL)  /* memory error */
        return "NULL"; //return 1;
    fread(buffer, sizeof(char), numbytes, infile);  /* copy all the text into the buffer */
    fclose(infile);
    //printf("The file called orig.csv contains this text\n\n%s", buffer); /* confirm we have read the file by outputing it to the console */
    //free(buffer);  /* free the memory we used for the buffer */ /*xxzreturn buffer;*/
    //return 0;
    return buffer;
}


double xmean(double a, double b) 
    { return (a+b) / 2;
    }



const char* myFunction(char* fname)
{  printf("working on %s\n\n",fname);
    return "done with My String";
}

//gcc fread.c -shared -fPIC -o fread.so; gcc fread.c; ./a.out
// x=unsafe_string(ccall((:xfread,"fread"), Cstring, (Cstring,), "tst.csv"))
// dfd = readtable(IOBuffer( x ) )
void main()
{
    double  numbytes;
    numbytes =  xmean(2.0,5.0);
    printf("xmean is :  %d\n", numbytes);
    
    printf("finished and returned %s\n\n", myFunction("tst.csv") )   ; 
    //const char* szSomeString = xfread(); // fraught with problems
    //printf("%s", szSomeString);    
    printf("file contents :  \n%s\n\n", xfread("tst.csv"));
}



/*
using DataFrames
push!(Base.DL_LOAD_PATH, pwd())
root=pwd();
function loadDF()
    x=unsafe_string(ccall((:xfread,"fread"), Cstring, (Cstring,), "orig.csv"));
    dfd = readtable(IOBuffer(x))
    #df_data = readtable(root*"/orig.csv",header=false);
    df_h = readtable(root*"/origHead.csv",header=false);  
    names!(dfd, convert(Array{Symbol}, df_h[:x1]) )
end
@time df_in=loadDF();

*/

