#include <stdlib.h>
#include <stdio.h>
#include <mkl_lapacke.h>
#include <math.h>
#include <time.h>


struct matrix
{
    long number_of_rows;
    long number_of_columns;
    long upper_bands;
    long lower_bands;
    long number_of_rows_inv;
    double *m;
    double *m_inv;
    int *infoinv;
};

int intialise(struct matrix *ma, long lb , long ub, long nc)
{
    ma -> number_of_rows = 1+ub+lb;
    ma-> number_of_columns = nc;
    ma -> upper_bands = ub;
    ma -> lower_bands = lb;
    ma -> number_of_rows_inv = ma->upper_bands*2 + ma -> lower_bands +1;
    ma -> m = (double *)calloc(ma->number_of_rows * ma->number_of_columns,sizeof(double));
    ma -> m_inv = (double *) calloc((ma -> number_of_rows + ma -> lower_bands)*(ma -> number_of_columns),sizeof(double));
    ma -> infoinv = (int *) malloc(sizeof(int)*ma->number_of_columns);

    if (!ma->m || !ma->m_inv || !ma->infoinv) {
        free(ma->m); // Safe to call free on NULL
        free(ma->m_inv);
        free(ma->infoinv);
        return 0;
    }

    /*
    for(long c = 0; c < ma->number_of_rows * ma->number_of_columns; c++ )
    {
        ma -> m[c] = 0.0;
    }
    */
    return 1;
};

void Free_m(struct matrix *ma)
{
    free(ma -> m);
    free(ma -> m_inv);
    free(ma -> infoinv);
}

void MemSwap(double** a, double** b){
  //function to swap two double arrays
  double* temp;
  temp = *a;
  *a = *b;
  *b = temp;
};


double *getp(struct matrix *ma, long r, long c)
{
    int bandnumber = ma -> upper_bands +r - c;
    if (r < 0 || c < 0 || r >= ma -> number_of_columns || c >= ma -> number_of_columns )
    {
        printf("index out of bound in getp: %ld %ld %ld \n", r,c,ma -> number_of_columns);
        exit(1);

    }

    return &ma -> m[ma -> number_of_rows*c + bandnumber];

}
double getv(struct matrix *ma, long r, long c)
{
    return *getp(ma,r,c);
}

double setv(struct matrix *ma, long r, long c, double value)
{
    *getp(ma,r,c) = value;
    return value;
}
 int solve_eq(struct matrix *ma, double *x, double *b)
 {
     int i,band_number;
     for(i =0; i < ma -> number_of_columns; i++)
    {
        for(band_number = 0; band_number < ma -> number_of_rows; band_number++)
        {
            ma ->m_inv[ma -> number_of_rows_inv*i + (band_number+ma->lower_bands)]=ma->m[ma->number_of_rows*i+band_number];
        }
        x[i] = b[i];
    }
    long ldab = ma ->lower_bands*2 + ma -> upper_bands +1;
    int info = LAPACKE_dgbsv(LAPACK_COL_MAJOR,ma->number_of_columns,ma -> lower_bands, ma->upper_bands,1,ma->m_inv,ldab,ma->infoinv, x, ma->number_of_columns);
    return info;
 }
 int print_matrix(struct matrix *ma)
 {
    long i,j;
    for(i = 0; i < ma -> number_of_columns;i++)
    {
        for(j = 0; j < ma -> number_of_rows;j++)
        {
            long row = i+j-ma->upper_bands;
            if (row <0 ){row = 0;}
            printf("%ld %ld || %ld %ld | %lg \n ",row,i,j,i,ma->m[ma ->number_of_rows*i+j]);
        }
    }
    return 0;
 }

int main()
{
    double x_l;
    double a_l;
    double b_l;
    double x_r;
    double a_r;
    double b_r;
    double E_0;
    long int n;
    long int I;
    struct matrix ma;
    long i,j;
    FILE *inp;
    inp = fopen("input.txt","r");
    int p = fscanf(inp,"%lf %lf %lu %lf %lf %lf %lf %lf %lu", &x_l,&x_r,&n,&a_l,&a_r,&b_l,&b_r,&E_0,&I);
    if(p != 9){
        fprintf(stderr, "Error in reading parameters");
        fclose(inp);
        exit(1);
    }
    fclose(inp);
    double diff_x = (x_r - x_l)/(n-1);

    FILE *pot;
    pot = fopen("potential.txt","r");
    double *l = malloc(sizeof(double)*n);


    if (pot == NULL) {
      fprintf(stderr, "Error in opening potential.txt\n");
      exit(1);
    }
    if (l == NULL) {
      fprintf(stderr, "Memory allocation failed for potentials\n");
      fclose(pot);
      exit(1);
    }
    for (i = 0; i<n; i++)
    {
        int r = fscanf(pot,"%lf", &l[i]);
        if(r != 1)
        {
            fprintf(stderr, "Error in reading potentials");
            free(l);
            fclose(pot);
            exit(1);
        }
    }
    fclose(pot);

    long ncolumns = n;

    long lower_bands = 1;
    long upper_bands = 1;

    //intialise(&ma,lower_bands,upper_bands,ncolumns);

    if (!intialise(&ma,lower_bands,upper_bands,ncolumns)) {
        fprintf(stderr, "Failed to initialize the matrix\n");
        free(l);
        return 1;
    }
    double *x = (double *)calloc(ncolumns,sizeof(double));
    double *b = (double *)calloc(ncolumns,sizeof(double));
    if (!x || !b) {
        fprintf(stderr, "Memory allocation failed for x or b\n");
        Free_m(&ma);
        free(l);
        free(x);
        free(b);
        return 1;
    }



    srand(0);
    for(i = 0; i < ncolumns; i++)
    {
       if (i > 0){setv(&ma, i, i-1,-1/(diff_x *diff_x));}
       setv(&ma,i,i,((2+l[i]*diff_x*diff_x)/(diff_x*diff_x))-E_0);
       if(i < ncolumns -1){setv(&ma ,i,i+1,-1/(diff_x*diff_x));}
       b[i] = rand();
    }
    if (b_l == 0){
        setv(&ma,0,0,1);
        setv(&ma,0,1,0);
        b[0] = 0;
    }else
    {
        setv(&ma,0,0, ((-(2*a_l*diff_x/b_l) + (2+l[0]*diff_x*diff_x))/(diff_x*diff_x)) - E_0);
        setv(&ma,0,1,-2/(diff_x*diff_x));
    }
    if(b_r == 0)
    {
        setv(&ma,ncolumns-1,ncolumns-1,1);
        setv(&ma,ncolumns-1,ncolumns-2,0);
        b[ncolumns-1] = 0;
    }else
    {
        setv(&ma,n-1,n-2,-2/(diff_x*diff_x));
        setv(&ma,n-1,n-1, ((-(2*a_r*diff_x/b_r) + (2+l[n-1]*diff_x*diff_x))/(diff_x*diff_x)) - E_0);
    }



    FILE *f = fopen("output.txt","w");
    if(f == NULL)
    {
        printf("Error opening file");
        exit(1);
    }
    for (i = 0; i < ncolumns; i++){fprintf(f,"%ld %.6f %.6f \n", 0,x_l + i*(diff_x),b[i]);}
    fprintf(f,"\n");
    int info;
    for(i=1;i<=I;i++){
        info = solve_eq(&ma,x,b);
        if(info == 0){
          MemSwap(&x,&b);
          if (b_l == 0.0){
            b[0] = 0.0;
            }

          if (b_r == 0.0){
            b[n-1] = 0.0;
          }

          for (j = 0; j < n; j++) {
              fprintf(f, "%ld, %.21f, %.21f \n", i,x_l + j * diff_x, b[j]);
          }
          //fprintf(f,"\n");
        }else{
          printf("Error Inverting matrix");
          fclose(f);
          Free_m(&ma);
          free(l);
          free(x);
          free(b);
          exit(4);
        }
    }
    fclose(f);

    Free_m(&ma);
    free(l);
    free(x);
    free(b);

    return 0;
}


