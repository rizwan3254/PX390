#include <stdio.h>
#include <stdlib.h>
#include <mkl_lapacke.h>
#include <time.h>
#include <math.h>
#include <stdbool.h>

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

 int printfull_matrix(struct matrix *ma)
 {
    int i,j;
    for(i=0; i<ma -> number_of_columns;i++) {
        for(j=0; j<ma -> number_of_columns; j++) {
            // printf("%ld %ld %g \n",i,j,getv(ma, i, j));
            printf("%g\t",getv(ma, i, j));
    }
    printf("\n");
  }
  return 0;
 }
int printV(FILE *f,double **v,int NumActive, int **position, int Numy,double t) {
  if(f == NULL){
    printf("Error in output file\n");
    return 1;
  }else{
    for(int i = 0; i < NumActive; i++) {
      fprintf(f, "%.6f, %d, %d, %.6f\n", t, (*position)[i]/Numy, (*position)[i]%Numy , (*v)[i]);
    }
    fprintf(f, "\n");
    return 0;
  }
}

int main()
{
    int i;
    int j;

    int numberX;
    int numberY;
    int numberActive;
    double LengthX;
    double LengthY;
    double final_time;
    double lambda;
    double diag_timestep;

    FILE *inp;
    inp = fopen("input.txt","r");
    if(inp == NULL){
        printf("Error opening the input file\n");
        exit(1);
    }
    else{
        int p = fscanf(inp,"%d %d %d %lf %lf %lf %lf %lf", &numberX,&numberY,&numberActive,&LengthX,&LengthY,&final_time,&lambda,&diag_timestep);
        if(p != 8){
            fprintf(stderr, "Error in reading parameters");
            fclose(inp);
            exit(2);
        }
        fclose(inp);

    }
    double diff_x = LengthX/numberX;
    double diff_y = LengthY/numberY;

    int x;
    int y;
    double *v = malloc(sizeof(double)*numberActive);
    int *gridI = malloc(sizeof(int)*numberActive);
    int *vI = malloc(sizeof(int)*numberX*numberY);

    for(i=0;i<numberX*numberY;i++){
        vI[i] = -1;
    }

    FILE *coeff;
    coeff = fopen("coefficients.txt","r");
    double max_c = sqrt(fabs(lambda));

    if (coeff == NULL) {
      fprintf(stderr, "Error in opening coefficients.txt\n");
      free(v);
      free(gridI);
      free(vI);
      exit(1);
    }
    else
    {
        for (i = 0; i < numberActive; i++)
        {
            int r = fscanf(coeff,"%d %d %lf",&x, &y, &v[i]);
            if(r != 3)
            {
                fprintf(stderr, "Error reading from coefficients file\n");
                fclose(coeff);
                free(v);
                free(gridI);
                exit(2);
            }
            else
            {

                if(fabs(v[i]) > max_c){
                    max_c = fabs(v[i]);
                }
                gridI[i] = x + y*numberY;
                vI[numberY*y + x] = i;

            }
        }
        fclose(coeff);
    }
    double bound1 = ((0.5*diff_x*diff_x*diff_y*diff_y)/( diff_x*diff_x + diff_y*diff_y));
    double bound2;
    bound2 = 1/(max_c*(sqrt(fabs(lambda))+max_c));

    double minim = fmin(bound1,bound2);

    int f = 1;
    double diff_t = diag_timestep;
    while(diff_t > minim){
        f = f + 1;
        diff_t = diag_timestep/f;
    }







    struct matrix m;
    intialise(&m, numberY, numberY, numberActive);
    bool t;
    bool b;
    bool l;
    bool r;

    for(i=0; i < numberActive; i++)
    {
        t = true;
        b = true;
        l = true;
        r = true;

        if(gridI[i]%numberY == numberY - 1){
            t = false;
        }else if (i <(numberActive-1) && gridI[i+1] - gridI[i] != 1 ){
            t = false;
        }else{
            t = true;
        }

        if(gridI[i]%numberY == 0){
            b = false;
        }else if((i > 0) && gridI[i] - gridI[i-1] != 1){
            b = false;
        }else{
            b = true;
        }

        if(gridI[i] < numberY){
            l = false;
        }else if(vI[gridI[i] - numberY] == -1){
            l = false;
        }else{
            l = true;
        }

        if(gridI[i] >= (numberX-1)*numberY  ){
            r = false;
        }else if(vI[gridI[i] + numberY] == -1){
            r = false;
        }else{
            r = true;
        }
        setv(&m, i, i, 1);
        if(b){
            setv(&m, i , i - 1  , -(diff_t)/(diff_y*diff_y));
            setv(&m, i , i  , getv(&m,i,i) +(diff_t)/(diff_y*diff_y));
        }
        if(t){
            setv(&m, i , i + 1  , -(diff_t)/(diff_y*diff_y));
            setv(&m, i , i  , getv(&m,i,i) +(diff_t)/(diff_y*diff_y));
        }
        if(r){
            setv(&m, i , vI[gridI[i] + numberY]  , -(diff_t)/(diff_x*diff_x));
            setv(&m, i , i  , getv(&m,i,i) +(diff_t)/(diff_x*diff_x));
        }
        if(l){
            setv(&m, i , vI[gridI[i] - numberY]  , -(diff_t)/(diff_x*diff_x));
            setv(&m, i , i  , getv(&m,i,i) +(diff_t)/(diff_x*diff_x));
        }

    }

    free(vI);

    FILE *f3 = fopen("output.txt", "w");
    if(printV(f3, &v, numberActive, &gridI,numberY,0.0) == 1){
        fprintf(stderr, "Error printing to output file\n");
        Free_m(&m);
        free(gridI);
        free(vI);
        free(v);
        exit(3);
    }

    int step = floor(final_time/diff_t);
    double *nV = malloc(sizeof(double)*numberActive);

    //int content;
    //double t;
    for(i = 1;i<=step;i++){
        for(j = 0; j < numberActive; j++) {
            v[j] = v[j] *((1 + diff_t*lambda) - diff_t * v[j]* v[j]);
    }
    solve_eq(&m,nV,v);
    MemSwap(&v,&nV);
    if(printV(f3, &v, numberActive, &gridI,numberY,i*diff_t) == 1){
      fprintf(stderr, "Error printing to output file\n");
      fclose(f3);
      Free_m(&m);
      free(gridI);
      free(v);
      //free(vI);
      free(nV);
      exit(3);
    }
  }

  fclose(f3);
  Free_m(&m);
  free(gridI);
  free(v);
  free(nV);
  //free(vI);

}
