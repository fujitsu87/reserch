void vec_2_mat()
{
    int i,j,t;
    for(i=0;i<K/DK;++i)
    {
        for(j=0;j<K/DK;++j)
        {
            gsl_complex z;
            z = gsl_vector_complex_get(g,i*K/DK+j);
            gsl_matrix_complex_set(G,j,i,z);
        }
    }
}
void culc_z(int t)
{
    int i;
    gsl_vector_complex* tmp;
    tmp = gsl_vector_complex_alloc(K/DK);
    for(i=0;i<K/DK;i++)
    {
        gsl_vector_complex_set(tmp,i,gsl_matrix_complex_get(y_sub,i,t));
    }
    AHx(G,tmp,Z);
}

void culc_x_h_sub(int t)
{
    int i;
    double tmp = sqrt(Pk/(2*rho_k));
    gsl_complex z;
    for(i=0;i<K/DK;i++){
        if(gsl_vector_complex_get(Z,i).dat[0] > 0)z.dat[0] = tmp;
        else z.dat[0] = (-1) * tmp;
        if(gsl_vector_complex_get(Z,i).dat[1] > 0)z.dat[1] = tmp;
        else z.dat[1] = (-1) * tmp;
        gsl_matrix_complex_set(x_h_sub,i,t,z);
    }
}


void culc_inno(int t)
{
    int i;
    gsl_matrix_complex* I_K;
    gsl_vector_complex *tmp;
    gsl_complex z;

    I_K = gsl_matrix_complex_calloc(K/DK,K/DK);
    tmp = gsl_vector_complex_calloc(K/DK);

    for(i=0;i<K/DK;++i){
        z = gsl_matrix_complex_get(y_sub,i,t);
        gsl_vector_complex_set(inno,i,z);
    }
    ATx(x_kro,g,tmp);
    gsl_vector_complex_sub(inno,tmp);

}
gsl_matrix_complex *ConjuMatrix(gsl_matrix_complex *a) {
    int i, j;
    int m, p;
    m = a->size1;
    p = a->size2;

    gsl_matrix_complex *c = gsl_matrix_complex_alloc(m, p);
    gsl_complex z;

     for (i = 0; i < m; i++)    {
          for (j = 0; j < p; j++)   {
            z = gsl_matrix_complex_get (a, i, j);
            z = gsl_complex_conjugate(z);
            gsl_matrix_complex_set (c,i,j,z);
          }
      }

    return c;
}

void culc_x_kro(int t)
{
    int i;
    
    gsl_matrix_complex *x_tmp;
    gsl_matrix_complex* I_K;
    gsl_complex z;
    x_tmp = gsl_matrix_complex_calloc(K/DK,1);
    I_K = gsl_matrix_complex_calloc(K/DK,K/DK);
    gsl_matrix_complex_set_identity(I_K);

    for(i=0;i<K/DK;++i){
        z = gsl_matrix_complex_get(x_h_sub,i,t);
        gsl_matrix_complex_set(x_tmp,i,0,z);
    }

    x_kro = KPro(x_tmp,I_K);
    x_kro_conju = ConjuMatrix(x_kro);
}

void culc_R_a(int t)
{
    int i;
    gsl_matrix_complex* I_K;
    gsl_matrix_complex *tmp;
    gsl_complex z;

    I_K = gsl_matrix_complex_calloc(K/DK,K/DK);
    tmp = gsl_matrix_complex_calloc((K/DK),(K/DK)*(K/DK));


    gsl_matrix_complex_set_identity(I_K);

    ATB(x_kro,P,tmp);
    AB(tmp,x_kro_conju,R_a);

    tmp = GSLMatrixFree(tmp);
    tmp = gsl_matrix_complex_calloc(K/DK,K/DK);

    z.dat[0] = N0;z.dat[1] = 0;
    for(i=0;i<DK-1;i++)
    {
        AB(D_other,R_x,tmp);
        z = gsl_complex_add(z,trace(tmp));
    }

    gsl_matrix_complex_memcpy(tmp,I_K);
    gsl_matrix_complex_scale(tmp,z);
    gsl_matrix_complex_add(R_a,tmp);
}
void culc_g()
{
    gsl_matrix_complex* tmp,* tmp2,*R_rev;
    gsl_vector_complex *right;
    tmp = gsl_matrix_complex_calloc((K/DK*K/DK),K/DK);
    tmp2 = gsl_matrix_complex_calloc((K/DK*K/DK),K/DK);
    right = gsl_vector_complex_calloc((K/DK)*(K/DK));
    R_rev = gsl_matrix_complex_calloc(K/DK,K/DK);

    AB(P,x_kro_conju,tmp);
    InverseA(R_a,R_rev);

    AB(tmp,R_rev,tmp2);
    Ax(tmp2,inno,right);

    gsl_vector_complex_add(g,right);


}

void culc_P()
{
    gsl_matrix_complex* tmp,*R_rev;
    gsl_matrix_complex *tmp2;
    tmp = gsl_matrix_complex_calloc((K/DK*K/DK),K/DK);
    tmp2 = gsl_matrix_complex_calloc((K/DK)*(K/DK),(K/DK)*(K/DK));
    R_rev = gsl_matrix_complex_calloc(K/DK,K/DK);

    AB(P,x_kro_conju,tmp);
    InverseA(R_a,R_rev);
    AB(tmp,R_rev,tmp);


    ABT(tmp,x_kro,tmp2);
    AB(tmp2,P,tmp2);

    gsl_matrix_complex_add(P,tmp2);
}