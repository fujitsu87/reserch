void input_y_sub_vec(gsl_vector_complex* y_sub_vec)
{
    int i,k=0,t=0;
    for(i=0;i<(K/DK)*Tp;++i)
    {
        if(k>(K/DK-1))
        {
            k = 0;
            t++;
        }
        gsl_complex z = gsl_matrix_complex_get(y_sub,k,t);
        gsl_vector_complex_set(y_sub_vec,i,z);
        k++;
    }
    return;
}

void input_S_h(gsl_matrix_complex* S_h)
{
    int i,j,k=0,t;
    gsl_matrix_complex* x_expand;
    x_expand = gsl_matrix_complex_calloc((K/DK)*(K/DK),K/DK);
    for(t=0;t<Tp;++t)
    {
        for(i=0;i<(K/DK)*(K/DK);i=i+K/DK)
        {
            gsl_complex z = gsl_matrix_complex_get(x_h,k,t);
            for(j=0;j<(K/DK);++j)
            {
                gsl_matrix_complex_set(x_expand,i+j,j,z);   
            }
            ++k;
        }
        k=0;
        for(i=0;i<(K/DK);++i)
        {
            for(j=0;j<(K/DK)*(K/DK);++j)
            {
                gsl_complex z = gsl_matrix_complex_get(x_expand,j,i);
                gsl_matrix_complex_set(S_h,t*(K/DK)+i,j,z);
            }  
        }
    }

    x_expand =  GSLMatrixFree(x_expand);
    return;
}


gsl_matrix_complex *KPro(gsl_matrix_complex *a, gsl_matrix_complex *b) {
    int i, j, k, l;
    int m, p, n, q;
    m = a->size1;
    p = a->size2;
    n = b->size1;
    q = b->size2;

    gsl_matrix_complex *c = gsl_matrix_complex_alloc(m*n, p*q);
    gsl_complex da, db;

     for (i = 0; i < m; i++)    {
          for (j = 0; j < p; j++)   {
              da =gsl_matrix_complex_get (a, i, j);
              for (k = 0; k < n; k++)   {
                  for (l = 0; l < q; l++)   {
                      db = gsl_matrix_complex_get (b, k, l);
                      gsl_matrix_complex_set (c, n*i+k, q*j+l, gsl_complex_mul(da , db));                
                  }
              }
          }
      }

    return c;
}

gsl_matrix_complex*  input_R(gsl_matrix_complex* R)
{
    gsl_matrix_complex* I_k;
    I_k = gsl_matrix_complex_calloc(K/DK,K/DK);
    gsl_matrix_complex_set_identity( I_k );
    R = KPro(D,I_k);
    I_k =  GSLMatrixFree(I_k);

    return R;
}

gsl_matrix_complex*  input_R_x(gsl_matrix_complex* R)
{
    int i;
    gsl_complex z;
    for(i=0;i<K/DK;++i){
        z.dat[0] = Pk;
        z.dat[1] = 0.0;
        gsl_matrix_complex_set(R,i,i,z); 
    }
    return R;
}


void input_I_KTp(gsl_matrix_complex* I_KTp)
{
    gsl_matrix_complex_set_identity( I_KTp );
    return;
}

gsl_complex trace(gsl_matrix_complex* a)
{
    gsl_complex ans;
    int i;
    ans.dat[0]=0.0;ans.dat[1]=0.0;
    for(i=0;i<K/DK;i++)
    {
        ans.dat[0] += gsl_matrix_complex_get(a,i,i).dat[0];
        ans.dat[1] += gsl_matrix_complex_get(a,i,i).dat[1]; 
    }
    return ans;
}
gsl_matrix_complex* culc_lmmse_1st_term(
    gsl_matrix_complex* S_h,
    gsl_matrix_complex* R_g)
{
    gsl_matrix_complex* ans;
    ans = gsl_matrix_complex_calloc((K/DK)*(K/DK),(K/DK)*Tp);

    ABH(R_g,S_h,ans);

    return ans;
}
gsl_matrix_complex* culc_lmmse_2nd_term(
    gsl_matrix_complex* S_h,
    gsl_matrix_complex* R_g,
    gsl_matrix_complex* R_x,
    gsl_matrix_complex* I_KTp)
{
    int i,j;
    gsl_matrix_complex* tmp,*tmp2,*ans;
    gsl_complex z;
    tmp = gsl_matrix_complex_calloc((K/DK)*Tp,(K/DK)*(K/DK));
    tmp2 = gsl_matrix_complex_calloc((K/DK)*Tp,(K/DK)*Tp);
    ans = gsl_matrix_complex_calloc((K/DK)*Tp,(K/DK)*Tp);

    AB(S_h,R_g,tmp);
    ABH(tmp,S_h,ans);

    tmp =  GSLMatrixFree(tmp);
    for(i=0;i<DK-1;i++)
    {
        tmp = gsl_matrix_complex_calloc(K/DK,K/DK);
        gsl_matrix_complex_memcpy(tmp2,I_KTp); 

        AB(D_other,R_x,tmp);
        z = trace(tmp);
        gsl_matrix_complex_scale(tmp2,z);

        gsl_matrix_complex_add(ans,tmp2);
    }

    gsl_matrix_complex_memcpy(tmp2,I_KTp);
    z.dat[0] = N0;
    z.dat[1] = 0;
    gsl_matrix_complex_scale(tmp2,z);

    gsl_matrix_complex_add(ans,tmp2);

    InverseA(ans,ans);

    return ans;
}

void init_P(
    gsl_matrix_complex* S_h,
    gsl_matrix_complex* R_g,
    gsl_matrix_complex* R_x,
    gsl_matrix_complex* I_KTp,
    gsl_matrix_complex* term2)
{
    // gsl_matrix_complex *ans;
    gsl_matrix_complex* tmp,*tmp2;
    tmp = gsl_matrix_complex_calloc((K/DK)*(K/DK),(K/DK)*Tp);
    tmp2 = gsl_matrix_complex_calloc((K/DK)*(K/DK),(K/DK)*(K/DK));

    gsl_matrix_complex_memcpy(P,R_g);

    ABH(R_g,S_h,tmp);
    AB(tmp,term2,tmp);
    AB(tmp,S_h,tmp2);
    AB(tmp2,R_g,tmp2);

    gsl_matrix_complex_sub(P,tmp2); 

    return;
}
void lmmse()
{
    gsl_vector_complex* y_sub_vec;
    gsl_matrix_complex* S_h;
    gsl_matrix_complex* I_KTp;
    gsl_matrix_complex* R_g;
    gsl_matrix_complex* term1;
    gsl_matrix_complex* term2;
    gsl_matrix_complex* tmp;
    
    y_sub_vec = gsl_vector_complex_calloc((K/DK)*Tp);
    S_h = gsl_matrix_complex_calloc(K/DK*Tp,(K/DK)*(K/DK));
    I_KTp = gsl_matrix_complex_calloc(K/DK*Tp,(K/DK)*Tp);
    R_g = gsl_matrix_complex_calloc((K/DK)*(K/DK),(K/DK)*(K/DK));
    term1 = gsl_matrix_complex_calloc((K/DK)*(K/DK),(K/DK)*Tp);
    term2 = gsl_matrix_complex_calloc(K/DK*Tp,(K/DK)*Tp);
    tmp = gsl_matrix_complex_calloc((K/DK)*(K/DK),(K/DK)*Tp);

    input_y_sub_vec(y_sub_vec);
    input_S_h(S_h);
    input_I_KTp(I_KTp);
    R_g = input_R(R_g);
    R_x = input_R_x(R_x);
    
    term1 = culc_lmmse_1st_term(S_h,R_g);
    term2 = culc_lmmse_2nd_term(S_h,R_g,R_x,I_KTp);
    AB(term1,term2,tmp);
    Ax(tmp,y_sub_vec,g);
    init_P(S_h,R_g,R_x,I_KTp,term2);

    // gsl_vector_memcpy(g,);

// PrintVector(stdout,(K/DK)*Tp,y_sub_vec);
// PrintMatrix(stdout,(K/DK)*Tp,(K/DK)*Tp,I_KTp);

    
    y_sub_vec =  GSLVectorFree(y_sub_vec);
    I_KTp =  GSLMatrixFree(I_KTp);
    S_h =  GSLMatrixFree(S_h);
    R_g =  GSLMatrixFree(R_g);
}