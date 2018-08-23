#include <stdio.h>
#include <math.h>
#include <stdlib.h>

double culc_channel_se(double N0,double alpha,double P,double xi,double xi_tr,double eta,double gamma,double gammap)
{
    int L = 2;
    double zeta_tr = N0 + alpha*(double)L*P*eta + alpha*((double)L-1)*xi_tr*(1-eta);
    double zeta    = N0 + alpha*(double)L*P*eta + alpha*(double)L*xi*(1-eta);
    double tmp     = (P*gammap +(P-xi_tr)*((double)L-1)*gammap)/zeta_tr
                    +((P-xi)*(gamma-(double)L*gammap))/zeta;
    tmp = 1.0/tmp;
    eta = tmp/(1+tmp);

    printf("%g\n",eta);
    return eta;
}
double culc_xi_tr_se(double N0,double alpha,double P,double xi,double xi_tr,double eta,double gamma,double gammap)
{
    int L = 2;
    double zeta_tr = N0 + alpha*(double)L*P*eta + alpha*((double)L-1)*(1-eta)*xi_tr;
    double xi_b     = zeta_tr / (1-eta);
    xi_tr = tanh(xi_b);
    return xi_tr;
}
double culc_xi_se(double N0,double alpha,double P,double xi,double xi_tr,double eta,double gamma,double gammap)
{
    int L = 2;
    double zeta = N0 + alpha*(double)L*P*eta + alpha*(double)L*(1-eta)*xi;
    double xi_b     = zeta / (1-eta);
    xi = tanh(xi_b);
    printf("%g\n",xi);
    return xi;
}

int main(int argc, char *argv[])
{
    if(argc != 6)
    {
        printf("Input error. Please input as below \n" );
        printf("./a.out N K T Tp SNRdB\n");
        return 0;
    }
    int N,K,T,Tp,i,j;
    int x_loop[100],h_loop[100],out_loop;
    double alpha,gamma,gammap,N0,P,eta,xi,xi_tr;
    double pre,pre_out[2],e = 0.000001;
    N = atof(argv[1]);
    K = atof(argv[2]);
    T = atof(argv[3]);
    Tp = atof(argv[4]);
    N0 = pow(10.0,-1.0*atof(argv[5])/10.0);
    alpha = (double)K/(double)N;
    gamma = (double)T/(double)N;
    gammap = (double)Tp/(double)N;
    P = 1.0;
    eta = 1.0;
    xi = P;
    xi_tr = P;
    for(j=0;j<100;++j)
    {
        pre_out[0] = eta;
        pre_out[1] = xi;
        for(i=0;i<100;++i)
        {
            pre = eta;
            eta = culc_channel_se(N0,alpha,P,xi,xi_tr,eta,gamma,gammap);
            if(fabs(eta - pre) < e){
                printf("%d Hloop %d\n",j+1,i+1);
                h_loop[j]=i+1;
                break;
            }
        }
        for(i=0;i<100;++i){
            pre = xi;
            xi_tr = culc_xi_tr_se(N0,alpha,P,xi,xi_tr,eta,gamma,gammap);
            xi = culc_xi_se(N0,alpha,P,xi,xi_tr,eta,gamma,gammap);
            if(fabs(xi - pre) < e){
                printf("%d Xloop %d\n",j+1,i+1);
                x_loop[j]=i+1;
                break;
            }
        }
        if(fabs(eta - pre_out[0]) < e && fabs(xi - pre_out[1])){
            printf("out loop %d\n",j+1);
            out_loop = j+1;
            break;
        }
    }
    for(i=0;i<out_loop;++i){
        printf("H loop %d = %d\n",i,h_loop[i]);
    }
        for(i=0;i<out_loop;++i){
        printf("X loop %d = %d\n",i,x_loop[i]);
    }
    printf("out loop %d\n",out_loop);
    return 0;
}