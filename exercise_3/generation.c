# include <stdio.h>
# include <stdlib.h>
# include <math.h>


# define IA 16807
# define IM 2147483647
# define IM1 2147483563
# define IM2 2147483399
# define AM (1.0 / IM)
# define IMM1 (IM1-1)
# define IA1 40014
# define IA2 40692
# define IQ1 53668
# define IQ2 52774
# define IR1 12211
# define IR2 3791
# define IQ 127773
# define IR 2836
# define MASK 123459876
# define NTAB 32
# define NDIV (1+(IM-1)/NTAB)
# define EPS 1.2e-7
# define RNMX (1.0-EPS)
# define MBIG 1000000000
# define MSEED 161803398
# define MZ 0
# define FAC (1.0/MBIG)


float ran_0(long *idum){
    long k;
    float ans;

    *idum ^= MASK;
    k = (*idum) / IQ;
    *idum = IA * (*idum - k * IQ) - IR * k;

    if(*idum < 0) *idum += IM;
    *idum ^= MASK;

    return ans;
}

float ran_1(long *idum){
    int j;
    long k;
    static long iy = 0;
    static long iv[NTAB];
    float temp;

    if (*idum <= 0 || !iy){
        if (-(*idum) < 1) *idum = 1;
        else *idum = -(*idum);
        
        for (j = NTAB + 7; j >= 0; j--){
            k = (*idum) / IQ;
            *idum = IA * (*idum - k * IQ) - IR * k;
            if (*idum < 0) *idum += IM;
            if (j < NTAB) iv[j] = *idum;
        }
        iy = iv[0];
    }
    k = (*idum) / IQ;
    *idum = IA * (*idum - k * IQ) - IR * k;
    if (*idum < 0) *idum += IM;
    j = iy / NDIV;
    iy = iv[j];
    iv[j] = *idum;
    if ((temp = AM * iy) > RNMX) return RNMX;
    else return temp;
}

float ran_2(long *idum){
    int j;
    long k;
    static long idum2 = 123456789;
    static long iy = 0;
    static long iv[NTAB];
    float temp;

    if (*idum <= 0) {
        if (-(*idum) < 1) *idum = 1;
        else *idum = -(*idum);
        idum2 = (*idum);
        
        for (j = NTAB + 7;j >= 0; j--){
        k = (*idum) / IQ1;
        *idum = IA1 * (*idum - k * IQ1) - k * IR1;
        if (*idum < 0) *idum += IM1;
        if (j < NTAB) iv[j] = *idum;
        }
        iy=iv[0];
    }
    k = (*idum) / IQ1;
    *idum = IA1 * (*idum - k * IQ1) - k * IR1;
    if (*idum < 0) *idum += IM1;
    k = idum2/IQ2;
    idum2 = IA2 * (idum2 - k * IQ2) - k * IR2;
    if (idum2 < 0) idum2 += IM2;
    j = iy / NDIV;
    iv[j] = *idum;
    if (iy < 1) iy += IMM1;
    if ((temp = AM * iy) > RNMX) return RNMX;
    else return temp;
}

/*float ran_3(long *idum){
    static int inext,inextp;
    static long ma[56];
    static int iff=0;
    long mj,mk;
    int i,ii,k;

    if(*idum < 0 || iff == 0){
        iff = 1;
        mj = labs(MSEED - labs(*idum));
        mj %= MBIG;
        ma[55] = mj;
        mk = 1;
        
        for(i = 1; i <= 54; i++){
            ii = (21 * i) % 55;
            ma[ii] = mk;
            mk = mj - mk;
            if(mk < MZ) mk += MBIG;
            mj=ma[ii];
        }
        
        for(k = 1; k <= 4; k++){
            for(i = 1; i <= 55; i++){
                ma[i] -= ma[1 + (i + 30) % 55];
                if (ma[i] < MZ) ma[i] += MBIG;
            }
            inext = 0;
            inextp = 31;
            *idum = 1;
        }
        if(++inext == 56) inext = 1;
        if(++inextp == 56) inextp = 1;
        if(mj < MZ) mj += MBIG;
    return mj * FAC;
}*/

float gasdev(long *idum){
    float ran_1(long *idum);
    static int iset = 0;
    static float gset;
    float fac, rsq, v1, v2;

    if (*idum < 0) iset = 0;
    if (iset == 0){
        do {
            v1 = 2.0*ran_1(idum) - 1.0;
            v2 = 2.0*ran_1(idum) - 1.0;
            rsq = v1 * v1 + v2 * v2;
        } while(rsq >= 1.0 || rsq == 0.0);
            fac = sqrt(-2.0 * log(rsq) / rsq);
            gset = v1 * fac;
            iset = 1;
            return v2 * fac;
    } else{
        iset = 0;
        return gset;
    }
}
