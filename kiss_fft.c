/*
 *  Copyright (c) 2003-2010, Mark Borgerding. All rights reserved.
 *  This file is part of KISS FFT - https://github.com/mborgerding/kissfft
 *
 *  SPDX-License-Identifier: BSD-3-Clause
 *  See COPYING file for more information.
 */


#include "_kiss_fft_guts.h"
/* The guts header contains all the multiplication and addition macros that are defined for
 fixed or floating point complex numbers.  It also delares the kf_ internal functions.
 */

static void kf_bfly2(
        kiss_fft_cpx * Fout,
        const size_t fstride,
        const kiss_fft_cfg st,
        int m
        )
{
    kiss_fft_cpx * Fout2;
    kiss_fft_cpx * tw1 = st->twiddles;
    kiss_fft_cpx t;
    Fout2 = Fout + m;
    do{
        C_FIXDIV(*Fout,2); C_FIXDIV(*Fout2,2);

        C_MUL (t,  *Fout2 , *tw1);
        tw1 += fstride;
        C_SUB( *Fout2 ,  *Fout , t );
        C_ADDTO( *Fout ,  t );
        ++Fout2;
        ++Fout;
    }while (--m);
}

static void kf_bfly4(
        kiss_fft_cpx * Fout,
        const size_t fstride,
        const kiss_fft_cfg st,
        const size_t m
        )
{
    kiss_fft_cpx *tw1,*tw2,*tw3;
    kiss_fft_cpx scratch[6];
    size_t k=m;
    const size_t m2=2*m;
    const size_t m3=3*m;


    tw3 = tw2 = tw1 = st->twiddles;

    do {
        C_FIXDIV(*Fout,4); C_FIXDIV(Fout[m],4); C_FIXDIV(Fout[m2],4); C_FIXDIV(Fout[m3],4);

        C_MUL(scratch[0],Fout[m] , *tw1 );
        C_MUL(scratch[1],Fout[m2] , *tw2 );
        C_MUL(scratch[2],Fout[m3] , *tw3 );

        C_SUB( scratch[5] , *Fout, scratch[1] );
        C_ADDTO(*Fout, scratch[1]);
        C_ADD( scratch[3] , scratch[0] , scratch[2] );
        C_SUB( scratch[4] , scratch[0] , scratch[2] );
        C_SUB( Fout[m2], *Fout, scratch[3] );
        tw1 += fstride;
        tw2 += fstride*2;
        tw3 += fstride*3;
        C_ADDTO( *Fout , scratch[3] );

        if(st->inverse) {
            Fout[m].r = scratch[5].r - scratch[4].i;
            Fout[m].i = scratch[5].i + scratch[4].r;
            Fout[m3].r = scratch[5].r + scratch[4].i;
            Fout[m3].i = scratch[5].i - scratch[4].r;
        }else{
            Fout[m].r = scratch[5].r + scratch[4].i;
            Fout[m].i = scratch[5].i - scratch[4].r;
            Fout[m3].r = scratch[5].r - scratch[4].i;
            Fout[m3].i = scratch[5].i + scratch[4].r;
        }
        ++Fout;
    }while(--k);
}

/*
static void kf_bfly4_m4(
    kiss_fft_cpx* Fout,
    const size_t fstride,
    const kiss_fft_cfg st
)
{
    kiss_fft_cpx* tw1, * tw2, * tw3;
    kiss_fft_cpx scratch[6];
    
    const size_t m = 4;
    size_t k = m;
    const size_t m2 = 8; //  2 * 4;
    const size_t m3 = 12; // 3 * 4;


    tw3 = tw2 = tw1 = st->twiddles;

    

    // do {
        // C_FIXDIV(*Fout, 4); C_FIXDIV(Fout[m], 4); C_FIXDIV(Fout[m2], 4); C_FIXDIV(Fout[m3], 4);

        // 1
        C_MUL(scratch[0], Fout[4], *tw1);
        C_MUL(scratch[1], Fout[8], *tw2);
        C_MUL(scratch[2], Fout[12], *tw3);

        C_SUB(scratch[5], *Fout, scratch[1]);
        C_ADDTO(Fout[0], scratch[1]);
        C_ADD(scratch[3], scratch[0], scratch[2]);
        C_SUB(scratch[4], scratch[0], scratch[2]);
        C_SUB(Fout[8], *Fout, scratch[3]);
        tw1 += fstride;
        tw2 += fstride * 2;
        tw3 += fstride * 3;
        C_ADDTO(*Fout, scratch[3]);

        if (st->inverse) {
            Fout[4].r = scratch[5].r - scratch[4].i;
            Fout[4].i = scratch[5].i + scratch[4].r;
            Fout[12].r = scratch[5].r + scratch[4].i;
            Fout[12].i = scratch[5].i - scratch[4].r;
        }
        else {
            Fout[4].r = scratch[5].r + scratch[4].i;
            Fout[4].i = scratch[5].i - scratch[4].r;
            Fout[12].r = scratch[5].r - scratch[4].i;
            Fout[12].i = scratch[5].i + scratch[4].r;
        }
        // ++Fout;

        // 2
        C_MUL(scratch[0], Fout[4+1], *tw1);
        C_MUL(scratch[1], Fout[8 + 1], *tw2);
        C_MUL(scratch[2], Fout[12 + 1], *tw3);

        C_SUB(scratch[5], Fout[1], scratch[1]);
        C_ADDTO(Fout[1], scratch[1]);
        C_ADD(scratch[3], scratch[0], scratch[2]);
        C_SUB(scratch[4], scratch[0], scratch[2]);
        C_SUB(Fout[8+1], Fout[1], scratch[3]);
        tw1 += fstride;
        tw2 += fstride * 2;
        tw3 += fstride * 3;
        C_ADDTO(Fout[1], scratch[3]);

        if (st->inverse) {
            Fout[5].r = scratch[5].r - scratch[4].i;
            Fout[5].i = scratch[5].i + scratch[4].r;
            Fout[13].r = scratch[5].r + scratch[4].i;
            Fout[13].i = scratch[5].i - scratch[4].r;
        }
        else {
            Fout[5].r = scratch[5].r + scratch[4].i;
            Fout[5].i = scratch[5].i - scratch[4].r;
            Fout[13].r = scratch[5].r - scratch[4].i;
            Fout[13].i = scratch[5].i + scratch[4].r;
        }
        // ++Fout;

        // 3
        C_MUL(scratch[0], Fout[4+2], *tw1);
        C_MUL(scratch[1], Fout[8 + 2], *tw2);
        C_MUL(scratch[2], Fout[12 + 2], *tw3);

        C_SUB(scratch[5], Fout[2], scratch[1]);
        C_ADDTO(Fout[2], scratch[1]);
        C_ADD(scratch[3], scratch[0], scratch[2]);
        C_SUB(scratch[4], scratch[0], scratch[2]);
        C_SUB(Fout[8+2], Fout[2], scratch[3]);
        tw1 += fstride;
        tw2 += fstride * 2;
        tw3 += fstride * 3;
        C_ADDTO(Fout[2], scratch[3]);

        if (st->inverse) {
            Fout[4+2].r = scratch[5].r - scratch[4].i;
            Fout[4 + 2].i = scratch[5].i + scratch[4].r;
            Fout[12 + 2].r = scratch[5].r + scratch[4].i;
            Fout[12 + 2].i = scratch[5].i - scratch[4].r;
        }
        else {
            Fout[4 + 2].r = scratch[5].r + scratch[4].i;
            Fout[4 + 2].i = scratch[5].i - scratch[4].r;
            Fout[12 + 2].r = scratch[5].r - scratch[4].i;
            Fout[12 + 2].i = scratch[5].i + scratch[4].r;
        }
        // ++Fout;

        // 3
        C_MUL(scratch[0], Fout[4+3], *tw1);
        C_MUL(scratch[1], Fout[8 + 3], *tw2);
        C_MUL(scratch[2], Fout[12 + 3], *tw3);

        C_SUB(scratch[5], Fout[3], scratch[1]);
        C_ADDTO(Fout[3], scratch[1]);
        C_ADD(scratch[3], scratch[0], scratch[2]);
        C_SUB(scratch[4], scratch[0], scratch[2]);
        C_SUB(Fout[8 + 3], Fout[3], scratch[3]);
        tw1 += fstride;
        tw2 += fstride * 2;
        tw3 += fstride * 3;
        C_ADDTO(Fout[3], scratch[3]);

        if (st->inverse) {
            Fout[4+3].r = scratch[5].r - scratch[4].i;
            Fout[4 + 3].i = scratch[5].i + scratch[4].r;
            Fout[12 + 3].r = scratch[5].r + scratch[4].i;
            Fout[12 + 3].i = scratch[5].i - scratch[4].r;
        }
        else {
            Fout[4 + 3].r = scratch[5].r + scratch[4].i;
            Fout[4 + 3].i = scratch[5].i - scratch[4].r;
            Fout[12 + 3].r = scratch[5].r - scratch[4].i;
            Fout[12 + 3].i = scratch[5].i + scratch[4].r;
        }
        // ++Fout;



    // } while (--k);
}


static void kf_bfly4_m4_remapped(
    kiss_fft_cpx* fast_Fout,
    const size_t fstride,
    const kiss_fft_cfg st
)
{
    kiss_fft_cpx* tw1, * tw2, * tw3;
    kiss_fft_cpx scratch[6];

    const size_t m = 4;
    size_t k = m;
    const size_t m2 = 8; //  2 * 4;
    const size_t m3 = 12; // 3 * 4;


    tw3 = tw2 = tw1 = st->twiddles;



    // do {
        // C_FIXDIV(*Fout, 4); C_FIXDIV(Fout[m], 4); C_FIXDIV(Fout[m2], 4); C_FIXDIV(Fout[m3], 4);

        // 1
    C_MUL(scratch[0], fast_Fout[0], *tw1);
    C_MUL(scratch[1], fast_Fout[1], *tw2);
    C_MUL(scratch[2], fast_Fout[2], *tw3);

    C_SUB(scratch[5], fast_Fout[3], scratch[1]);
    C_ADDTO(fast_Fout[3], scratch[1]);
    C_ADD(scratch[3], scratch[0], scratch[2]);
    C_SUB(scratch[4], scratch[0], scratch[2]);
    C_SUB(fast_Fout[1], fast_Fout[3], scratch[3]);
    tw1 += fstride;
    tw2 += fstride * 2;
    tw3 += fstride * 3;
    C_ADDTO(fast_Fout[3], scratch[3]);

    if (st->inverse) {
        fast_Fout[0].r = scratch[5].r - scratch[4].i;
        fast_Fout[0].i = scratch[5].i + scratch[4].r;
        fast_Fout[2].r = scratch[5].r + scratch[4].i;
        fast_Fout[2].i = scratch[5].i - scratch[4].r;
    }
    else {
        fast_Fout[0].r = scratch[5].r + scratch[4].i;
        fast_Fout[0].i = scratch[5].i - scratch[4].r;
        fast_Fout[2].r = scratch[5].r - scratch[4].i;
        fast_Fout[2].i = scratch[5].i + scratch[4].r;
    }
    // ++Fout;

    // 2
    C_MUL(scratch[0], fast_Fout[4], *tw1);
    C_MUL(scratch[1], fast_Fout[5], *tw2);
    C_MUL(scratch[2], fast_Fout[6], *tw3);

    C_SUB(scratch[5], fast_Fout[7], scratch[1]);
    C_ADDTO(fast_Fout[7], scratch[1]);
    C_ADD(scratch[3], scratch[0], scratch[2]);
    C_SUB(scratch[4], scratch[0], scratch[2]);
    C_SUB(fast_Fout[5], fast_Fout[7], scratch[3]);
    tw1 += fstride;
    tw2 += fstride * 2;
    tw3 += fstride * 3;
    C_ADDTO(fast_Fout[7], scratch[3]);

    if (st->inverse) {
        fast_Fout[4].r = scratch[5].r - scratch[4].i;
        fast_Fout[4].i = scratch[5].i + scratch[4].r;
        fast_Fout[6].r = scratch[5].r + scratch[4].i;
        fast_Fout[6].i = scratch[5].i - scratch[4].r;
    }
    else {
        fast_Fout[4].r = scratch[5].r + scratch[4].i;
        fast_Fout[4].i = scratch[5].i - scratch[4].r;
        fast_Fout[6].r = scratch[5].r - scratch[4].i;
        fast_Fout[6].i = scratch[5].i + scratch[4].r;
    }
    // ++Fout;

    // 3
    C_MUL(scratch[0], fast_Fout[8], *tw1);
    C_MUL(scratch[1], fast_Fout[9], *tw2);
    C_MUL(scratch[2], fast_Fout[10], *tw3);

    C_SUB(scratch[5], fast_Fout[11], scratch[1]);
    C_ADDTO(fast_Fout[11], scratch[1]);
    C_ADD(scratch[3], scratch[0], scratch[2]);
    C_SUB(scratch[4], scratch[0], scratch[2]);
    C_SUB(fast_Fout[9], fast_Fout[11], scratch[3]);
    tw1 += fstride;
    tw2 += fstride * 2;
    tw3 += fstride * 3;
    C_ADDTO(fast_Fout[11], scratch[3]);

    if (st->inverse) {
        fast_Fout[8].r = scratch[5].r - scratch[4].i;
        fast_Fout[8].i = scratch[5].i + scratch[4].r;
        fast_Fout[10].r = scratch[5].r + scratch[4].i;
        fast_Fout[10].i = scratch[5].i - scratch[4].r;
    }
    else {
        fast_Fout[8].r = scratch[5].r + scratch[4].i;
        fast_Fout[8].i = scratch[5].i - scratch[4].r;
        fast_Fout[10].r = scratch[5].r - scratch[4].i;
        fast_Fout[10].i = scratch[5].i + scratch[4].r;
    }
    // ++Fout;

    // 3
    C_MUL(scratch[0], fast_Fout[12], *tw1);
    C_MUL(scratch[1], fast_Fout[13], *tw2);
    C_MUL(scratch[2], fast_Fout[14], *tw3);

    C_SUB(scratch[5], fast_Fout[15], scratch[1]);
    C_ADDTO(fast_Fout[15], scratch[1]);
    C_ADD(scratch[3], scratch[0], scratch[2]);
    C_SUB(scratch[4], scratch[0], scratch[2]);
    C_SUB(fast_Fout[13], fast_Fout[15], scratch[3]);
    tw1 += fstride;
    tw2 += fstride * 2;
    tw3 += fstride * 3;
    C_ADDTO(fast_Fout[15], scratch[3]);

    if (st->inverse) {
        fast_Fout[12].r = scratch[5].r - scratch[4].i;
        fast_Fout[12].i = scratch[5].i + scratch[4].r;
        fast_Fout[14].r = scratch[5].r + scratch[4].i;
        fast_Fout[14].i = scratch[5].i - scratch[4].r;
    }
    else {
        fast_Fout[12].r = scratch[5].r + scratch[4].i;
        fast_Fout[12].i = scratch[5].i - scratch[4].r;
        fast_Fout[14].r = scratch[5].r - scratch[4].i;
        fast_Fout[14].i = scratch[5].i + scratch[4].r;
    }
    // ++Fout;



// } while (--k);
}

static void f_map_Fout_to_fast_Fout(const kiss_fft_cpx* Fin, kiss_fft_cpx* fast_Fout)
{
    fast_Fout[0].r = Fin[4].r;
    fast_Fout[0].i = Fin[4].i;
    fast_Fout[1].r = Fin[8].r;
    fast_Fout[1].i = Fin[8].i;
    fast_Fout[2].r = Fin[12].r;
    fast_Fout[2].i = Fin[12].i;
    fast_Fout[3].r = Fin[0].r;
    fast_Fout[3].i = Fin[0].i;

    fast_Fout[4].r = Fin[5].r;
    fast_Fout[4].i = Fin[5].i;
    fast_Fout[5].r = Fin[9].r;
    fast_Fout[5].i = Fin[9].i;
    fast_Fout[6].r = Fin[13].r;
    fast_Fout[6].i = Fin[13].i;
    fast_Fout[7].r = Fin[1].r;
    fast_Fout[7].i = Fin[1].i;

    fast_Fout[8].r = Fin[6].r;
    fast_Fout[8].i = Fin[6].i;
    fast_Fout[9].r = Fin[10].r;
    fast_Fout[9].i = Fin[10].i;
    fast_Fout[10].r = Fin[14].r;
    fast_Fout[10].i = Fin[14].i;
    fast_Fout[11].r = Fin[2].r;
    fast_Fout[11].i = Fin[2].i;

    fast_Fout[12].r = Fin[7].r;
    fast_Fout[12].i = Fin[7].i;
    fast_Fout[13].r = Fin[11].r;
    fast_Fout[13].i = Fin[11].i;
    fast_Fout[14].r = Fin[15].r;
    fast_Fout[14].i = Fin[15].i;
    fast_Fout[15].r = Fin[3].r;
    fast_Fout[15].i = Fin[3].i;

}

static void f_map_fast_Fout_to_Fout(const kiss_fft_cpx* fast_Fout, kiss_fft_cpx* Fout)
{
    Fout[4].r = fast_Fout[0].r; // = Fin[4].r;
    Fout[4].i = fast_Fout[0].i; // = Fin[4].i;
    Fout[8].r = fast_Fout[1].r; // = Fin[8].r;
    Fout[8].i = fast_Fout[1].i; // = Fin[8].i;
    Fout[12].r = fast_Fout[2].r; // = Fin[12].r;
    Fout[12].i = fast_Fout[2].i; // = Fin[12].i;
    Fout[0].r = fast_Fout[3].r; // = Fin[0].r;
    Fout[0].i = fast_Fout[3].i; // = Fin[0].i;

    Fout[5].r = fast_Fout[4].r; //  = Fin[5].r;
    Fout[5].i = fast_Fout[4].i; //  = Fin[5].i;
    Fout[9].r = fast_Fout[5].r; //  = Fin[9].r;
    Fout[9].i = fast_Fout[5].i; //  = Fin[9].i;
    Fout[13].r = fast_Fout[6].r; //  = Fin[13].r;
    Fout[13].i = fast_Fout[6].i; //  = Fin[13].i;
    Fout[1].r = fast_Fout[7].r; //  = Fin[1].r;
    Fout[1].i = fast_Fout[7].i; //  = Fin[1].i;

    Fout[6].r = fast_Fout[8].r; //   = Fin[6].r;
    Fout[6].i = fast_Fout[8].i; //   = Fin[6].i;
    Fout[10].r = fast_Fout[9].r; //   = Fin[10].r;
    Fout[10].i = fast_Fout[9].i; //   = Fin[10].i;
    Fout[14].r = fast_Fout[10].r; //   = Fin[14].r;
    Fout[14].i = fast_Fout[10].i; //   = Fin[14].i;
    Fout[2].r = fast_Fout[11].r; //  = Fin[2].r;
    Fout[2].i = fast_Fout[11].i; //  = Fin[2].i;

    Fout[7].r = fast_Fout[12].r; // = Fin[7].r;
    Fout[7].i = fast_Fout[12].i; // = Fin[7].i;
    Fout[11].r = fast_Fout[13].r; // = Fin[11].r;
    Fout[11].i = fast_Fout[13].i; // = Fin[11].i;
    Fout[15].r = fast_Fout[14].r; //= Fin[15].r;
    Fout[15].i = fast_Fout[14].i; // = Fin[15].i;
    Fout[3].r = fast_Fout[15].r; //= Fin[3].r;
    Fout[3].i = fast_Fout[15].i; // = Fin[3].i;

}

*/

static void kf_bfly4_m4_multiple_scratches(
    kiss_fft_cpx* Fout,
    const size_t fstride,
    const kiss_fft_cfg st
)
{
    kiss_fft_cpx* tw1, * tw2, * tw3;

    kiss_fft_cpx scratch_interleaved[6*4];
    kiss_fft_cpx* p_scratch_interleaved_a;
    kiss_fft_cpx* p_scratch_interleaved_b;
    kiss_fft_cpx* p_Fout;
    kiss_fft_cpx* p_Fout2;

    const size_t m = 4;
    size_t k = m;
    const size_t m2 = 8; //  2 * 4;
    const size_t m3 = 12; // 3 * 4;


    tw3 = tw2 = tw1 = st->twiddles;


#if 1
    C_MUL(scratch_interleaved[0],   Fout[4], *(tw1 + 0));              //  C_MUL(scratch1[0], Fout[4], *(tw1 + 0));
    C_MUL(scratch_interleaved[0+1], Fout[4+1], *(tw1 + fstride));    //  C_MUL(scratch2[0], Fout[4 + 1], *(tw1 + fstride));
    C_MUL(scratch_interleaved[0+2], Fout[4+2], *(tw1 + 2 * fstride));//   C_MUL(scratch3[0], Fout[4 + 2], *(tw1 + 2 * fstride));
    C_MUL(scratch_interleaved[0+3], Fout[4+3], *(tw1 + 3 * fstride)); //  C_MUL(scratch4[0], Fout[4 + 3], *(tw1 + 3 * fstride));

    C_MUL(scratch_interleaved[1*4],    Fout[8], *(tw2 + 0));
    C_MUL(scratch_interleaved[1*4 +1], Fout[8+1], *(tw2 + 2 * fstride));
    C_MUL(scratch_interleaved[1*4 +2], Fout[8+2], *(tw2 + 2 * 2 * fstride));
    C_MUL(scratch_interleaved[1*4 +3], Fout[8+3], *(tw2 + 3 * 2 * fstride));

    C_MUL(scratch_interleaved[2*4],    Fout[12], *(tw3 + 0));
    C_MUL(scratch_interleaved[2*4 +1], Fout[12+1], *(tw3 + 3 * fstride));
    C_MUL(scratch_interleaved[2*4 +2], Fout[12+2], *(tw3 + 2 * 3 * fstride));
    C_MUL(scratch_interleaved[2*4 +3], Fout[12+3], *(tw3 + 3 * 3 * fstride));

    C_SUB(scratch_interleaved[5*4],   *Fout,    scratch_interleaved[1*4]);
    C_SUB(scratch_interleaved[5*4 +1], Fout[1], scratch_interleaved[1*4 +1]);
    C_SUB(scratch_interleaved[5*4 +2], Fout[2], scratch_interleaved[1*4 +2]);
    C_SUB(scratch_interleaved[5*4 +3], Fout[3], scratch_interleaved[1*4 +3]);

    C_ADDTO(Fout[0], scratch_interleaved[1*4]);
    C_ADDTO(Fout[1], scratch_interleaved[1*4 +1]);
    C_ADDTO(Fout[2], scratch_interleaved[1*4 +2]);
    C_ADDTO(Fout[3], scratch_interleaved[1*4 +3]);

    C_ADD(scratch_interleaved[3*4],    scratch_interleaved[0],   scratch_interleaved[2*4]);
    C_ADD(scratch_interleaved[3*4 +1], scratch_interleaved[0+1], scratch_interleaved[2*4 +1]);
    C_ADD(scratch_interleaved[3*4 +2], scratch_interleaved[0+2], scratch_interleaved[2*4 +2]);
    C_ADD(scratch_interleaved[3*4 +3], scratch_interleaved[0+3], scratch_interleaved[2*4 +3]);

    C_SUB(scratch_interleaved[4*4],   scratch_interleaved[0],    scratch_interleaved[2*4]);
    C_SUB(scratch_interleaved[4*4 +1], scratch_interleaved[0+1], scratch_interleaved[2*4 +1]);
    C_SUB(scratch_interleaved[4*4 +2], scratch_interleaved[0+2], scratch_interleaved[2*4 +2]);
    C_SUB(scratch_interleaved[4*4 +3], scratch_interleaved[0+3], scratch_interleaved[2*4 +3]);

    p_scratch_interleaved_a = &scratch_interleaved[3 * 4];
    p_Fout = &Fout[8];
    p_Fout2 = &Fout[0];
    (*p_Fout).r = (*p_Fout2).r - (*p_scratch_interleaved_a).r;
    (*p_Fout).i = (*p_Fout2).i - (*p_scratch_interleaved_a).i; p_scratch_interleaved_a++; p_Fout++; p_Fout2++;
    (*p_Fout).r = (*p_Fout2).r - (*p_scratch_interleaved_a).r;
    (*p_Fout).i = (*p_Fout2).i - (*p_scratch_interleaved_a).i; p_scratch_interleaved_a++; p_Fout++; p_Fout2++;
    (*p_Fout).r = (*p_Fout2).r - (*p_scratch_interleaved_a).r;
    (*p_Fout).i = (*p_Fout2).i - (*p_scratch_interleaved_a).i; p_scratch_interleaved_a++; p_Fout++; p_Fout2++;
    (*p_Fout).r = (*p_Fout2).r - (*p_scratch_interleaved_a).r;
    (*p_Fout).i = (*p_Fout2).i - (*p_scratch_interleaved_a).i; 

    // C_SUB(Fout[8],    *Fout,    scratch_interleaved[3*4]);
    // C_SUB(Fout[8 + 1], Fout[1], scratch_interleaved[3*4 +1]);
    //  C_SUB(Fout[8 + 2], Fout[2], scratch_interleaved[3*4 +2]);
    // C_SUB(Fout[8 + 3], Fout[3], scratch_interleaved[3*4 +3]);

    p_scratch_interleaved_a = &scratch_interleaved[3 * 4];
    p_Fout = &Fout[0];
    (*p_Fout).r += (*p_scratch_interleaved_a).r;  
    (*p_Fout).i += (*p_scratch_interleaved_a).i; p_scratch_interleaved_a++; p_Fout++;
    (*p_Fout).r += (*p_scratch_interleaved_a).r;
    (*p_Fout).i += (*p_scratch_interleaved_a).i; p_scratch_interleaved_a++; p_Fout++;
    (*p_Fout).r += (*p_scratch_interleaved_a).r;
    (*p_Fout).i += (*p_scratch_interleaved_a).i; p_scratch_interleaved_a++; p_Fout++;
    (*p_Fout).r += (*p_scratch_interleaved_a).r;
    (*p_Fout).i += (*p_scratch_interleaved_a).i;


    // C_ADDTO(*Fout,   scratch_interleaved[3*4]);
    // C_ADDTO(Fout[1], scratch_interleaved[3*4 +1]);
    // C_ADDTO(Fout[2], scratch_interleaved[3*4 +2]);
    // C_ADDTO(Fout[3], scratch_interleaved[3*4 +3]);

    if (st->inverse) {

#if 1
        p_scratch_interleaved_a = &scratch_interleaved[5 * 4];
        p_scratch_interleaved_b = &scratch_interleaved[4 * 4];
        p_Fout                  = &Fout[4];

        (*p_Fout).r   = (*p_scratch_interleaved_a).r - (*p_scratch_interleaved_b).i; // ->
        (*p_Fout).i   = (*p_scratch_interleaved_a).i + (*p_scratch_interleaved_b).r;
        p_scratch_interleaved_a++; p_scratch_interleaved_b++; p_Fout++;

        (*p_Fout).r = (*p_scratch_interleaved_a).r - (*p_scratch_interleaved_b).i;
        (*p_Fout).i = (*p_scratch_interleaved_a).i + (*p_scratch_interleaved_b).r;
        p_scratch_interleaved_a++; p_scratch_interleaved_b++; p_Fout++;

        (*p_Fout).r = (*p_scratch_interleaved_a).r - (*p_scratch_interleaved_b).i;
        (*p_Fout).i = (*p_scratch_interleaved_a).i + (*p_scratch_interleaved_b).r;
        p_scratch_interleaved_a++; p_scratch_interleaved_b++; p_Fout++;

        (*p_Fout).r = (*p_scratch_interleaved_a).r - (*p_scratch_interleaved_b).i;
        (*p_Fout).i = (*p_scratch_interleaved_a).i + (*p_scratch_interleaved_b).r;



        p_scratch_interleaved_a = &scratch_interleaved[5 * 4];
        p_scratch_interleaved_b = &scratch_interleaved[4 * 4];
        p_Fout                  = &Fout[12];

        (*p_Fout).r = (*p_scratch_interleaved_a).r + (*p_scratch_interleaved_b).i;
        (*p_Fout).i = (*p_scratch_interleaved_a).i - (*p_scratch_interleaved_b).r;
        p_scratch_interleaved_a++;  p_scratch_interleaved_b++; p_Fout++;

        (*p_Fout).r = (*p_scratch_interleaved_a).r + (*p_scratch_interleaved_b).i;
        (*p_Fout).i = (*p_scratch_interleaved_a).i - (*p_scratch_interleaved_b).r;
        p_scratch_interleaved_a++;  p_scratch_interleaved_b++; p_Fout++;

        (*p_Fout).r = (*p_scratch_interleaved_a).r + (*p_scratch_interleaved_b).i;
        (*p_Fout).i = (*p_scratch_interleaved_a).i - (*p_scratch_interleaved_b).r;
        p_scratch_interleaved_a++;  p_scratch_interleaved_b++; p_Fout++;

        (*p_Fout).r = (*p_scratch_interleaved_a).r + (*p_scratch_interleaved_b).i;
        (*p_Fout).i = (*p_scratch_interleaved_a).i - (*p_scratch_interleaved_b).r;
#else
        Fout[4].r = scratch_interleaved[5 * 4].r - scratch_interleaved[4 * 4].i;
        Fout[4].i = scratch_interleaved[5 * 4].i + scratch_interleaved[4 * 4].r;
        Fout[4 + 1].r = scratch_interleaved[5 * 4 + 1].r - scratch_interleaved[4 * 4 + 1].i;
        Fout[4 + 1].i = scratch_interleaved[5 * 4 + 1].i + scratch_interleaved[4 * 4 + 1].r;
        Fout[4 + 2].r = scratch_interleaved[5 * 4 + 2].r - scratch_interleaved[4 * 4 + 2].i;
        Fout[4 + 2].i = scratch_interleaved[5 * 4 + 2].i + scratch_interleaved[4 * 4 + 2].r;
        Fout[4 + 3].r = scratch_interleaved[5 * 4 + 3].r - scratch_interleaved[4 * 4 + 3].i;
        Fout[4 + 3].i = scratch_interleaved[5 * 4 + 3].i + scratch_interleaved[4 * 4 + 3].r;

        Fout[12].r = scratch_interleaved[5 * 4].r + scratch_interleaved[4 * 4].i;
        Fout[12].i = scratch_interleaved[5 * 4].i - scratch_interleaved[4 * 4].r;
        Fout[12 + 1].r = scratch_interleaved[5 * 4 + 1].r + scratch_interleaved[4 * 4 + 1].i;
        Fout[12 + 1].i = scratch_interleaved[5 * 4 + 1].i - scratch_interleaved[4 * 4 + 1].r;
        Fout[12 + 2].r = scratch_interleaved[5 * 4 + 2].r + scratch_interleaved[4 * 4 + 2].i;
        Fout[12 + 2].i = scratch_interleaved[5 * 4 + 2].i - scratch_interleaved[4 * 4 + 2].r;
        Fout[12 + 3].r = scratch_interleaved[5 * 4 + 3].r + scratch_interleaved[4 * 4 + 3].i;
        Fout[12 + 3].i = scratch_interleaved[5 * 4 + 3].i - scratch_interleaved[4 * 4 + 3].r;

#endif

    }
    else {

        p_scratch_interleaved_a = &scratch_interleaved[5 * 4];
        p_scratch_interleaved_b = &scratch_interleaved[4 * 4];
        p_Fout                  = &Fout[4];

        (*p_Fout).r = (*p_scratch_interleaved_a).r + (*p_scratch_interleaved_b).i; // ->
        (*p_Fout).i = (*p_scratch_interleaved_a).i - (*p_scratch_interleaved_b).r;
        p_scratch_interleaved_a++; p_scratch_interleaved_b++; p_Fout++;

        (*p_Fout).r = (*p_scratch_interleaved_a).r + (*p_scratch_interleaved_b).i;
        (*p_Fout).i = (*p_scratch_interleaved_a).i - (*p_scratch_interleaved_b).r;
        p_scratch_interleaved_a++; p_scratch_interleaved_b++; p_Fout++;

        (*p_Fout).r = (*p_scratch_interleaved_a).r + (*p_scratch_interleaved_b).i;
        (*p_Fout).i = (*p_scratch_interleaved_a).i - (*p_scratch_interleaved_b).r;
        p_scratch_interleaved_a++; p_scratch_interleaved_b++; p_Fout++;

        (*p_Fout).r = (*p_scratch_interleaved_a).r + (*p_scratch_interleaved_b).i;
        (*p_Fout).i = (*p_scratch_interleaved_a).i - (*p_scratch_interleaved_b).r;



        p_scratch_interleaved_a = &scratch_interleaved[5 * 4];
        p_scratch_interleaved_b = &scratch_interleaved[4 * 4];
        p_Fout                  = &Fout[12];

        (*p_Fout).r = (*p_scratch_interleaved_a).r - (*p_scratch_interleaved_b).i;
        (*p_Fout).i = (*p_scratch_interleaved_a).i + (*p_scratch_interleaved_b).r;
        p_scratch_interleaved_a++;  p_scratch_interleaved_b++; p_Fout++;

        (*p_Fout).r = (*p_scratch_interleaved_a).r - (*p_scratch_interleaved_b).i;
        (*p_Fout).i = (*p_scratch_interleaved_a).i + (*p_scratch_interleaved_b).r;
        p_scratch_interleaved_a++;  p_scratch_interleaved_b++; p_Fout++;

        (*p_Fout).r = (*p_scratch_interleaved_a).r - (*p_scratch_interleaved_b).i;
        (*p_Fout).i = (*p_scratch_interleaved_a).i + (*p_scratch_interleaved_b).r;
        p_scratch_interleaved_a++;  p_scratch_interleaved_b++; p_Fout++;

        (*p_Fout).r = (*p_scratch_interleaved_a).r - (*p_scratch_interleaved_b).i;
        (*p_Fout).i = (*p_scratch_interleaved_a).i + (*p_scratch_interleaved_b).r;


    }
#else

    kiss_fft_cpx scratch1[6];
    kiss_fft_cpx scratch2[6];
    kiss_fft_cpx scratch3[6];
    kiss_fft_cpx scratch4[6];

    C_MUL(scratch1[0], Fout[4],      *(tw1+0));
    C_MUL(scratch2[0], Fout[4 + 1],  *(tw1 + fstride));
    C_MUL(scratch3[0], Fout[4 + 2],  *(tw1 + 2*fstride));
    C_MUL(scratch4[0], Fout[4 + 3],  *(tw1 + 3 * fstride));

    C_MUL(scratch1[1], Fout[8],      *(tw2 + 0));
    C_MUL(scratch2[1], Fout[8 + 1],  *(tw2 + 2*fstride));
    C_MUL(scratch3[1], Fout[8 + 2],  *(tw2 + 2* 2 * fstride));
    C_MUL(scratch4[1], Fout[8 + 3],  *(tw2 + 3 * 2 * fstride));

    C_MUL(scratch1[2], Fout[12],     *(tw3 + 0));
    C_MUL(scratch2[2], Fout[12 + 1], *(tw3 + 3*fstride));
    C_MUL(scratch3[2], Fout[12 + 2], *(tw3 + 2* 3 * fstride));
    C_MUL(scratch4[2], Fout[12 + 3], *(tw3 + 3 * 3 * fstride));

    C_SUB(scratch1[5], *Fout,   scratch1[1]);
    C_SUB(scratch2[5], Fout[1], scratch2[1]);
    C_SUB(scratch3[5], Fout[2], scratch3[1]);
    C_SUB(scratch4[5], Fout[3], scratch4[1]);

    C_ADDTO(Fout[0], scratch1[1]);
    C_ADDTO(Fout[1], scratch2[1]);
    C_ADDTO(Fout[2], scratch3[1]);
    C_ADDTO(Fout[3], scratch4[1]);

    C_ADD(scratch1[3], scratch1[0], scratch1[2]);
    C_ADD(scratch2[3], scratch2[0], scratch2[2]);
    C_ADD(scratch3[3], scratch3[0], scratch3[2]);
    C_ADD(scratch4[3], scratch4[0], scratch4[2]);

    C_SUB(scratch1[4], scratch1[0], scratch1[2]);
    C_SUB(scratch2[4], scratch2[0], scratch2[2]);
    C_SUB(scratch3[4], scratch3[0], scratch3[2]);
    C_SUB(scratch4[4], scratch4[0], scratch4[2]);

    C_SUB(Fout[8],     *Fout,   scratch1[3]);
    C_SUB(Fout[8 + 1], Fout[1], scratch2[3]);
    C_SUB(Fout[8 + 2], Fout[2], scratch3[3]);
    C_SUB(Fout[8 + 3], Fout[3], scratch4[3]);

    C_ADDTO(*Fout,   scratch1[3]);
    C_ADDTO(Fout[1], scratch2[3]);
    C_ADDTO(Fout[2], scratch3[3]);
    C_ADDTO(Fout[3], scratch4[3]);

    if (st->inverse) {

        Fout[4].r   = scratch1[5].r - scratch1[4].i;
        Fout[4].i   = scratch1[5].i + scratch1[4].r;
        Fout[4+1].r = scratch2[5].r - scratch2[4].i;
        Fout[4+1].i = scratch2[5].i + scratch2[4].r;
        Fout[4+2].r = scratch3[5].r - scratch3[4].i;
        Fout[4+2].i = scratch3[5].i + scratch3[4].r;
        Fout[4+3].r = scratch4[5].r - scratch4[4].i;
        Fout[4+3].i = scratch4[5].i + scratch4[4].r;

        Fout[12].r   = scratch1[5].r + scratch1[4].i;
        Fout[12].i   = scratch1[5].i - scratch1[4].r;
        Fout[12+1].r = scratch2[5].r + scratch2[4].i;
        Fout[12+1].i = scratch2[5].i - scratch2[4].r;
        Fout[12+2].r = scratch3[5].r + scratch3[4].i;
        Fout[12+2].i = scratch3[5].i - scratch3[4].r;
        Fout[12+3].r = scratch4[5].r + scratch4[4].i;
        Fout[12+3].i = scratch4[5].i - scratch4[4].r;

    }
    else {
        Fout[4].r   = scratch1[5].r + scratch1[4].i;
        Fout[4].i   = scratch1[5].i - scratch1[4].r;
        Fout[4+1].r = scratch2[5].r + scratch2[4].i;
        Fout[4+1].i = scratch2[5].i - scratch2[4].r;
        Fout[4+2].r = scratch3[5].r + scratch3[4].i;
        Fout[4+2].i = scratch3[5].i - scratch3[4].r;
        Fout[4+3].r = scratch4[5].r + scratch4[4].i;
        Fout[4+3].i = scratch4[5].i - scratch4[4].r;

        Fout[12].r   = scratch1[5].r - scratch1[4].i;
        Fout[12].i   = scratch1[5].i + scratch1[4].r;
        Fout[12+1].r = scratch2[5].r - scratch2[4].i;
        Fout[12+1].i = scratch2[5].i + scratch2[4].r;
        Fout[12+2].r = scratch3[5].r - scratch3[4].i;
        Fout[12+2].i = scratch3[5].i + scratch3[4].r;
        Fout[12+3].r = scratch4[5].r - scratch4[4].i;
        Fout[12+3].i = scratch4[5].i + scratch4[4].r;
    }
#endif

}





static void kf_bfly3(
         kiss_fft_cpx * Fout,
         const size_t fstride,
         const kiss_fft_cfg st,
         size_t m
         )
{
     size_t k=m;
     const size_t m2 = 2*m;
     kiss_fft_cpx *tw1,*tw2;
     kiss_fft_cpx scratch[5];
     kiss_fft_cpx epi3;
     epi3 = st->twiddles[fstride*m];

     tw1=tw2=st->twiddles;

     do{
         C_FIXDIV(*Fout,3); C_FIXDIV(Fout[m],3); C_FIXDIV(Fout[m2],3);

         C_MUL(scratch[1],Fout[m] , *tw1);
         C_MUL(scratch[2],Fout[m2] , *tw2);

         C_ADD(scratch[3],scratch[1],scratch[2]);
         C_SUB(scratch[0],scratch[1],scratch[2]);
         tw1 += fstride;
         tw2 += fstride*2;

         Fout[m].r = Fout->r - HALF_OF(scratch[3].r);
         Fout[m].i = Fout->i - HALF_OF(scratch[3].i);

         C_MULBYSCALAR( scratch[0] , epi3.i );

         C_ADDTO(*Fout,scratch[3]);

         Fout[m2].r = Fout[m].r + scratch[0].i;
         Fout[m2].i = Fout[m].i - scratch[0].r;

         Fout[m].r -= scratch[0].i;
         Fout[m].i += scratch[0].r;

         ++Fout;
     }while(--k);
}

static void kf_bfly5(
        kiss_fft_cpx * Fout,
        const size_t fstride,
        const kiss_fft_cfg st,
        int m
        )
{
    kiss_fft_cpx *Fout0,*Fout1,*Fout2,*Fout3,*Fout4;
    int u;
    kiss_fft_cpx scratch[13];
    kiss_fft_cpx * twiddles = st->twiddles;
    kiss_fft_cpx *tw;
    kiss_fft_cpx ya,yb;
    ya = twiddles[fstride*m];
    yb = twiddles[fstride*2*m];

    Fout0=Fout;
    Fout1=Fout0+m;
    Fout2=Fout0+2*m;
    Fout3=Fout0+3*m;
    Fout4=Fout0+4*m;

    tw=st->twiddles;
    for ( u=0; u<m; ++u ) {
        C_FIXDIV( *Fout0,5); C_FIXDIV( *Fout1,5); C_FIXDIV( *Fout2,5); C_FIXDIV( *Fout3,5); C_FIXDIV( *Fout4,5);
        scratch[0] = *Fout0;

        C_MUL(scratch[1] ,*Fout1, tw[u*fstride]);
        C_MUL(scratch[2] ,*Fout2, tw[2*u*fstride]);
        C_MUL(scratch[3] ,*Fout3, tw[3*u*fstride]);
        C_MUL(scratch[4] ,*Fout4, tw[4*u*fstride]);

        C_ADD( scratch[7],scratch[1],scratch[4]);
        C_SUB( scratch[10],scratch[1],scratch[4]);
        C_ADD( scratch[8],scratch[2],scratch[3]);
        C_SUB( scratch[9],scratch[2],scratch[3]);

        Fout0->r += scratch[7].r + scratch[8].r;
        Fout0->i += scratch[7].i + scratch[8].i;

        scratch[5].r = scratch[0].r + S_MUL(scratch[7].r,ya.r) + S_MUL(scratch[8].r,yb.r);
        scratch[5].i = scratch[0].i + S_MUL(scratch[7].i,ya.r) + S_MUL(scratch[8].i,yb.r);

        scratch[6].r =  S_MUL(scratch[10].i,ya.i) + S_MUL(scratch[9].i,yb.i);
        scratch[6].i = -S_MUL(scratch[10].r,ya.i) - S_MUL(scratch[9].r,yb.i);

        C_SUB(*Fout1,scratch[5],scratch[6]);
        C_ADD(*Fout4,scratch[5],scratch[6]);

        scratch[11].r = scratch[0].r + S_MUL(scratch[7].r,yb.r) + S_MUL(scratch[8].r,ya.r);
        scratch[11].i = scratch[0].i + S_MUL(scratch[7].i,yb.r) + S_MUL(scratch[8].i,ya.r);
        scratch[12].r = - S_MUL(scratch[10].i,yb.i) + S_MUL(scratch[9].i,ya.i);
        scratch[12].i = S_MUL(scratch[10].r,yb.i) - S_MUL(scratch[9].r,ya.i);

        C_ADD(*Fout2,scratch[11],scratch[12]);
        C_SUB(*Fout3,scratch[11],scratch[12]);

        ++Fout0;++Fout1;++Fout2;++Fout3;++Fout4;
    }
}

/* perform the butterfly for one stage of a mixed radix FFT */
static void kf_bfly_generic(
        kiss_fft_cpx * Fout,
        const size_t fstride,
        const kiss_fft_cfg st,
        int m,
        int p
        )
{
    int u,k,q1,q;
    kiss_fft_cpx * twiddles = st->twiddles;
    kiss_fft_cpx t;
    int Norig = st->nfft;

    kiss_fft_cpx * scratch = (kiss_fft_cpx*)KISS_FFT_TMP_ALLOC(sizeof(kiss_fft_cpx)*p);
    if (scratch == NULL){
        KISS_FFT_ERROR("Memory allocation failed.");
        return;
    }

    for ( u=0; u<m; ++u ) {
        k=u;
        for ( q1=0 ; q1<p ; ++q1 ) {
            scratch[q1] = Fout[ k  ];
            C_FIXDIV(scratch[q1],p);
            k += m;
        }

        k=u;
        for ( q1=0 ; q1<p ; ++q1 ) {
            int twidx=0;
            Fout[ k ] = scratch[0];
            for (q=1;q<p;++q ) {
                twidx += fstride * k;
                if (twidx>=Norig) twidx-=Norig;
                C_MUL(t,scratch[q] , twiddles[twidx] );
                C_ADDTO( Fout[ k ] ,t);
            }
            k += m;
        }
    }
    KISS_FFT_TMP_FREE(scratch);
}


static
void kf_work(
        kiss_fft_cpx * Fout,
        const kiss_fft_cpx * f,
        const size_t fstride,
        int in_stride,
        int * factors,
        const kiss_fft_cfg st
        )
{
    kiss_fft_cpx * Fout_beg=Fout;
    const int p=*factors++; /* the radix  */
    const int m=*factors++; /* stage's fft length/p */
    const kiss_fft_cpx * Fout_end = Fout + p*m;

#ifdef _OPENMP
    // use openmp extensions at the
    // top-level (not recursive)
    if (fstride==1 && p<=5 && m!=1)
    {
        int k;

        // execute the p different work units in different threads
#       pragma omp parallel for
        for (k=0;k<p;++k)
            kf_work( Fout +k*m, f+ fstride*in_stride*k,fstride*p,in_stride,factors,st);
        // all threads have joined by this point

        switch (p) {
            case 2: kf_bfly2(Fout,fstride,st,m); break;
            case 3: kf_bfly3(Fout,fstride,st,m); break;
            case 4: kf_bfly4(Fout,fstride,st,m); break;
            case 5: kf_bfly5(Fout,fstride,st,m); break;
            default: kf_bfly_generic(Fout,fstride,st,m,p); break;
        }
        return;
    }
#endif

    if (m==1) {
        do{
            *Fout = *f;
            f += fstride*in_stride;
        }while(++Fout != Fout_end );
    }else{
        do{
            // recursive call:
            // DFT of size m*p performed by doing
            // p instances of smaller DFTs of size m,
            // each one takes a decimated version of the input
            kf_work( Fout , f, fstride*p, in_stride, factors,st);
            f += fstride*in_stride;
        }while( (Fout += m) != Fout_end );
    }

    Fout=Fout_beg;

    // recombine the p smaller DFTs
    switch (p) {
        case 2: kf_bfly2(Fout,fstride,st,m); break;
        case 3: kf_bfly3(Fout,fstride,st,m); break;
        case 4: 
            if ((m == 4) && (fstride==16))
            {
                
                
                kf_bfly4_m4_multiple_scratches(Fout, fstride, st);

                // kiss_fft_cpx* fast_Fout = (kiss_fft_cpx*)malloc(16 * sizeof(kiss_fft_cpx));
                // f_map_Fout_to_fast_Fout(Fout, fast_Fout);
                // kf_bfly4_m4_remapped(fast_Fout, fstride, st);
                // f_map_fast_Fout_to_Fout(fast_Fout, Fout);
               
            }
            else
            {
                kf_bfly4(Fout, fstride, st, m);
            }
            break;
        case 5: kf_bfly5(Fout,fstride,st,m); break;
        default: kf_bfly_generic(Fout,fstride,st,m,p); break;
    }
}

/*  facbuf is populated by p1,m1,p2,m2, ...
    where
    p[i] * m[i] = m[i-1]
    m0 = n                  */
static
void kf_factor(int n,int * facbuf)
{
    int p=4;
    double floor_sqrt;
    floor_sqrt = floor( sqrt((double)n) );

    /*factor out powers of 4, powers of 2, then any remaining primes */
    do {
        while (n % p) {
            switch (p) {
                case 4: p = 2; break;
                case 2: p = 3; break;
                default: p += 2; break;
            }
            if (p > floor_sqrt)
                p = n;          /* no more factors, skip to end */
        }
        n /= p;
        *facbuf++ = p;
        *facbuf++ = n;
    } while (n > 1);
}

/*
 *
 * User-callable function to allocate all necessary storage space for the fft.
 *
 * The return value is a contiguous block of memory, allocated with malloc.  As such,
 * It can be freed with free(), rather than a kiss_fft-specific function.
 * */
kiss_fft_cfg kiss_fft_alloc(int nfft,int inverse_fft,void * mem,size_t * lenmem )
{
    KISS_FFT_ALIGN_CHECK(mem)

    kiss_fft_cfg st=NULL;
    size_t memneeded = KISS_FFT_ALIGN_SIZE_UP(sizeof(struct kiss_fft_state)
        + sizeof(kiss_fft_cpx)*(nfft-1)); /* twiddle factors*/

    if ( lenmem==NULL ) {
        st = ( kiss_fft_cfg)KISS_FFT_MALLOC( memneeded );
    }else{
        if (mem != NULL && *lenmem >= memneeded)
            st = (kiss_fft_cfg)mem;
        *lenmem = memneeded;
    }
    if (st) {
        int i;
        st->nfft=nfft;
        st->inverse = inverse_fft;

        for (i=0;i<nfft;++i) {
            const double pi=3.141592653589793238462643383279502884197169399375105820974944;
            double phase = -2*pi*i / nfft;
            if (st->inverse)
                phase *= -1;
            kf_cexp(st->twiddles+i, phase );
        }

        kf_factor(nfft,st->factors);
    }
    return st;
}


void kiss_fft_stride(kiss_fft_cfg st,const kiss_fft_cpx *fin,kiss_fft_cpx *fout,int in_stride)
{
    if (fin == fout) {
        //NOTE: this is not really an in-place FFT algorithm.
        //It just performs an out-of-place FFT into a temp buffer
        if (fout == NULL){
            KISS_FFT_ERROR("fout buffer NULL.");
        return;
        }

        kiss_fft_cpx * tmpbuf = (kiss_fft_cpx*)KISS_FFT_TMP_ALLOC( sizeof(kiss_fft_cpx)*st->nfft);
        if (tmpbuf == NULL){
            KISS_FFT_ERROR("Memory allocation error.");
        return;
        }



        kf_work(tmpbuf,fin,1,in_stride, st->factors,st);
        memcpy(fout,tmpbuf,sizeof(kiss_fft_cpx)*st->nfft);
        KISS_FFT_TMP_FREE(tmpbuf);
    }else{
        kf_work( fout, fin, 1,in_stride, st->factors,st );
    }
}

void kiss_fft(kiss_fft_cfg cfg,const kiss_fft_cpx *fin,kiss_fft_cpx *fout)
{
    kiss_fft_stride(cfg,fin,fout,1);
}


void kiss_fft_cleanup(void)
{
    // nothing needed any more
}

int kiss_fft_next_fast_size(int n)
{
    while(1) {
        int m=n;
        while ( (m%2) == 0 ) m/=2;
        while ( (m%3) == 0 ) m/=3;
        while ( (m%5) == 0 ) m/=5;
        if (m<=1)
            break; /* n is completely factorable by twos, threes, and fives */
        n++;
    }
    return n;
}
