// From bcftools,  bcftools/bam2bcf.c b46e7cd59811f87fded9f15db9416853034e74c6
// modified to remove assert calls

/*  bam2bcf.c -- variant calling.
    Copyright (C) 2010-2012 Broad Institute.
    Copyright (C) 2012-2022 Genome Research Ltd.
    Author: Heng Li <lh3@sanger.ac.uk>
Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:
The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
DEALINGS IN THE SOFTWARE.  */

#include <math.h>
#include <float.h>
#include <htslib/sam.h>
#include <htslib/kfunc.h>

/*
 *  calc_vdb() - returns value between zero (most biased) and one (no bias)
 *               on success, or HUGE_VAL when VDB cannot be calculated because
 *               of insufficient depth (<2x)
 *
 *  Variant Distance Bias tests if the variant bases are positioned within the
 *  reads with sufficient randomness. Unlike other tests, it looks only at
 *  variant reads and therefore gives different kind of information than Read
 *  Position Bias for instance. VDB was developed for detecting artefacts in
 *  RNA-seq calls where reads from spliced transcripts span splice site
 *  boundaries.  The current implementation differs somewhat from the original
 *  version described in supplementary material of PMID:22524474, but the idea
 *  remains the same. (Here the random variable tested is the average distance
 *  from the averaged position, not the average pairwise distance.)
 *
 *  For coverage of 2x, the calculation is exact but is approximated for the
 *  rest. The result is most accurate between 4-200x. For 3x or >200x, the
 *  reported values are slightly more favourable than those of a true random
 *  distribution.
 */


double calc_vdb(int *pos, int npos)
{
    // Note well: the parameters were obtained by fitting to simulated data of
    // 100bp reads. This assumes rescaling to 100bp in bcf_call_glfgen().
    const int readlen = 100;

    #define nparam 15
    const float param[nparam][3] = { {3,0.079,18}, {4,0.09,19.8}, {5,0.1,20.5}, {6,0.11,21.5},
        {7,0.125,21.6}, {8,0.135,22}, {9,0.14,22.2}, {10,0.153,22.3}, {15,0.19,22.8},
        {20,0.22,23.2}, {30,0.26,23.4}, {40,0.29,23.5}, {50,0.35,23.65}, {100,0.5,23.7},
        {200,0.7,23.7} };

    int i, dp = 0;
    float mean_pos = 0, mean_diff = 0;
    for (i=0; i<npos; i++)
    {
        if ( !pos[i] ) continue;
        dp += pos[i];
        mean_pos += pos[i]*i;
    }
    if ( dp<2 ) return HUGE_VAL;     // one or zero reads can be placed anywhere

    mean_pos /= dp;
    for (i=0; i<npos; i++)
    {
        if ( !pos[i] ) continue;
        mean_diff += pos[i] * fabs(i - mean_pos);
    }
    mean_diff /= dp;

    int ipos = mean_diff;   // tuned for float-to-int implicit conversion
    if ( dp==2 )
        return (2*readlen-2*(ipos+1)-1)*(ipos+1)/(readlen-1)/(readlen*0.5);

    if ( dp>=200 )
        i = nparam; // shortcut for big depths
    else
    {
        for (i=0; i<nparam; i++)
            if ( param[i][0]>=dp ) break;
    }
    float pshift, pscale;
    if ( i==nparam )
    {
        // the depth is too high, go with 200x
        pscale = param[nparam-1][1];
        pshift = param[nparam-1][2];
    }
    else if ( i>0 && param[i][0]!=dp )
    {
        // linear interpolation of parameters
        pscale = (param[i-1][1] + param[i][1])*0.5;
        pshift = (param[i-1][2] + param[i][2])*0.5;
    }
    else
    {
        pscale = param[i][1];
        pshift = param[i][2];
    }
    return 0.5*kf_erfc(-(mean_diff-pshift)*pscale);
}



// bam2bcf.c:596-612

static double mann_whitney_1947_(int n, int m, int U)
{
     if (U<0) return 0;
     if (n==0||m==0) return U==0 ? 1 : 0;
    return (double)n/(n+m)*mann_whitney_1947_(n-1,m,U-m) + (double)m/(n+m)*mann_whitney_1947_(n,m-1,U);
}

double mann_whitney_1947(int n, int m, int U)
{
    #include "mw.h"

    return (n < 8 && m < 8 && U < 50)
        ? mw[n-2][m-2][U]
        : mann_whitney_1947_(n,m,U);
}

// bam2bcf.c:725-793

// A Z-score version of the above function.
//
// See "Normal approximation and tie correction" at
// https://en.wikipedia.org/wiki/Mann%E2%80%93Whitney_U_test
//
// The Z score is the number of standard deviations above or below the mean
// with 0 being equality of the two distributions and +ve/-ve from there.
//
// This is a more robust score to filter on.
double calc_mwu_biasZ(int *a, int *b, int n, int left_only, int do_Z) {
    int i;
    int64_t t;

    // Optimisation
    for (i = 0; i < n; i++)
        if (b[i])
            break;
    int b_empty = (i == n);

    // Count equal (e), less-than (l) and greater-than (g) permutations.
    int e = 0, l = 0, na = 0, nb = 0;
    if (b_empty) {
        for (t = 0, i = n-1; i >= 0; i--) {
            na += a[i];
            t += (a[i]*a[i]-1)*a[i];  // adjustment score for ties
        }
    } else {
        for (t = 0, i = n-1; i >= 0; i--) {
            // Combinations of a[i] and b[j] for i==j
            e += a[i]*b[i];

            // nb is running total of b[i+1]..b[n-1].
            // Therefore a[i]*nb is the number of combinations of a[i] and b[j]
            // for all i < j.
            l += a[i]*nb;    // a<b

            na += a[i];
            nb += b[i];
            int p = a[i]+b[i];
            t += (p*p-1)*p;  // adjustment score for ties
        }
    }

    if (!na || !nb)
        return HUGE_VAL;

    double U, m;
    U = l + e*0.5; // Mann-Whitney U score
    m = na*nb / 2.0;

    // With ties adjustment
    double var2 = (na*nb)/12.0 * ((na+nb+1) - t/(double)((na+nb)*(na+nb-1)));
    // var = na*nb*(na+nb+1)/12.0; // simpler; minus tie adjustment
    if (var2 <= 0)
        return do_Z ? 0 : 1;

    if (do_Z) {
        // S.D. normalised Z-score
        //Z = (U - m - (U-m >= 0 ? 0.5 : -0.5)) / sd; // gatk method?
        return (U - m) / sqrt(var2);
    }

    // Else U score, which can be asymmetric for some data types.
    if (left_only && U > m)
        return HUGE_VAL; // one-sided, +ve bias is OK, -ve is not.

    if (na >= 8 || nb >= 8) {
        // Normal approximation, very good for na>=8 && nb>=8 and
        // reasonable if na<8 or nb<8
        return exp(-0.5*(U-m)*(U-m)/var2);
    }

    // Exact calculation
    if (na==1 || nb == 1)
        return mann_whitney_1947_(na, nb, U) * sqrt(2*M_PI*var2);
    else
        return mann_whitney_1947(na, nb, U) * sqrt(2*M_PI*var2);
}

