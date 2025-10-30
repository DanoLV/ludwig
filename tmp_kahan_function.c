/*****************************************************************************
 *
 *  d3q19_mode2f_chunk_kahan
 *
 *  New version with Kahan compensated summation to reduce roundoff error
 *  propagation in long simulations.
 *
 *****************************************************************************/

__device__ void d3q19_mode2f_chunk_kahan(double* mode, double* fchunk) {

  double ftmp[NSIMDVL];
  double ftmp_c[NSIMDVL];  /* Kahan compensation */

  int iv;

  /* Helper macro for Kahan summation */
  #define KAHAN_ADD(coeff, mode_idx) \
    for_simd_v(iv, NSIMDVL) { \
      volatile double val = (coeff) * mode[(mode_idx) * NSIMDVL + iv]; \
      volatile double y = val - ftmp_c[iv]; \
      volatile double t = ftmp[iv] + y; \
      ftmp_c[iv] = (t - ftmp[iv]) - y; \
      ftmp[iv] = t; \
    }

  /* p=0 */
  for_simd_v(iv, NSIMDVL) { ftmp[iv] = 0.0; ftmp_c[iv] = 0.0; }
  KAHAN_ADD(w0, 0);
  KAHAN_ADD(c0, 1);
  KAHAN_ADD(c0, 2);
  KAHAN_ADD(c0, 3);
  KAHAN_ADD(-r2, 4);
  KAHAN_ADD(c0, 5);
  KAHAN_ADD(c0, 6);
  KAHAN_ADD(-r2, 7);
  KAHAN_ADD(c0, 8);
  KAHAN_ADD(-r2, 9);
  KAHAN_ADD(c0, 10);
  KAHAN_ADD(c0, 11);
  KAHAN_ADD(c0, 12);
  KAHAN_ADD(c0, 13);
  KAHAN_ADD(c0, 14);
  KAHAN_ADD(c0, 15);
  KAHAN_ADD(c0, 16);
  KAHAN_ADD(c0, 17);
  KAHAN_ADD(r6, 18);
  for_simd_v(iv, NSIMDVL) fchunk[0 * NSIMDVL + iv] = ftmp[iv];

  /* p=1 */
  for_simd_v(iv, NSIMDVL) { ftmp[iv] = 0.0; ftmp_c[iv] = 0.0; }
  KAHAN_ADD(w2, 0);
  KAHAN_ADD(wa, 1);
  KAHAN_ADD(wa, 2);
  KAHAN_ADD(c0, 3);
  KAHAN_ADD(wa, 4);
  KAHAN_ADD(r4, 5);
  KAHAN_ADD(c0, 6);
  KAHAN_ADD(wa, 7);
  KAHAN_ADD(c0, 8);
  KAHAN_ADD(-wb, 9);
  KAHAN_ADD(-wb, 10);
  KAHAN_ADD(-wa, 11);
  KAHAN_ADD(-wa, 12);
  KAHAN_ADD(c0, 13);
  KAHAN_ADD(c0, 14);
  KAHAN_ADD(c0, 15);
  KAHAN_ADD(c0, 16);
  KAHAN_ADD(c0, 17);
  KAHAN_ADD(wc, 18);
  for_simd_v(iv, NSIMDVL) fchunk[1 * NSIMDVL + iv] = ftmp[iv];

  /* p=2 */
  for_simd_v(iv, NSIMDVL) { ftmp[iv] = 0.0; ftmp_c[iv] = 0.0; }
  KAHAN_ADD(w2, 0);
  KAHAN_ADD(wa, 1);
  KAHAN_ADD(c0, 2);
  KAHAN_ADD(wa, 3);
  KAHAN_ADD(wa, 4);
  KAHAN_ADD(c0, 5);
  KAHAN_ADD(r4, 6);
  KAHAN_ADD(-wb, 7);
  KAHAN_ADD(c0, 8);
  KAHAN_ADD(wa, 9);
  KAHAN_ADD(wd, 10);
  KAHAN_ADD(wb, 11);
  KAHAN_ADD(c0, 12);
  KAHAN_ADD(wb, 13);
  KAHAN_ADD(-we, 14);
  KAHAN_ADD(-r8, 15);
  KAHAN_ADD(c0, 16);
  KAHAN_ADD(-r8, 17);
  KAHAN_ADD(wc, 18);
  for_simd_v(iv, NSIMDVL) fchunk[2 * NSIMDVL + iv] = ftmp[iv];

  /* p=3 */
  for_simd_v(iv, NSIMDVL) { ftmp[iv] = 0.0; ftmp_c[iv] = 0.0; }
  KAHAN_ADD(w1, 0);
  KAHAN_ADD(r6, 1);
  KAHAN_ADD(c0, 2);
  KAHAN_ADD(c0, 3);
  KAHAN_ADD(r6, 4);
  KAHAN_ADD(c0, 5);
  KAHAN_ADD(c0, 6);
  KAHAN_ADD(-wa, 7);
  KAHAN_ADD(c0, 8);
  KAHAN_ADD(-wa, 9);
  KAHAN_ADD(wb, 10);
  KAHAN_ADD(wa, 11);
  KAHAN_ADD(c0, 12);
  KAHAN_ADD(c0, 13);
  KAHAN_ADD(r8, 14);
  KAHAN_ADD(r4, 15);
  KAHAN_ADD(c0, 16);
  KAHAN_ADD(c0, 17);
  KAHAN_ADD(-w1, 18);
  for_simd_v(iv, NSIMDVL) fchunk[3 * NSIMDVL + iv] = ftmp[iv];

  /* p=4 */
  for_simd_v(iv, NSIMDVL) { ftmp[iv] = 0.0; ftmp_c[iv] = 0.0; }
  KAHAN_ADD(w2, 0);
  KAHAN_ADD(wa, 1);
  KAHAN_ADD(c0, 2);
  KAHAN_ADD(-wa, 3);
  KAHAN_ADD(wa, 4);
  KAHAN_ADD(c0, 5);
  KAHAN_ADD(-r4, 6);
  KAHAN_ADD(-wb, 7);
  KAHAN_ADD(c0, 8);
  KAHAN_ADD(wa, 9);
  KAHAN_ADD(wd, 10);
  KAHAN_ADD(wb, 11);
  KAHAN_ADD(c0, 12);
  KAHAN_ADD(-wb, 13);
  KAHAN_ADD(-we, 14);
  KAHAN_ADD(-r8, 15);
  KAHAN_ADD(c0, 16);
  KAHAN_ADD(r8, 17);
  KAHAN_ADD(wc, 18);
  for_simd_v(iv, NSIMDVL) fchunk[4 * NSIMDVL + iv] = ftmp[iv];

  /* p=5 */
  for_simd_v(iv, NSIMDVL) { ftmp[iv] = 0.0; ftmp_c[iv] = 0.0; }
  KAHAN_ADD(w2, 0);
  KAHAN_ADD(wa, 1);
  KAHAN_ADD(-wa, 2);
  KAHAN_ADD(c0, 3);
  KAHAN_ADD(wa, 4);
  KAHAN_ADD(-r4, 5);
  KAHAN_ADD(c0, 6);
  KAHAN_ADD(wa, 7);
  KAHAN_ADD(c0, 8);
  KAHAN_ADD(-wb, 9);
  KAHAN_ADD(-wb, 10);
  KAHAN_ADD(-wa, 11);
  KAHAN_ADD(wa, 12);
  KAHAN_ADD(c0, 13);
  KAHAN_ADD(c0, 14);
  KAHAN_ADD(c0, 15);
  KAHAN_ADD(c0, 16);
  KAHAN_ADD(c0, 17);
  KAHAN_ADD(wc, 18);
  for_simd_v(iv, NSIMDVL) fchunk[5 * NSIMDVL + iv] = ftmp[iv];

  /* p=6 */
  for_simd_v(iv, NSIMDVL) { ftmp[iv] = 0.0; ftmp_c[iv] = 0.0; }
  KAHAN_ADD(w2, 0);
  KAHAN_ADD(c0, 1);
  KAHAN_ADD(wa, 2);
  KAHAN_ADD(wa, 3);
  KAHAN_ADD(-wb, 4);
  KAHAN_ADD(c0, 5);
  KAHAN_ADD(c0, 6);
  KAHAN_ADD(wa, 7);
  KAHAN_ADD(r4, 8);
  KAHAN_ADD(wa, 9);
  KAHAN_ADD(wd, 10);
  KAHAN_ADD(c0, 11);
  KAHAN_ADD(wb, 12);
  KAHAN_ADD(wb, 13);
  KAHAN_ADD(we, 14);
  KAHAN_ADD(c0, 15);
  KAHAN_ADD(r8, 16);
  KAHAN_ADD(r8, 17);
  KAHAN_ADD(wc, 18);
  for_simd_v(iv, NSIMDVL) fchunk[6 * NSIMDVL + iv] = ftmp[iv];

  /* p=7 */
  for_simd_v(iv, NSIMDVL) { ftmp[iv] = 0.0; ftmp_c[iv] = 0.0; }
  KAHAN_ADD(w1, 0);
  KAHAN_ADD(c0, 1);
  KAHAN_ADD(r6, 2);
  KAHAN_ADD(c0, 3);
  KAHAN_ADD(-wa, 4);
  KAHAN_ADD(c0, 5);
  KAHAN_ADD(c0, 6);
  KAHAN_ADD(r6, 7);
  KAHAN_ADD(c0, 8);
  KAHAN_ADD(-wa, 9);
  KAHAN_ADD(wb, 10);
  KAHAN_ADD(c0, 11);
  KAHAN_ADD(wa, 12);
  KAHAN_ADD(c0, 13);
  KAHAN_ADD(-r8, 14);
  KAHAN_ADD(c0, 15);
  KAHAN_ADD(-r4, 16);
  KAHAN_ADD(c0, 17);
  KAHAN_ADD(-w1, 18);
  for_simd_v(iv, NSIMDVL) fchunk[7 * NSIMDVL + iv] = ftmp[iv];

  /* p=8 */
  for_simd_v(iv, NSIMDVL) { ftmp[iv] = 0.0; ftmp_c[iv] = 0.0; }
  KAHAN_ADD(w2, 0);
  KAHAN_ADD(c0, 1);
  KAHAN_ADD(wa, 2);
  KAHAN_ADD(-wa, 3);
  KAHAN_ADD(-wb, 4);
  KAHAN_ADD(c0, 5);
  KAHAN_ADD(c0, 6);
  KAHAN_ADD(wa, 7);
  KAHAN_ADD(-r4, 8);
  KAHAN_ADD(wa, 9);
  KAHAN_ADD(wd, 10);
  KAHAN_ADD(c0, 11);
  KAHAN_ADD(wb, 12);
  KAHAN_ADD(-wb, 13);
  KAHAN_ADD(we, 14);
  KAHAN_ADD(c0, 15);
  KAHAN_ADD(r8, 16);
  KAHAN_ADD(-r8, 17);
  KAHAN_ADD(wc, 18);
  for_simd_v(iv, NSIMDVL) fchunk[8 * NSIMDVL + iv] = ftmp[iv];

  /* p=9 */
  for_simd_v(iv, NSIMDVL) { ftmp[iv] = 0.0; ftmp_c[iv] = 0.0; }
  KAHAN_ADD(w1, 0);
  KAHAN_ADD(c0, 1);
  KAHAN_ADD(c0, 2);
  KAHAN_ADD(r6, 3);
  KAHAN_ADD(-wa, 4);
  KAHAN_ADD(c0, 5);
  KAHAN_ADD(c0, 6);
  KAHAN_ADD(-wa, 7);
  KAHAN_ADD(c0, 8);
  KAHAN_ADD(r6, 9);
  KAHAN_ADD(-wa, 10);
  KAHAN_ADD(c0, 11);
  KAHAN_ADD(c0, 12);
  KAHAN_ADD(-r6, 13);
  KAHAN_ADD(c0, 14);
  KAHAN_ADD(c0, 15);
  KAHAN_ADD(c0, 16);
  KAHAN_ADD(c0, 17);
  KAHAN_ADD(-w1, 18);
  for_simd_v(iv, NSIMDVL) fchunk[9 * NSIMDVL + iv] = ftmp[iv];

  /* p=10 */
  for_simd_v(iv, NSIMDVL) { ftmp[iv] = 0.0; ftmp_c[iv] = 0.0; }
  KAHAN_ADD(w1, 0);
  KAHAN_ADD(c0, 1);
  KAHAN_ADD(c0, 2);
  KAHAN_ADD(-r6, 3);
  KAHAN_ADD(-wa, 4);
  KAHAN_ADD(c0, 5);
  KAHAN_ADD(c0, 6);
  KAHAN_ADD(-wa, 7);
  KAHAN_ADD(c0, 8);
  KAHAN_ADD(r6, 9);
  KAHAN_ADD(-wa, 10);
  KAHAN_ADD(c0, 11);
  KAHAN_ADD(c0, 12);
  KAHAN_ADD(r6, 13);
  KAHAN_ADD(c0, 14);
  KAHAN_ADD(c0, 15);
  KAHAN_ADD(c0, 16);
  KAHAN_ADD(c0, 17);
  KAHAN_ADD(-w1, 18);
  for_simd_v(iv, NSIMDVL) fchunk[10 * NSIMDVL + iv] = ftmp[iv];

  /* p=11 */
  for_simd_v(iv, NSIMDVL) { ftmp[iv] = 0.0; ftmp_c[iv] = 0.0; }
  KAHAN_ADD(w2, 0);
  KAHAN_ADD(c0, 1);
  KAHAN_ADD(-wa, 2);
  KAHAN_ADD(wa, 3);
  KAHAN_ADD(-wb, 4);
  KAHAN_ADD(c0, 5);
  KAHAN_ADD(c0, 6);
  KAHAN_ADD(wa, 7);
  KAHAN_ADD(-r4, 8);
  KAHAN_ADD(wa, 9);
  KAHAN_ADD(wd, 10);
  KAHAN_ADD(c0, 11);
  KAHAN_ADD(-wb, 12);
  KAHAN_ADD(wb, 13);
  KAHAN_ADD(we, 14);
  KAHAN_ADD(c0, 15);
  KAHAN_ADD(-r8, 16);
  KAHAN_ADD(r8, 17);
  KAHAN_ADD(wc, 18);
  for_simd_v(iv, NSIMDVL) fchunk[11 * NSIMDVL + iv] = ftmp[iv];

  /* p=12 */
  for_simd_v(iv, NSIMDVL) { ftmp[iv] = 0.0; ftmp_c[iv] = 0.0; }
  KAHAN_ADD(w1, 0);
  KAHAN_ADD(c0, 1);
  KAHAN_ADD(-r6, 2);
  KAHAN_ADD(c0, 3);
  KAHAN_ADD(-wa, 4);
  KAHAN_ADD(c0, 5);
  KAHAN_ADD(c0, 6);
  KAHAN_ADD(r6, 7);
  KAHAN_ADD(c0, 8);
  KAHAN_ADD(-wa, 9);
  KAHAN_ADD(wb, 10);
  KAHAN_ADD(c0, 11);
  KAHAN_ADD(-wa, 12);
  KAHAN_ADD(c0, 13);
  KAHAN_ADD(-r8, 14);
  KAHAN_ADD(c0, 15);
  KAHAN_ADD(r4, 16);
  KAHAN_ADD(c0, 17);
  KAHAN_ADD(-w1, 18);
  for_simd_v(iv, NSIMDVL) fchunk[12 * NSIMDVL + iv] = ftmp[iv];

  /* p=13 */
  for_simd_v(iv, NSIMDVL) { ftmp[iv] = 0.0; ftmp_c[iv] = 0.0; }
  KAHAN_ADD(w2, 0);
  KAHAN_ADD(c0, 1);
  KAHAN_ADD(-wa, 2);
  KAHAN_ADD(-wa, 3);
  KAHAN_ADD(-wb, 4);
  KAHAN_ADD(c0, 5);
  KAHAN_ADD(c0, 6);
  KAHAN_ADD(wa, 7);
  KAHAN_ADD(r4, 8);
  KAHAN_ADD(wa, 9);
  KAHAN_ADD(wd, 10);
  KAHAN_ADD(c0, 11);
  KAHAN_ADD(-wb, 12);
  KAHAN_ADD(-wb, 13);
  KAHAN_ADD(we, 14);
  KAHAN_ADD(c0, 15);
  KAHAN_ADD(-r8, 16);
  KAHAN_ADD(-r8, 17);
  KAHAN_ADD(wc, 18);
  for_simd_v(iv, NSIMDVL) fchunk[13 * NSIMDVL + iv] = ftmp[iv];

  /* p=14 */
  for_simd_v(iv, NSIMDVL) { ftmp[iv] = 0.0; ftmp_c[iv] = 0.0; }
  KAHAN_ADD(w2, 0);
  KAHAN_ADD(-wa, 1);
  KAHAN_ADD(wa, 2);
  KAHAN_ADD(c0, 3);
  KAHAN_ADD(wa, 4);
  KAHAN_ADD(-r4, 5);
  KAHAN_ADD(c0, 6);
  KAHAN_ADD(wa, 7);
  KAHAN_ADD(c0, 8);
  KAHAN_ADD(-wb, 9);
  KAHAN_ADD(-wb, 10);
  KAHAN_ADD(wa, 11);
  KAHAN_ADD(-wa, 12);
  KAHAN_ADD(c0, 13);
  KAHAN_ADD(c0, 14);
  KAHAN_ADD(c0, 15);
  KAHAN_ADD(c0, 16);
  KAHAN_ADD(c0, 17);
  KAHAN_ADD(wc, 18);
  for_simd_v(iv, NSIMDVL) fchunk[14 * NSIMDVL + iv] = ftmp[iv];

  /* p=15 */
  for_simd_v(iv, NSIMDVL) { ftmp[iv] = 0.0; ftmp_c[iv] = 0.0; }
  KAHAN_ADD(w2, 0);
  KAHAN_ADD(-wa, 1);
  KAHAN_ADD(c0, 2);
  KAHAN_ADD(wa, 3);
  KAHAN_ADD(wa, 4);
  KAHAN_ADD(c0, 5);
  KAHAN_ADD(-r4, 6);
  KAHAN_ADD(-wb, 7);
  KAHAN_ADD(c0, 8);
  KAHAN_ADD(wa, 9);
  KAHAN_ADD(wd, 10);
  KAHAN_ADD(-wb, 11);
  KAHAN_ADD(c0, 12);
  KAHAN_ADD(wb, 13);
  KAHAN_ADD(-we, 14);
  KAHAN_ADD(r8, 15);
  KAHAN_ADD(c0, 16);
  KAHAN_ADD(-r8, 17);
  KAHAN_ADD(wc, 18);
  for_simd_v(iv, NSIMDVL) fchunk[15 * NSIMDVL + iv] = ftmp[iv];

  /* p=16 */
  for_simd_v(iv, NSIMDVL) { ftmp[iv] = 0.0; ftmp_c[iv] = 0.0; }
  KAHAN_ADD(w1, 0);
  KAHAN_ADD(-r6, 1);
  KAHAN_ADD(c0, 2);
  KAHAN_ADD(c0, 3);
  KAHAN_ADD(r6, 4);
  KAHAN_ADD(c0, 5);
  KAHAN_ADD(c0, 6);
  KAHAN_ADD(-wa, 7);
  KAHAN_ADD(c0, 8);
  KAHAN_ADD(-wa, 9);
  KAHAN_ADD(wb, 10);
  KAHAN_ADD(-wa, 11);
  KAHAN_ADD(c0, 12);
  KAHAN_ADD(c0, 13);
  KAHAN_ADD(r8, 14);
  KAHAN_ADD(-r4, 15);
  KAHAN_ADD(c0, 16);
  KAHAN_ADD(c0, 17);
  KAHAN_ADD(-w1, 18);
  for_simd_v(iv, NSIMDVL) fchunk[16 * NSIMDVL + iv] = ftmp[iv];

  /* p=17 */
  for_simd_v(iv, NSIMDVL) { ftmp[iv] = 0.0; ftmp_c[iv] = 0.0; }
  KAHAN_ADD(w2, 0);
  KAHAN_ADD(-wa, 1);
  KAHAN_ADD(c0, 2);
  KAHAN_ADD(-wa, 3);
  KAHAN_ADD(wa, 4);
  KAHAN_ADD(c0, 5);
  KAHAN_ADD(r4, 6);
  KAHAN_ADD(-wb, 7);
  KAHAN_ADD(c0, 8);
  KAHAN_ADD(wa, 9);
  KAHAN_ADD(wd, 10);
  KAHAN_ADD(-wb, 11);
  KAHAN_ADD(c0, 12);
  KAHAN_ADD(-wb, 13);
  KAHAN_ADD(-we, 14);
  KAHAN_ADD(r8, 15);
  KAHAN_ADD(c0, 16);
  KAHAN_ADD(r8, 17);
  KAHAN_ADD(wc, 18);
  for_simd_v(iv, NSIMDVL) fchunk[17 * NSIMDVL + iv] = ftmp[iv];

  /* p=18 */
  for_simd_v(iv, NSIMDVL) { ftmp[iv] = 0.0; ftmp_c[iv] = 0.0; }
  KAHAN_ADD(w2, 0);
  KAHAN_ADD(-wa, 1);
  KAHAN_ADD(-wa, 2);
  KAHAN_ADD(c0, 3);
  KAHAN_ADD(wa, 4);
  KAHAN_ADD(r4, 5);
  KAHAN_ADD(c0, 6);
  KAHAN_ADD(wa, 7);
  KAHAN_ADD(c0, 8);
  KAHAN_ADD(-wb, 9);
  KAHAN_ADD(-wb, 10);
  KAHAN_ADD(wa, 11);
  KAHAN_ADD(wa, 12);
  KAHAN_ADD(c0, 13);
  KAHAN_ADD(c0, 14);
  KAHAN_ADD(c0, 15);
  KAHAN_ADD(c0, 16);
  KAHAN_ADD(c0, 17);
  KAHAN_ADD(wc, 18);
  for_simd_v(iv, NSIMDVL) fchunk[18 * NSIMDVL + iv] = ftmp[iv];

  #undef KAHAN_ADD
}
