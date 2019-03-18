/* Downsample.cs
 *
 * Copyright 2017 Outloud Oy
 * mikko.k.koivisto@outloud.fi
 * 
 * This C# implementation is derived from the C-language implementation of RAPT distributed 
 * as part of the ESPS toolkit downloaded from http://www.speech.kth.se/speech/esps/esps.zip
 * ESPS toolkit has a BSD style license (see below)
 * 
 ******* BSD.txt in ESPS toolkit *******
 * Copyright (c) 2002 Microsoft Corp.
 *
 * The following terms apply to all files associated
 * with the software unless explicitly disclaimed in individual files.
 *
 * The authors hereby grant permission to use, copy, modify, distribute,
 * and license this software and its documentation for any purpose, provided
 * that existing copyright notices are retained in all copies and that this
 * notice is included verbatim in any distributions. No written agreement,
 * license, or royalty fee is required for any of the authorized uses.
 * Modifications to this software may be copyrighted by their authors
 * and need not follow the licensing terms described here, provided that
 * the new terms are clearly indicated on the first page of each file where
 * they apply.

 * IN NO EVENT SHALL THE AUTHORS OR DISTRIBUTORS BE LIABLE TO ANY PARTY
 * FOR DIRECT, INDIRECT, SPECIAL, INCIDENTAL, OR CONSEQUENTIAL DAMAGES
 * ARISING OUT OF THE USE OF THIS SOFTWARE, ITS DOCUMENTATION, OR ANY
 * DERIVATIVES THEREOF, EVEN IF THE AUTHORS HAVE BEEN ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 * 
 * THE AUTHORS AND DISTRIBUTORS SPECIFICALLY DISCLAIM ANY WARRANTIES,
 * INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE, AND NON-INFRINGEMENT.  THIS SOFTWARE
 * IS PROVIDED ON AN "AS IS" BASIS, AND THE AUTHORS AND DISTRIBUTORS HAVE
 * NO OBLIGATION TO PROVIDE MAINTENANCE, SUPPORT, UPDATES, ENHANCEMENTS, OR
 * MODIFICATIONS.
 * 
*/

using System;
using System.Collections;
using System.Collections.Generic;
using UnityEngine;

namespace PitchDetector {

    public class Downsample {

        public const float filterLength = 0.005f;

        int FFIR_fsize = 0;
        float[] FFIR_co = null;
        float[] FFIR_mem = null;
        float[] FFIR_state = new float[1000];
        public float[] b = new float[2048];

        public int ncoeff = 127;
        public int ncoefft = 0;

        public float[] GetDownsampledData(float[] fdata, int buffsize, int nframes, int decimate, float toFreq, Params par, bool first, bool last_time) {
            float[] dsdata;
            int samsds;
            samsds = ((nframes - 1) * par.step + par.ncomp) / decimate;
            if (samsds < 1)
                return null;

            if (!Apply(fdata, buffsize, out dsdata, par.sdstep, par.freq, ref samsds, decimate, first, last_time)) {
                Debug.Log("downsample failed");
                return null;
            }
            return dsdata;
        }


        private bool Apply(float[] input, int buff_size, out float[] foutput, int state_idx, double freq, ref int samsds,
            int decimate, bool first, bool last_time) {

            float beta;
            int init;
            foutput = null;
            if (input != null && buff_size > 0 && decimate > 0) {
                if (decimate == 1) {
                    foutput = input;
                    return true;
                }

                foutput = ArrayPool.Take(buff_size / decimate);
                if (first) {
                    ncoeff = ((int)(freq * filterLength)) | 1;
                    beta = .5f / decimate;

                    lc_lin_fir(beta, ref ncoeff, b);

                    ncoefft = (ncoeff / 2) + 1;
                }

                if (first) {
                    init = 1;
                } else if (last_time) {
                    init = 2;
                } else {
                    init = 0;
                }

                doFFIR(input, buff_size, foutput, ref samsds, state_idx, ncoefft, b, decimate, init);
                return true;
            }

            return false;
        }

        /**
        * create the coefficients for a symmetric FIR lowpass filter using the
        * window technique with a Hanning window.
        */
        private void lc_lin_fir(float fc, ref int nf, float[] coef) {

            if ((nf % 2) != 1) {
                nf = nf + 1;
            }

            int n = (nf + 1) / 2;

            /*  Compute part of the ideal impulse response (the sin(x)/x kernel). */
            const float twopi = Mathf.PI * 2.0f;
            coef[0] = (float)(2.0 * fc);
            float fn = twopi * fc;
            for (int i = 1; i < n; i++) {
                coef[i] = (Mathf.Sin(i * fn) / (Mathf.PI * i));
            }

            /* Now apply a Hanning window to the (infinite) impulse response. */
            fn = twopi / (float)(nf);
            for (int i = 0; i < n; i++) {
                coef[n - i - 1] *= (0.5f - (.5f * Mathf.Cos(fn * ((float)i + 0.5f))));
            }
        }

        /**
        * From jkGetF0.c
        * <p/>
        * fc contains 1/2 the coefficients of a symmetric FIR filter with unity
        * passband gain.
        * This filter is convolved with the signal in input.
        * The output is placed in output.
        * <p/>
        * If(invert), the filter magnitude
        * response will be inverted.
        * <p/>
        * If(init&1), beginning of signal is in input;
        * if(init&2), end of signal is in input.
        * outsize is set to the number of
        *
        * @param input     input data
        * @param insize    size of the input data
        * @param output    output data
        * @param outsize   size of the output data
        * @param state_idx index into the state array (storing previous input data)
        * @param ncoef     number of fir coefficients
        * @param fc        filter ceofficients
        * @param skip      number of samples to skip when downsampling
        * @param init      Is the beginning of the signal is already loaded? Is the end?
        */
        private void doFFIR(float[] input, int insize, float[] output, ref int outsize, int state_idx, int ncoef, float[] fc,
            int skip, int init) {
            // Reallocate static FIR filter parameters.
            if (ncoef > FFIR_fsize) {
                int ful_ffir_size = (ncoef + 1) * 2;

                FFIR_co = new float[ful_ffir_size];
                FFIR_mem = new float[ful_ffir_size];
                FFIR_fsize = ncoef;
            }

            // Fill the second half of mem with input data.
            int dp1_memidx = ncoef - 1;
            System.Array.Copy(input, 0, FFIR_mem, dp1_memidx, ncoef);
            int in_idx = ncoef;

            float sum;

            if ((init & 1) != 0) { /* Is the beginning of the signal in buf? */
                                   /* Copy the half-filter and its mirror image into the coefficient array. */
                int dp3_fcidx = ncoef - 1;
                int dp2_coidx = 0;
                int dp1_coidx = (ncoef - 1) * 2;

                for (int i = 0; i < ncoef - 1; i++) {
                    FFIR_co[dp1_coidx--] = FFIR_co[dp2_coidx++] = fc[dp3_fcidx--];
                }
                // set point of symmetry
                FFIR_co[dp1_coidx] = fc[dp3_fcidx];

                Array.Clear(FFIR_mem, 0, ncoef - 1);
            } else {
                System.Array.Copy(FFIR_state, 0, FFIR_mem, 0, ncoef - 1);
            }

            int resid;
            int k = (ncoef << 1) - 1;

            int out_idx = 0;

            for (int l = 0; l < outsize; l++) {
                dp1_memidx = 0;
                int dp2_coidx = 0;
                int dp3_memidx = skip;
                sum = 0.0f;
                for (int j = k - skip; j-- > 0;) {

                    sum += FFIR_co[dp2_coidx++] * FFIR_mem[dp1_memidx];
                    FFIR_mem[dp1_memidx++] = FFIR_mem[dp3_memidx++];
                }

                for (int j = skip; j-- > 0;) {
                    sum += FFIR_co[dp2_coidx++] * FFIR_mem[dp1_memidx];
                    FFIR_mem[dp1_memidx++] = input[in_idx++];
                }
                output[out_idx++] = (sum < 0.0) ? sum - 0.5f : sum + 0.5f;
            }

            if ((init & 2) != 0) {  // at the end
                resid = insize - outsize * skip;

                for (int l = resid / skip; l-- > 0;) {
                    dp1_memidx = 0;
                    int dp2_coidx = 0;
                    int dp3_memidx = skip;
                    sum = 0.0f;
                    for (int j = k - skip; j-- > 0;) {
                        sum += FFIR_co[dp2_coidx++] * FFIR_mem[dp1_memidx];
                        FFIR_mem[dp1_memidx++] = FFIR_mem[dp3_memidx++];
                    }
                    for (int j = skip; j-- > 0; FFIR_mem[dp1_memidx++] = 0.0f)
                        sum += FFIR_co[dp2_coidx++] * FFIR_mem[dp1_memidx];
                    output[out_idx++] = (sum < 0.0) ? sum - 0.5f : sum + 0.5f;
                    outsize++;
                }
            } else {  // not end end
                Array.Copy(input, state_idx - ncoef + 1, FFIR_state, 0, ncoef - 1);
            }
        }
    }
}