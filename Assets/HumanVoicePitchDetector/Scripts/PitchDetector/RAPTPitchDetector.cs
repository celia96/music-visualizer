/* RAPTPitchDetector.cs
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


using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using System;
using System.Linq;

#if UNSAFE_OPTIMISED
#warning UNSAFE_OPTIMISED defined: remember to add mcs.rsp file with text '-unsafe' to Asset folder
#endif
namespace PitchDetector {
    public class RAPTPitchDetector {
        public const int InputScale = (1 << (16 - 1));

        private Params par;
        private Downsample downsampler;

        public class Statistics {
            public int findBestPath;
            public int bestPathFound;
            public int framesWithFreqFound;
            public int totalFramesAnalyzed;
            public int voicedTransition;
            public int unvoicedTransition;
            public void Reset(){
                findBestPath = 0;
                bestPathFound = 0;
                framesWithFreqFound = 0;
                totalFramesAnalyzed = 0;
                voicedTransition = 0;
                unvoicedTransition = 0;
            }

            override public string ToString() {
                float pathFound = findBestPath > 0 ? (float)bestPathFound / findBestPath : 0;
                float freqFound = totalFramesAnalyzed > 0 ? (float)framesWithFreqFound / totalFramesAnalyzed : 0;
                int totalTransitions = voicedTransition + unvoicedTransition;
                float voicedTransitions = totalTransitions > 0 ? (float)voicedTransition / totalTransitions : 0;
                return "path found:" + pathFound.ToString() +
                                                " freq found:" + freqFound +
                                                " voiced transitions:" + voicedTransitions;
            }
        }

        public static Statistics statistics = new Statistics();

        public RAPTPitchDetector() {
            downsampler = new Downsample();
            par = new Params();
        }

        public RAPTPitchDetector(float minF0, float maxF0) {
            downsampler = new Downsample();
            par = new Params(minF0, maxF0);
        }

        public RAPTPitchDetector(float freq, float minF0, float maxF0) {
            downsampler = new Downsample();
            par = new Params(freq, minF0, maxF0);
        }


        public class WindowStat {
            public float err;
            public float rms;
            public float[] rho;

            public WindowStat(float err, float rms, float[] rho) {
                this.err = err;
                this.rms = rms;
                this.rho = rho;
            }

            public void Get(out float err, out float rms, out float[] rho) {
                err = this.err;
                rms = this.rms;
                rho = this.rho;
            }
        }

        private Stack<WindowStat> windStats = new Stack<WindowStat>();

        public struct PeakData {
            public static double ln2 = Math.Log(2);
            public float value;
            public int location;
            public void Set(float val, int loc) {
                value = val;
                location = loc;
            }

            public void Set(PeakData other) {
                value = other.value;
                location = other.location;
            }

            public float TransitionCost(ref PeakData prev, ref Stat stationarity, Params parameters) {
                if (location > 0 && prev.location > 0) { // current and prev voiced
                    statistics.voicedTransition += 1;
                    return VoicedTransition(ref prev, parameters);
                } else if (location < 0 && prev.location < 0) { // current and prev unvoiced
                    statistics.unvoicedTransition += 1;
                    return 0f;
                } else {
                    statistics.unvoicedTransition += 1;
                    return UnvoicedTransition(ref stationarity, parameters); // other unvoiced
                }
            }

            private float UnvoicedTransition(ref Stat stationarity, Params parameters) {
                return parameters.trans_cost
                    + (parameters.trans_spec * stationarity.stat)
                    + (parameters.trans_amp * stationarity.rms_ratio);
            }

            private float VoicedTransition(ref PeakData prev, Params parameters) {
                double ftemp = Math.Log((double)location / prev.location);
                double ft1;
                if (Math.Abs(ftemp) < Double.Epsilon) {
                    return 0f;
                } else if (ftemp > 0) {
                    ft1 = parameters.double_cost + Math.Abs(ftemp - ln2);
                    return (float)((ftemp > ft1) ? ft1 * parameters.freqwt : ftemp * parameters.freqwt);
                } else {
                    ft1 = parameters.double_cost + Math.Abs(ftemp + ln2);
                    return (float)((-ftemp > ft1) ? ft1 * parameters.freqwt : -ftemp * parameters.freqwt);
                }
            }

        }

        int pool_nlags;
        int pool_ncands;
        Stack<DPFrame> freeDPFrames = new Stack<DPFrame>();
        Stack<DPFrame> usedDPFrames = new Stack<DPFrame>();

        public DPFrame CreateDPFrame(int nlags, int ncands) {
            DPFrame r;
            if (nlags == pool_nlags && ncands == pool_ncands) {
                if (freeDPFrames.Count > 0) {
                    r = freeDPFrames.Pop();
                    r.Clear();
                } else {
                    r = new DPFrame(nlags, ncands);
                }
            } else {
                pool_nlags = nlags;
                pool_ncands = ncands;
                freeDPFrames.Clear();
                r = new DPFrame(nlags, ncands);
            }
            usedDPFrames.Push(r);
            return r;
        }

        public void FreeUsedDPFrames() {
            while (usedDPFrames.Count > 0) {
                freeDPFrames.Push(usedDPFrames.Pop());
            }
        }

        public class DPFrame {
            public DPRecord dp;
            public Cross cp;
            public float rms;

            public DPFrame next;
            public DPFrame prev;

            public DPFrame(int nlags, int ncands) {
                this.dp = new DPRecord(ncands);
                this.cp = new Cross(nlags);
                this.next = null;
                this.prev = null;
            }

            public void Clear() {
                this.rms = 0f;
                this.dp.Clear();
                this.cp.Clear();
                this.next = null;
                this.prev = null;
            }

        }



        public class DPRecord {
            public int ncands;
            public PeakData[] peaks;
            public float[] mpvals;
            public int[] prept;
            public float[] dpvals;

            public DPRecord(int ncands) {
                this.ncands = 0;
                peaks = new PeakData[ncands];
                prept = new int[ncands];
                mpvals = new float[ncands];
                dpvals = new float[ncands];
            }

            public void Clear() {
                Array.Clear(prept, 0, prept.Length);
                Array.Clear(peaks, 0, peaks.Length);
                Array.Clear(mpvals, 0, mpvals.Length);
                Array.Clear(dpvals, 0, dpvals.Length);
                ncands = 0;
            }

            public void CalculateMinimumCost(DPRecord prev, ref Stat stationarity, Params par, bool isFirst) {
                for (int k = 0; k < ncands; k++) {

                    if (isFirst) {
                        dpvals[k] = mpvals[k];
                        prept[k] = 0;
                        continue;
                    }

                    int minloc = 0;
                    float errmin = float.MaxValue;
                    for (int j = 0; j < prev.ncands; j++) {
                        float ferr = peaks[k].TransitionCost(ref prev.peaks[j], ref stationarity, par);
                        float err = ferr + prev.dpvals[j];
                        if (err < errmin) {
                            errmin = err;
                            minloc = j;
                        }
                    }
                    dpvals[k] = errmin + mpvals[k];
                    prept[k] = minloc;
                }
            }

            public void GetBestCandidate(out float minErr, out int idx) {
                dpvals.MinAndIndex(out idx, out minErr, ncands);
            }

            public int GetBestCandidate() {
                return dpvals.MinIndex(ncands);
            }

            public void ApplyWeighting(int count, Params par) {
                for (int i = 0; i < count; i++) {
                    float tmp = 1.0f - ((float)peaks[i].location * par.lagwt);
                    mpvals[i] = 1.0f - (peaks[i].value * tmp);
                }

            }
        }

        public class Cross {
            public float[] correl;
            public int maxloc;
            public float maxval;
            public float rms;
            public int firstlag;

            public Cross(int nlags) {
                correl = new float[nlags];
            }

            public void Clear() {
                Array.Clear(correl, 0, correl.Length);
            }
        }


        public struct Stat {
            public float stat;
            public float rms;
            public float rms_ratio;

            public void Set(float stat, float rms, float rms_ratio) {
                this.stat = stat;
                this.rms = rms;
                this.rms_ratio = rms_ratio;
            }
        }

        /*
        * DP_CIRCULAR: determines the initial size of DP circular buffer in sec
        * DP_HIST: stored frame history in second before checking for common path
        *      DP_CIRCULAR > READ_SIZE, DP_CIRCULAR at least 2 times of DP_HIST
        * DP_LIMIT: in case no convergence is found, DP frames of DP_LIMIT secs
        *      are kept before output is forced by simply picking the lowest cost
        *      path
        */
        private static float DP_CIRCULAR = 0.7f;
        private static float DP_HIST = 0.15f;
        private static float DP_LIMIT = 0.3f; // 1.0

        // Member variables for the RAPT algorithm
        protected DPFrame headF = null;   // Current frame in the circular buffer
        protected DPFrame tailF = null;   // Frame where tracks start
        protected DPFrame cmpthF = null;  // the starting frame of converged paths to backtrack

        int[] pcands = null;

        int size_frame_hist;
        int size_frame_out;
        int num_active_frames;


        PeakData[] peaks = null;

        bool first = true;

        public const int BIGSORD = 100;


        void initDP() {

            float frame_int = (float)((float)par.step / par.freq);

            size_frame_hist = (int)(DP_HIST / frame_int);
            size_frame_out = (int)(DP_LIMIT / frame_int);

            if (par.READ_SIZE >= DP_CIRCULAR ||
                DP_CIRCULAR <= 2 * DP_HIST) {
                throw new ArgumentException("invalid buffer parameters");
            }

            if (first) {
                statistics.Reset();
                // Allocate a circular buffer of DPFrames
                int size_circ_buf = (int)(DP_CIRCULAR / frame_int);

                FreeUsedDPFrames();

                tailF = CreateDPFrame(par.nlags, par.n_cands);
                headF = tailF;
                for (int j = 0; j < size_circ_buf; j++) {
                    headF.next = CreateDPFrame(par.nlags, par.n_cands);
                    headF.next.prev = headF;
                    headF = headF.next;
                }
                headF.next = tailF;
                tailF.prev = headF;
                headF = tailF;

                pcands = new int[par.n_cands];

                int maxpeaks = (short)(2 + (par.nlags / 2));
                peaks = new PeakData[maxpeaks];
                num_active_frames = 0;
            }
        }


        public DPFrame TraceBackBestPath(ref bool checkpathDone, ref int bestCand) {
            DPFrame frm = headF;
            int numPaths = frm.dp.ncands;
            Array.Copy(frm.dp.prept, pcands, numPaths);
            /* Starting from the most recent frame, trace back each candidate's
            best path until reaching a common candidate at some past frame. */
            statistics.findBestPath += 1;
            while (true) {
                frm = frm.prev;
                checkpathDone = pcands.AllSame(numPaths);
                if (!checkpathDone) { // prepare for checking at previous frame
                    for (int k = 0; k < numPaths; k++) {
                        pcands[k] = frm.dp.prept[pcands[k]];
                    }
                } else { // all paths have converged
                    bestCand = pcands[0];
                    statistics.bestPathFound += 1;
                    break;
                }
                if (frm == tailF) { // used all available data
                    frm = (num_active_frames >= size_frame_out) ? headF : null;
                    checkpathDone = (num_active_frames >= size_frame_out);
                    break;
                }
            }
            return frm;
        }


        public bool dpF0(float[] fdata, int buffsize, List<float> f0p, bool last_time) {
            const float downsampleFreq = 2000f;
            int decimate = (int)(par.freq / downsampleFreq);
            int nframes = (buffsize < par.pad) ? 0 : (int)((buffsize - par.pad) / par.step);
            float[] dsdata = (decimate > 1) ? downsampler.GetDownsampledData(fdata, buffsize, nframes, decimate, downsampleFreq, par, first, last_time) : fdata;
            if (dsdata == null) {
                return false;
            }

            // Get a function of the stationarity of the speech signal
            Stat[] stats = getStationarity(fdata, 0, buffsize, nframes, par.step, first);

            /***********************************************************************/
            /* MAIN FUNDAMENTAL FREQUENCY ESTIMATION LOOP */
            /***********************************************************************/

            // advance to the next frame
            if (!first && nframes > 0) {
                headF = headF.next;
            }

            int ncand;
            for (int i = 0; i < nframes; i++) {
                headF.rms = stats[i].rms;

                getFastCands(fdata, dsdata, i, decimate, headF.cp, peaks, out ncand, par);

                for (int j = 0; j < ncand; ++j) {
                    headF.dp.peaks[j].Set(peaks[j]);
                }

                headF.dp.peaks[ncand].Set(headF.cp.maxval, -1);
                headF.dp.mpvals[ncand] = par.voice_bias + headF.cp.maxval;

                /* Apply a lag-dependent weight to the peaks to encourage the selection
                of the first major peak.  Translate the modified peak values into
                costs (high peak ==> low cost). */
                headF.dp.ApplyWeighting(ncand, par);

                ncand++;            /* include the unvoiced candidate */
                headF.dp.ncands = ncand;
                headF.dp.CalculateMinimumCost(headF.prev.dp, ref stats[i], par, (first && i == 0));


                if (i < nframes - 1) {
                    headF = headF.next;
                }
            }
            if (decimate > 1) {
                ArrayPool.PutBack(dsdata);
            }
            // done propagating dp structure for the current set of sampled data.

            num_active_frames += nframes;

            if (num_active_frames >= size_frame_hist || last_time) {

                int best_cand = 0;
                bool checkpath_done = true;
                if (last_time) {
                    cmpthF = headF;
                    best_cand = headF.dp.GetBestCandidate();
                } else {
                    cmpthF = TraceBackBestPath(ref checkpath_done, ref best_cand);
                }

                if (checkpath_done) {
                    int count = 0;
                    // Backtracking from cmpthF (best cand) to tailf
                    DPFrame frm = cmpthF;  // start where convergence was found
                    while (frm != tailF.prev) {
                        int loc1 = frm.dp.peaks[best_cand].location;
                        best_cand = frm.dp.prept[best_cand];
                        float ftemp = loc1;

                        if (loc1 > 0) { // was f0 actually estimated for this frame
                            if (loc1 > par.start && loc1 < par.stop) {
                                ftemp += peak(frm.cp.correl, loc1 - par.start);
                            }
                            f0p.Add(par.freq / ftemp);
                            statistics.framesWithFreqFound += 1;
                        } else {
                            f0p.Add(0);
                        }
                        frm = frm.prev;
                        ++count;
                    }
                    statistics.totalFramesAnalyzed += count;
                    tailF = cmpthF.next;
                    num_active_frames -= count;
                    f0p.PartialReverse(f0p.Count - count, count);
                }
            }

            if (first) {
                first = false;
            }
            return true;
        }

        public void getFastCands(float[] fdata, float[] fdsdata, int ind, int dec, Cross cp, PeakData[] peaks, out int ncand, Params par) {

            float lag_wt = par.lag_weight / par.nlags;
            int decnlags = 1 + (par.nlags / dec);
            int decstart = Mathf.Max(par.start / dec, 1);

            int decind = (ind * par.step) / dec;
            int decsize = 1 + (par.size / dec);

            crossf(fdsdata, decind, decsize, decstart, decnlags, cp);

            get_cand(cp, peaks, decnlags, out ncand, par.cand_thresh); /* return high peaks in xcorr */

            /* Interpolate to estimate peak locations and values at high sample rate. */
            for (int i = 0; i < ncand; i++) {
                int j = peaks[i].location - decstart - 1;
                float xp, yp;
                peak(cp.correl, j, out xp, out yp);
                peaks[i].location = (peaks[i].location * dec) + (int)(0.5 + (xp * dec)); /* refined lag */
                peaks[i].value = yp * (1.0f - (lag_wt * peaks[i].location)); /* refined amplitude */
            }

            if (ncand >= par.n_cands) { /* need to prune candidates? */
                peaks.SortByVal(ncand);
                ncand = par.n_cands - 1;  /* leave room for the unvoiced hypothesis */
            }
            crossfi(fdata, (ind * par.step), 7, cp, peaks, ncand);

            get_cand(cp, peaks, par.nlags, out ncand, par.cand_thresh); /* return high peaks in xcorr */
            if (ncand >= par.n_cands) { /* need to prune candidates again? */
                peaks.SortByVal(ncand);
                ncand = par.n_cands - 1;  /* leave room for the unvoiced hypothesis */
            }

        }

        float[] dbdata = null;

        /* Return a sequence based on the normalized crosscorrelation of the
        * signal in data.  This is similar to crossf(), but is designed to
        * compute only small patches of the correlation sequence.  The length of
        * each patch is determined by nlags; the number of patches by nlocs, and
        * the locations of the patches is specified by the array locs.  Regions
        * of the CCF that are not computed are set to 0.
        * <p/>
        * data is the input speech array
        * doff is offset to input array
        * size is the number of samples in each correlation
        * start0 is the first (virtual) lag to compute (governed by highest F0)
        * nlags0 is the number of lags (virtual+actual) in the correlation sequence
        * nlags is the number of cross correlations to compute at each location
        * cp contains correlation data
        * locs is an array of indices pointing to the center of a patches where the
        * cross correlation is to be computed.
        * nlocs is the number of correlation patches to compute.
        */
        void crossfi(float[] data, int doff, int nlags, Cross cp, PeakData[] locs, int nlocs) {

            /* Compute mean in reference window and subtract this from the entire sequence. */
            int size = par.size; // cache size property
            int total = size + par.start + par.nlags;
            if (dbdata == null) {
                dbdata = ArrayPool.Take(total);
            } else if (total > dbdata.Length) {
                ArrayPool.PutBack(dbdata);
                dbdata = ArrayPool.Take(total);
            }

#if UNSAFE_OPTIMISED
            unsafe {
                fixed (float* ptrDbdata = dbdata, ptrData = data) {
                    float engr = data.Mean(size, doff);

                    float* a = ptrDbdata;
                    float* b = ptrData + doff;
                    for (int j = total; j > 0; --j) {
                        *a++ = *b++ - engr;
                    }

                    Array.Clear(cp.correl, 0, par.nlags);
                    /* compute energy in reference window */
                    engr = dbdata.SumOfSquared(size);

                    cp.rms = Mathf.Sqrt(engr / size);
                    cp.firstlag = par.start;

                    if (engr > 0.0) {
                        cp.maxloc = -1;
                        cp.maxval = 0f;
                        for (int li = 0; li < nlocs; li++) {
                            int start = Mathf.Max(locs[li].location - (nlags / 2), par.start);
                            /* compute energy at first requested lag */
                            float engc = dbdata.SumOfSquared(size, start);

                            /* COMPUTE CORRELATIONS AT ALL REQUESTED LAGS */
                            int i, ci;
                            for (i = 0, ci = start - par.start; i < nlags; i++, ci++) {
                                float sum = 0f;
                                a = ptrDbdata;
                                b = ptrDbdata + i + start;
                                for (int j = size; j > 0; --j)
                                    sum += (*a++) * (*b++);

                                engc = Mathf.Max(engc, 1f);
                                cp.correl[ci] = sum / Mathf.Sqrt(10000.0f + (engc * engr));
                                engc -= (dbdata[i + start] * dbdata[i + start]);
                                engc += (*b) * (*b);
                                if (cp.correl[ci] > cp.maxval) {
                                    cp.maxval = cp.correl[ci];
                                    cp.maxloc = i + start;
                                }
                            }
                        }
                    } else {
                        cp.maxloc = 0;
                        cp.maxval = 0.0f;
                    }
                }
            }
#else
		float engr = data.Mean (size, doff);

		for (int j = total - 1; j >= 0; --j) {
			dbdata[j] = data[j + doff] - engr;
		}

		Array.Clear (cp.correl, 0, par.nlags);
		/* compute energy in reference window */
		engr = dbdata.SumOfSquared (size);

		cp.rms = Mathf.Sqrt (engr / size);
		cp.firstlag = par.start;

		if (engr > 0.0) {
			cp.maxloc = -1;
			cp.maxval = 0f;
			for (int li = 0; li < nlocs; li++) {
				int start = Mathf.Max (locs [li].location - (nlags / 2), par.start);
				/* compute energy at first requested lag */
				float engc = dbdata.SumOfSquared (size, start);

				/* COMPUTE CORRELATIONS AT ALL REQUESTED LAGS */
				int i, ci;
				for (i = 0, ci = start - par.start; i < nlags; i++, ci++) {
					float sum = 0f;
					int ds, dbi;
					for (dbi = 0, ds = i + start; dbi < size; ++dbi, ++ds)
						sum += dbdata [dbi] * dbdata [ds];

					engc = Mathf.Max (engc, 1f);
					cp.correl [ci] = sum / Mathf.Sqrt (10000.0f + (engc * engr));
					engc -= (dbdata [i + start] * dbdata [i + start]);
					engc += dbdata[ds] * dbdata[ds];
					if (cp.correl [ci] > cp.maxval) {
						cp.maxval = cp.correl [ci];
						cp.maxloc = i + start;
					}
				}
			}
		} else {
			cp.maxloc = 0;
			cp.maxval = 0.0f;
		}
#endif
        }


        /**
        * +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        * Return a sequence based on the normalized crosscorrelation of the signal
        * in data.
        * <p/>
        * data is the input speech array
        * doff is offset of data
        * size is the number of samples in each correlation
        * start is the first lag to compute (governed by the highest expected F0)
        * nlags is the number of cross correlations to compute (set by lowest F0)
        * cp contains correlation data
        */
        public void crossf(float[] data, int doff, int size, int start, int nlags, Cross cp) {
            /* Compute mean in reference window and subtract this from the
            entire sequence.  This doesn't do too much damage to the data
            sequenced for the purposes of F0 estimation and removes the need for
            more principled (and costly) low-cut filtering. */
            int total = size + start + nlags;
            if (dbdata == null) {
                dbdata = ArrayPool.Take(total);
            } else if (total > dbdata.Length) {
                ArrayPool.PutBack(dbdata);
                dbdata = ArrayPool.Take(total);
            }
#if UNSAFE_OPTIMISED
            unsafe {
                fixed (float* ptrDbdata = dbdata, ptrData = data) {
                    float engr = data.Mean(size, doff);
                    float* a = ptrDbdata;
                    float* b = ptrData + doff;
                    for (int j = total; j > 0; --j) {
                        *a++ = (*b++) - engr;
                    }
                    engr = dbdata.SumOfSquared(size);

                    cp.firstlag = start;
                    cp.rms = Mathf.Sqrt(engr / size);
                    if (engr > 0.0) {    /* If there is any signal energy to work with... */
                                         /* Compute energy at the first requested lag. */
                        float engc = dbdata.SumOfSquared(size, start);

                        cp.maxval = 0f;
                        cp.maxloc = -1;
                        /* COMPUTE CORRELATIONS AT ALL OTHER REQUESTED LAGS. */
                        for (int i = 0; i < nlags; i++) {
                            float sum = 0.0f;
                            //int ds, dbi;
                            //for (dbi = 0, ds = i + start; dbi < size; ++dbi,++ds)
                            //	sum += dbdata [dbi] * dbdata [ds];
                            a = ptrDbdata;
                            b = ptrDbdata + i + start;
                            for (int j = size; j > 0; --j)
                                sum += (*a++) * (*b++);
                            cp.correl[i] = (sum / Mathf.Sqrt(engc * engr)); /* output norm. CC */
                            engc -= dbdata[i + start] * dbdata[i + start]; /* adjust norm. energy for next lag */
                            engc += (*b) * (*b);
                            if (engc < 1f) {
                                engc = 1f;
                            }
                            if (cp.correl[i] > cp.maxval) {     /* Find abs. max. as we go. */
                                cp.maxval = cp.correl[i];
                                cp.maxloc = i + start;
                            }
                        }
                    } else {    /* No energy in signal; fake reasonable return vals. */
                        cp.maxloc = 0;
                        cp.maxval = 0.0f;
                        Array.Clear(cp.correl, 0, nlags);
                    }
                }
            }
#else
		float engr = data.Mean (size, doff);

		for (int j = total - 1; j >= 0; --j) {
			dbdata[j] = data[j + doff] - engr;
		}
		engr = dbdata.SumOfSquared (size);

		cp.firstlag = start;
		cp.rms = Mathf.Sqrt (engr / size);
		if (engr > 0.0) {    /* If there is any signal energy to work with... */
			/* Compute energy at the first requested lag. */
			float engc = dbdata.SumOfSquared (size, start);

			cp.maxval = 0f;
			cp.maxloc = -1;
			/* COMPUTE CORRELATIONS AT ALL OTHER REQUESTED LAGS. */
			for (int i = 0; i < nlags; i++) {
				float sum = 0.0f;
				int ds, dbi;
				for (dbi = 0, ds = i + start; dbi < size; ++dbi,++ds)
					sum += dbdata [dbi] * dbdata [ds];

				cp.correl [i] = (sum / Mathf.Sqrt (engc * engr)); /* output norm. CC */
				engc -= dbdata [i + start] * dbdata [i + start]; /* adjust norm. energy for next lag */
				engc += dbdata[ds] * dbdata[ds];
				if (engc < 1f) {
					engc = 1f;
				}
				if (cp.correl [i] > cp.maxval) {		/* Find abs. max. as we go. */
					cp.maxval = cp.correl [i];
					cp.maxloc = i + start;
				}
			}
		} else {	/* No energy in signal; fake reasonable return vals. */
			cp.maxloc = 0;
			cp.maxval = 0.0f;
			Array.Clear (cp.correl, 0, nlags);
		}
#endif

        }

        /**
       * -----------------------------------------------------------------------
       * * Use parabolic interpolation over the three points defining the peak
       * vicinity to estimate the "true" peak.
       */
        public static void peak(float[] y, int idx, out float xp, out float yp) {
            float a = (float)((y[2 + idx] - y[1 + idx]) + (.5f * (y[idx] - y[2 + idx])));
            if (Mathf.Abs(a) > .000001f) {
                float c = (y[idx] - y[2 + idx]) / (4.0f * a);
                xp = c;
                yp = y[1 + idx] - (a * c * c);
            } else {
                xp = 0.0f;
                yp = y[1 + idx];
            }
        }

        public static float peak(float[] arr, int idx) {
            float cormax = arr[idx];
            float cprev = arr[idx + 1];
            float cnext = arr[idx - 1];

            float den = (float)(2.0 * (cprev + cnext - (2 * cormax)));
            // parabolic interpolation
            if (Mathf.Abs(den) > 0.000001f) {
                return 2.0f - ((((5.0f * cprev) + (3.0f * cnext) - (8.0f * cormax)) / den));
            } else {
                return 0f;
            }
        }


        /* Get likely candidates for F0 peaks. */
        public static void get_cand(Cross cross, PeakData[] peak, int nlags, out int ncand, float cand_thresh) {
            float clip = (cand_thresh * cross.maxval);
            ncand = 0;
            int peakIdx = 0;
            for (int i = 1; i < nlags - 2; i++) {
                if (cross.correl[i] > clip && cross.correl[i] >= cross.correl[i - 1] && cross.correl[i] >= cross.correl[i + 1]) {
                    //if (cross.correl.IsPeak(i, clip)) { 
                    peak[peakIdx].Set(cross.correl[i], i + cross.firstlag);
                    peakIdx++;
                    ncand++;            /* count number of peaks found */
                }
            }
        }


        public Stat[] statArray = null;
        public float[] mem = null;

        public Stat[] getStationarity(float[] fdata, int fidx, int buff_size, int nframes, int frame_step,
            bool first) {

            const float preemp = 0.4f;
            const float stab = 30.0f;

            // initialize stat and memory space
            if (statArray == null || statArray.Length < nframes || first) {
                statArray = new Stat[nframes];
                mem = new float[par.stat_wsize + par.agap];
            }

            if (nframes == 0) {
                return statArray;
            }

            int q = fidx + par.ind;
            int datend = fidx + buff_size;
            int order = Mathf.Min(BIGSORD, (int)(2.0f + (par.freq / 1000.0f)));

            Array.Copy(fdata, 0, mem, mem.Length / 2, mem.Length - mem.Length / 2);


            for (int j = 0, p = q - par.agap; j < nframes; j++, p += frame_step, q += frame_step) {
                if ((p >= fidx) && (q >= fidx) && (q + par.stat_wsize <= datend)) {
                    statArray[j] = getSimilarity(order, par.stat_wsize, fdata, p, q, preemp, stab, false);
                } else {
                    if (first) {
                        if ((p < fidx) && (q >= fidx) && (q + par.stat_wsize <= datend)) {
                            statArray[j] = getSimilarity(order, par.stat_wsize, fdata, -1, q, preemp, stab, true);
                        } else {
                            statArray[j].Set(0.1f * 0.2f, 0f, 1.0f);
                        }
                    } else {
                        if ((p < fidx) && (q + par.stat_wsize <= datend)) {
                            statArray[j] = getSimilarity(order, par.stat_wsize, mem, 0, (mem.Length / 2) + par.ind, preemp, stab, false);
                            // Prepare the next frame.
                            if (p + frame_step < fidx) {
                                // slide memory to the left
                                Array.Copy(mem, frame_step, mem, 0, mem.Length - frame_step);
                                Array.Copy(fdata, q + par.stat_wsize, mem, mem.Length - frame_step, frame_step);
                            }
                        }
                    }
                }
            }

            int copyLength = Mathf.Min(mem.Length / 2, (nframes * frame_step));
            Array.Copy(fdata, fidx, mem, 0, copyLength);
            return statArray;
        }


        public Stat getSimilarity(int order, int size, float[] data, int dprev, int dcur, float preemp,
            float stab, bool init) {

            float b0 = 0f;
            float t;
            Stat r = new Stat();
            float rms3 = windEnergy(data, dcur, size);
            float[] a2 = ArrayPool.Take(BIGSORD + 1);
            float[] rho3 = ArrayPool.Take(BIGSORD + 1);
            float err3;
            xlpc(order, stab, size - 1, data, dcur, a2, rho3, out err3, preemp);
            if (!init) {
                float[] rho1;
                float err1, rms1;
                if (windStats.Count > 0) {
                    var s = windStats.Pop();
                    s.Get(out err1, out rms1, out rho1);
                } else {
                    rho1 = ArrayPool.Take(BIGSORD + 1);
                    xlpc(order, stab, size - 1, data, dprev, null, rho1, out err1, preemp);
                    rms1 = windEnergy(data, dprev, size);
                }
                float[] b = ArrayPool.Take(BIGSORD + 1);
                xaToaca(a2, b, out b0, order);
                t = xitakura(order, b, b0, rho1, 1, err1) - 0.8f;
                ArrayPool.PutBack(rho1);
                ArrayPool.PutBack(b);
                if (rms1 > 0.0f) {
                    r.rms_ratio = (0.001f + rms3) / rms1;
                } else if (rms3 > 0.0f) {
                    r.rms_ratio = 2f;
                } else {
                    r.rms_ratio = 1f;
                }
            } else {
                t = 10.0f;
                r.rms_ratio = 1.0f;
            }
            ArrayPool.PutBack(a2);
            var savedStats = new WindowStat(err3, rms3, rho3);
            windStats.Push(savedStats);
            r.rms = rms3;
            r.stat = 0.2f / t;
            return r;
        }

        /**
        * Compute the autocorrelations of the p LP coefficients in a.
        * (a[0] is assumed to be = 1 and not explicitely accessed.)
        * The magnitude of a is returned in c.
        * 2* the other autocorrelation coefficients are returned in b.
        */
        void xaToaca(float[] a, float[] b, out float c, int p) {
            float s;
            int ap_idx = 0;
            int a0_idx;
            int b_idx = 0;

            int i, j;

            for (s = 1.0f, i = p; i-- > 0; ap_idx++)
                s += a[ap_idx] * a[ap_idx];

            c = s;
            for (i = 1; i <= p; i++) {
                s = a[i - 1];
                for (a0_idx = 0, ap_idx = i, j = p - i; j-- > 0;)
                    s += (a[a0_idx++] * a[ap_idx++]);
                b[b_idx++] = (float)(2.0 * s);
            }
        }

        /**
        * Compute the Itakura LPC distance between the model represented
        * by the signal autocorrelation (r) and its residual (gain) and
        * the model represented by an LPC autocorrelation (c, b).
        * Both models are of order p.
        * r is assumed normalized and r[0]=1 is not explicitly accessed.
        * Values returned by the function are >= 1.
        */
        float xitakura(int p, float[] b, float c, float[] r, int r_offset, float gain) {
            float s = c;
            for (int i = 0; i < p; ++i)
                s += r[r_offset + i] * b[i];

            return (s / gain);
        }

        /**
        * Compute the time-weighted RMS of a size segment of data.  The data
        * is weighted by a window of type w_type before RMS computation.  w_type
        * is decoded above in window().
        *
        * @param data   input data
        * @param didx   offset into the input data
        * @param size   size of the window overwhich to calculate energy
        * @param w_type window type
        * @return the energy in the window
        */
        public float windEnergy(float[] data, int didx, int size) {

            float[] window = Window.Hanning(size);
            float sum = 0.0f;
            for (int i = 0; i < size; i++) {
                float f = window[i] * data[i + didx];
                sum += f * f;
            }

            return Mathf.Sqrt(sum / size);
        }


        public void xlpc(int lpc_ord, float lpc_stabl, int wsize, float[] data, int didx, float[] lpca, float[] ar,
            out float normerr, float preemp) {

            if ((wsize <= 0) || (data == null) || (lpc_ord > BIGSORD))
                throw new ArgumentException("xlpc invalid parameters");

            float[] dwind = ArrayPool.Take(wsize);

            Window.ApplyHanning(data, didx, dwind, wsize, preemp);
            float[] rho = (ar != null) ? ar : ArrayPool.Take(BIGSORD + 1);
            float[] a = (lpca != null) ? lpca : ArrayPool.Take(BIGSORD + 1);

            xautoc(wsize, dwind, lpc_ord, rho);
            ArrayPool.PutBack(dwind);

            if (lpc_stabl > 1) {  // add a little to the diagonal
                float ffact = (float)(1 / (1 + Math.Exp((-lpc_stabl / 20) * Math.Log(10))));
                for (int i = 1; i <= lpc_ord; i++) {
                        rho[i] *= ffact;
                }
            }
            xdurbin(rho, a, lpc_ord, out normerr);

            a[0] = 1;
            if (ar == null) {
                ArrayPool.PutBack(rho);
            }
            if (lpca == null) {
                ArrayPool.PutBack(a);
            }
        }

        /**
        * Compute the pp+1 autocorrelation lags of the windowsize samples in s.
        * Return the normalized autocorrelation coefficients in r.
        * The rms is returned in e.
        */
        public void xautoc(int wsize, float[] s, int p, float[] r) {

            float sum, sum0 = 0;
#if UNSAFE_OPTIMISED
            unsafe {
                fixed (float* sptr = s) {
                    float* pa = sptr;
                    for (int j = wsize - 1; j >= 0; --j) {
                        sum = *pa++;
                        sum0 += sum * sum;
                    }
                    r[0] = 1f;
                    if (sum0 == 0.0f) {  // no energy.  fake low-energy white noise autocorr
                        Array.Clear(r, 1, p);
                        return;
                    }
                    sum0 = 1.0f / sum0;
                    for (int i = 1; i <= p; i++) {
                        sum = 0.0f;
                        pa = sptr;
                        float* pb = sptr + i;
                        for (int j = wsize - i; j >= 0; --j) {
                            sum += (*pa++) * (*pb++);
                        }
                    }
                }
            }
#else
    		for (int j = wsize - 1; j >= 0; --j) {
    			sum = s[j];
    			sum0 += sum * sum;
    		}
    		r [0] = 1f;
    		if (sum0 == 0.0f) {  // no energy.  fake low-energy white noise autocorr
    			Array.Clear (r, 1, p);
    			return;
    		}
    		sum0 = 1.0f / sum0;
    		for (int i = 1; i <= p; i++) {
    			sum = 0.0f;
                for (int j = 0; j < wsize - i; j++) {
                    sum += s[j] * s[j + i];
                }
                r[i] = sum * sum0;
    		}
#endif
        }

        /**
        * Using Durbin's recursion, convert the autocorrelation sequence in r
        * to reflection coefficients in k and predictor coefficients in a.
        * The prediction error energy (gain) is left in ex[0].
        * Note: durbin returns the coefficients in normal sign format.
        * (i.e. a[0] is assumed to be = +1.)
        */
        public void xdurbin(float[] r, float[] a, int p, out float ex) {

            float[] b = ArrayPool.Take(BIGSORD);
            float[] k = ArrayPool.Take(BIGSORD);
            float e = r[0];
            k[0] = -r[1] / e;
            a[0] = k[0];
            e *= (float)(1.0 - (k[0] * k[0]));

            for (int i = 1; i < p; i++) {
                float s = 0;
                for (int j = 0; j < i; j++) {
                    s -= a[j] * r[i - j];
                }
                k[i] = (s - r[i + 1]) / e;
                a[i] = k[i];
                Array.Copy(a, 0, b, 0, i + 1);
                for (int j = 0; j < i; j++) {
                    a[j] += k[i] * b[i - j - 1];
                }
                e *= (float)(1.0 - (k[i] * k[i]));
            }

            ex = e;
            ArrayPool.PutBack(b);
            ArrayPool.PutBack(k);
        }


        public int RequiredSampleCount(float sampleRate) {
            return (int)((par.frame_step * 2 + par.wind_dur) * sampleRate);
        }


        public List<float> getPitch(float[] wav, ref float db, float sampleRate) {
            int end = wav.Length;
            int start = 0;
            return getPitch(wav, start, ref end, ref db, sampleRate, false, false);
        }

        /**
        * Calculates the pitch of the wave file, wav, and returns a list of pitch values.
        * <p/>
        * Uses the spsk get_f0 algorithm.
        *
        * @param wav the wave data in -1..1 range
        * @param start start index of data in wav array
        * @param end end index of data in wav array
        * @param db of analyzed samples
        * @param sampleRate signal sample rate
        * @param streaming is streaming mode (i.e. more samples expected)
        * @param last is last sample of streaming
        * @return a pitch list
        */
        public List<float> getPitch(float[] wav, int start, ref int end, ref float db, float sampleRate, bool streaming, bool last) {

            par.freq = sampleRate;
            int total_samples;
            if (end == wav.Length && start == 0) {
                total_samples = wav.Length;
            } else {
                total_samples = (wav.Length + end - start) % wav.Length;
            }
            if (total_samples < RequiredSampleCount(sampleRate)) {
                var errorMsg = String.Format("Too short audio! sample count {0} < (({1} * 2) + {2}) * {3}",
                    total_samples, par.frame_step, par.wind_dur, sampleRate);
                throw new ArgumentException(errorMsg);
            }

            if (!streaming && !first) {
                Debug.LogWarning("Re-initialising while streaming!");
            }

            initDP();

            int buffsize = (par.buffsize > total_samples) ? total_samples : par.buffsize;
            int actsize = buffsize < total_samples ? buffsize : total_samples;

            bool done;
            float[] fdata = ArrayPool.Take(Mathf.Max(buffsize, par.sdstep));
            int length = total_samples;
            int ndone = 0;
            float sum_squared = 0f;

            List<float> f0p = new List<float>();
            while (true) {
                done = ((actsize < buffsize) || (total_samples == buffsize));

                if (!first && done && streaming && !last) {
                    break;
                }

                loadData(wav, start, fdata, ndone, actsize, ref sum_squared);

                if (!dpF0(fdata, actsize, f0p, done)) {
                    Debug.Log("dpF0 failed.");
                }

                if (done)
                    break;
                
                ndone += par.sdstep;
                actsize = Math.Min(buffsize, length - ndone);
                total_samples -= par.sdstep;

                if (actsize > total_samples) {
                    actsize = total_samples;
                }
            }
            end = (start + ndone) % wav.Length;
            if (ndone != 0) {
                sum_squared /= (start == 0 && end == wav.Length) ? wav.Length : ndone;
                float rms = Mathf.Sqrt(sum_squared);
                db = 20 * Mathf.Log10(rms);
            } else {
                db = -160f;
            }

            if (!streaming || streaming && last) {
                first = true;
                while (windStats.Count > 0) {
                    var p = windStats.Pop();
                    ArrayPool.PutBack(p.rho);
                }
            }
            ArrayPool.PutBack(fdata);
            return f0p;
        }

        /**
        * Copies and scales audio data from wav data into float based on an offset.
        *
        * @param wav      The wave data array.
        * @param start    offset in wav array.
        * @param fdata    a preallocated array to store the raw PCM data as a float.
        * @param pos      pcm frame offset
        * @param nsamples the number of samples to load
        */
        public void loadData(float[] wav, int start, float[] fdata, int pos, long nsamples, ref float sum) {
            int i = 0;
            int j = (start + pos) % wav.Length;
            float sample;
            while (i < nsamples && j < wav.Length) {
                sample = wav[j++];
                fdata[i++] = sample * (float)InputScale;
                sum += sample * sample;
            }
            j = 0;
            while (i < nsamples && j < wav.Length) {
                sample = wav[j++];
                fdata[i++] = sample * (float)InputScale;
                sum += sample * sample;
            }
            while (i < fdata.Length) {
                fdata[i++] = 0.0f;
            }
        }


    }
}