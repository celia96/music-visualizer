/* Params.cs
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

namespace PitchDetector {

    public class Params {
        /**
        * stationarity parameters -
        * STAT_WSIZE: window size in sec used in measuring frame energy/stationarity
        * STAT_AINT: analysis interval in sec in measuring frame energy/stationarity
        */
        public const float STAT_WSIZE = 0.030f;
        public const float STAT_AINT = 0.020f;

        public Params(float sampleRate, float minf0, float maxf0) {
            READ_SIZE = 0.2f;
            trans_cost = 0.005f;
            trans_amp = 0.5f;
            trans_spec = 0.5f;
            voice_bias = 0.0f;
            double_cost = 0.4;
            n_cands = 20;
            cand_thresh = 0.4f;

            _lagweight = 0.3f;
            _freqweight = 0.02;
            _winddur = 0.01125f; // original 0.0075f;
            _framestep = 0.015f; // original 0.01f

            _freq = sampleRate;
            _minf0 = minf0;
            _maxf0 = maxf0;

            RecalculateParams();
        }

        public Params(float minf0, float maxf0) : this(16000f, minf0, maxf0) {
        }

        public Params() : this(16000f, 50f, 800f) {
        }

        public int step { get { return _step; } }
        public int size { get { return _size; } }

        public int start { get { return _start; } }
        public int stop { get { return _stop; } }

        public int nlags { get { return _nlags; } }
        public int ncomp { get { return _ncomp; } }

        public double freqwt { get { return _freqwt; } }
        public float lagwt { get { return _lagwt; } }

        public int pad { get { return _pad; } }
        public int stat_wsize { get { return _statwsize; } }
        public int agap { get { return _agap; } }
        public int ind { get { return _ind; } }

        public int buffsize { get { return _buffsize; } }
        public int sdstep { get { return _sdstep; } }
        public int nframes { get { return _nframes; } }

        public float freq {
            get {
                return _freq;
            }
            set {
                _freq = value;
                RecalculateParams();
            }
        }

        public float maxF0 {
            get {
                return _maxf0;
            }
            set {
                _maxf0 = value;
                RecalculateParams();
            }
        }

        public float minF0 {
            get {
                return _minf0;
            }
            set {
                _minf0 = value;
                RecalculateParams();
            }
        }

        public float frame_step {
            get {
                return _framestep;
            }
            set {
                _framestep = value;
                RecalculateParams();
            }
        }

        public float wind_dur {
            get {
                return _winddur;
            }
            set {
                _winddur = value;
                RecalculateParams();
            }
        }

        public float lag_weight {
            get {
                return _lagweight;
            }
            set {
                _lagweight = value;
                RecalculateParams();
            }
        }

        public double freq_weight {
            get {
                return _freqweight;
            }
            set {
                _freqweight = value;
                RecalculateParams();
            }
        }

        // calculate parameters that have dependencies
        private void RecalculateParams() {
            _step = Mathf.RoundToInt(_framestep * _freq);
            _size = Mathf.RoundToInt(_winddur * _freq);
            _freqwt = freq_weight / ((double)_step / _freq);
            _start = Mathf.RoundToInt(_freq / _maxf0);
            _stop = Mathf.RoundToInt(_freq / _minf0);
            _nlags = _stop - _start + 1;
            _ncomp = _size + _stop + 1;
            _lagwt = _lagweight / _stop;

            _statwsize = (int)(STAT_WSIZE * _freq);
            _agap = (int)(STAT_AINT * _freq);
            _ind = (_agap - _statwsize) / 2;
            int i = _statwsize + _ind;
            int downpatch = (((int)(_freq * Downsample.filterLength)) + 1) / 2;
            _pad = downpatch + ((i > _ncomp) ? i : _ncomp);

            i = (int)(READ_SIZE * _freq);
            if (_ncomp >= _step) {
                _nframes = ((i - _ncomp) / _step) + 1;
            } else {
                _nframes = i / _step;
            }

            _buffsize = _nframes * _step + _pad;
            _sdstep = _nframes * _step;

            if (_buffsize <= 0 || _sdstep <= 0) {
                throw new System.ArgumentException("Invalid parameters");
            }
        }

        public float READ_SIZE;
        public float trans_cost;
        public float trans_amp;
        public float trans_spec;
        public float voice_bias;
        public double double_cost;

        public int n_cands;
        public float cand_thresh;

        private double _freqwt;
        private float _freq;
        private int _step;
        private int _size;
        private int _start;
        private int _stop;
        private int _nlags;
        private int _ncomp;
        private float _lagwt;
        private float _maxf0;
        private float _minf0;
        private float _framestep;
        private float _winddur;
        private float _lagweight;
        private double _freqweight;
        private int _pad;
        private int _agap;
        private int _ind;
        private int _statwsize;
        private int _nframes;
        private int _sdstep;
        private int _buffsize;
    }
}