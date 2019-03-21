/* Window.cs
*
* Copyright 2017 Outloud Oy
* Written by: mikko.k.koivisto@outloud.fi
* 
* float array containing Hanning window
*/
using System.Collections;
using System.Collections.Generic;
using UnityEngine;

namespace PitchDetector {
    public class Window {

        private static Dictionary<int, float[]> hanningList =
            new Dictionary<int, float[]>();

        private static readonly object listAccessLock = new object();

        public static float[] Hanning(int n) {
            float[] r = null;

            if (hanningList.ContainsKey(n)) {
                r = hanningList[n];
            } else {
                r = CreateHanning(n);
                lock (listAccessLock) {
                    hanningList[n] = r;
                }
            }
            return r;
        }

        /**
        * Apply a Hanning window to the signal in din.
        * <p/>
        *
        * @param din    input data
        * @param didx   offset into the input data
        * @param dout   output data
        * @param n      the size of the window
        * @param preemp the preemphasis to apply (if any)
        */
        public static void ApplyHanning(float[] din, int didx, float[] dout, int n, float preemp) {
            // Generate the window if it doesn't exist
            float[] window = Window.Hanning(n);

            // Apply window and any preemphasis.
            if (preemp != 0) {
                for (int i = 0; i < n; i++) {
                    dout[i] = window[i] * (din[didx + i + 1] - (preemp * din[didx + i]));
                }
            } else {
                for (int i = 0; i < n; i++) {
                    dout[i] = window[i] * din[didx + i];
                }
            }
        }

        static float[] CreateHanning(int n) {
            float[] r = new float[n];
            float arg = (Mathf.PI * 2f) / (float)n;
            for (int i = 0; i < n; i++) {
                r[i] = 0.5f - 0.5f * Mathf.Cos((0.5f + i) * arg);
            }
            return r;
        }

    }
}