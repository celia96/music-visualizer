/* ArrayPool.cs
*
* Copyright 2017 Outloud Oy
* Written by: mikko.k.koivisto@outloud.fi
* 
* float array pool
*/
using System;
using System.Collections;
using System.Collections.Generic;
using UnityEngine;

namespace PitchDetector {
    public class ArrayPool {

        private static Dictionary<int, Stack<float[]>> freeList =
            new Dictionary<int, Stack<float[]>>();

        private static readonly object listAccessLock = new object();

        public static float[] Take(int n) {
            float[] r = null;
            lock (listAccessLock) {
                if (freeList.ContainsKey(n)) {
                    Stack<float[]> l = freeList[n];
                    if (l.Count > 0) {
                        r = l.Pop();
                    }
                }
            }
            if (r == null) {
                r = new float[n];
            } else {
                Array.Clear(r, 0, r.Length);
            }
            return r;
        }

        public static void PutBack(float[] arr) {
            if (arr == null)
                return;
            lock (listAccessLock) {
                if (!freeList.ContainsKey(arr.Length))
                    freeList[arr.Length] = new Stack<float[]>();
                freeList[arr.Length].Push(arr);
            }

        }
    }
}