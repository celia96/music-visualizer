/* RAPTPitchDetectorExtensions.cs
*
* Copyright 2017 Outloud Oy
* Written by: mikko.k.koivisto@outloud.fi
* 
* various helper functions for RAPTPitchDetector
*/
using System;
using System.Collections;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using UnityEngine;

namespace PitchDetector {
    public static class RAPTPitchDetectorExtensions {
        public const int UndefinedMidi = int.MinValue;
        static string[] notes = new string[] { "C", "C#", "D", "D#", "E", "F", "F#", "G", "G#", "A", "A#", "B" };

        public static void PartialReverse(this List<float> data, int start, int length) {
            for (int i = 0; i < length / 2; ++i) {
                float tmp = data[start + i];
                data[start + i] = data[start + length - i - 1];
                data[start + length - i - 1] = tmp;
            }
        }

        public static float SumOfSquared(this float[] data, int length, int offset = 0) {
            float sum = 0f;
            if (length + offset > data.Length) {
                return sum;
            }
#if UNSAFE_OPTMISED
		unsafe {
			fixed (float* ptrData = data) {
				float* a = ptrData;
				for (int i = length; i > 0; --i) {
					var value = *a++;
					sum += value * value;
				}
			}
		}
#else
            for (int i = length - 1; i >= 0; --i) {
                var value = data[i];
                sum += value * value;
            }
#endif
            return sum;
        }

        public static float Mean(this float[] data, int length, int offset = 0) {
            float sum = 0f;
            for (int i = 0; i < length && i + offset < data.Length; ++i) {
                sum += data[i + offset];
            }
            return sum / length;
        }

        public static bool AllSame(this int[] data, int length) {
            for (int i = 1; i < data.Length && i < length; ++i) {
                if (data[0] != data[i]) {
                    return false;
                }
            }
            return true;
        }

        public static bool AllSame(this float[] data, int length) {
            for (int i = 1; i < data.Length && i < length; ++i) {
                if (data[0] != data[i]) {
                    return false;
                }
            }
            return true;
        }

        public static string PrintArray(this float[] data) {
            var strings = data.Select<float, string>(j => j.ToString()).ToArray();
            return string.Join(",", strings);
        }


        public static int MinIndex(this float[] data, int length = int.MaxValue) {
            length = Mathf.Min(length, data.Length);
            float minVal = float.MaxValue;
            int idx = 0;
            for (int i = 0; i < length; ++i) {
                if (data[i] < minVal) {
                    minVal = data[i];
                    idx = i;
                }
            }
            return idx;
        }

        public static void MinAndIndex(this float[] data, out int idx, out float minVal, int length = int.MaxValue) {
            length = Mathf.Min(length, data.Length);
            minVal = float.MaxValue;
            idx = -1;
            for (int i = 0; i < length; ++i) {
                if (data[i] < minVal) {
                    minVal = data[i];
                    idx = i;
                }
            }
        }

        public static bool IsPeak(this float[] data, int idx, float threshold) {
            return data[idx] > threshold && data[idx] >= data[idx - 1] && data[idx] >= data[idx + 1];
        }


        public static void SortByVal(this RAPTPitchDetector.PeakData[] peaks, int count) {
            for (int i = 1; i < count; i++) {
                var tmpPeak = peaks[i];
                int j = i;
                while (j > 0 && peaks[j - 1].value < tmpPeak.value) {
                    peaks[j] = peaks[j - 1];
                    j--;
                }
                peaks[j] = tmpPeak;
            }

        }

        public static void SortBoth(float[] peaks, int[] locations, int count) {
            for (int i = 1; i < count; i++) {
                float tmpPeak = peaks[i];
                int tmpLoc = locations[i];
                int j = i;
                while (j > 0 && peaks[j - 1] < tmpPeak) {
                    peaks[j] = peaks[j - 1];
                    locations[j] = locations[j - 1];
                    j--;
                }
                peaks[j] = tmpPeak;
                locations[j] = tmpLoc;
            }
        }

        public static void RemoveGlitches(this List<int> data, int threshold, int glitchWidth = 1) {
            int i = 0;
            while (i + glitchWidth + 1 < data.Count) {
                int prev = data[i];
                int glitchEnd = i + 1;
                while (glitchEnd < data.Count && Mathf.Abs(data[glitchEnd] - prev) >= threshold && glitchEnd - i <= glitchWidth + 1) {
                    ++glitchEnd;
                }
                if (glitchEnd >= data.Count) {
                    break;
                } else if (glitchEnd - i > 1 && glitchEnd - i <= glitchWidth + 1) {
                    while (i < glitchEnd) {
                        data[i] = prev;
                        ++i;
                    }
                }
                i = glitchEnd;
            }
        }

        public static List<int> GetDominantValues(this List<int> data, int count) {
            var mostDominant = data.GroupBy(i => i).OrderByDescending(grp => grp.Count())
                .Select(grp => grp.Key);
            count = Mathf.Min(count, mostDominant.Count());
            return mostDominant.Take(count).ToList();

        }

        public static List<int> GetAverageBuckets(this List<int> data, int numBuckets) {
            numBuckets = Mathf.Clamp(numBuckets, 1, data.Count);
            int bucketSize = data.Count / numBuckets;
            List<int> r = new List<int>(numBuckets);
            for (int i = 0; i < numBuckets; ++i) {
                int startIdx = i * bucketSize;
                int average = 0;
                for (int j = 0; j < bucketSize; ++j) {
                    average += data[startIdx + j];
                }
                r.Add(average / bucketSize);
            }
            return r;
        }

        public static float GetAverageDelta(this List<int> data) {
            if (data.Count == 0) {
                return int.MaxValue;
            }
            float avgDelta = 0;
            for (int i = 0; i + 1 < data.Count; ++i) {
                avgDelta += Mathf.Abs(data[i + 1] - data[i]);
            }
            return avgDelta /= data.Count;
        }

        public static float[] CombinedArray(this Queue<float[]> arrayList) {
            if (arrayList.Count == 1) {
                return arrayList.First();
            }
            float[] r = new float[arrayList.CombinedSize()];
            int dstIdx = 0;
            foreach (var arr in arrayList) {
                Array.Copy(arr, 0, r, dstIdx, arr.Length);
                dstIdx += arr.Length;
            }
            return r;
        }

        public static int CombinedSize(this Queue<float[]> arrayList) {
            return arrayList.Select(item => item.Length).Sum();
        }

        public static int HerzToMidi(float hz) {
            if (hz == 0f) {
                return UndefinedMidi;
            }
            return Mathf.RoundToInt(69f + 12f * Mathf.Log(hz / 440f, 2));
        }

        public static float HerzToFloatMidi(float hz) {
            return 69f + 12f * Mathf.Log(hz / 440f, 2);
        }

        public static float MidiToHerz(int midi) {
            return 440.0f * Mathf.Pow(2, ((float)midi - 69f) / 12f);
        }

        public static List<int> HerzToMidi(List<float> hzList) {
            return hzList.Select(item => HerzToMidi(item)).ToList();
        }

        public static string MidiToNote(int midi) {
            if (midi < 0)
                return "?";
            int octave = midi / 12 - 1;
            int note = midi % 12;
            return notes[note] + octave;

        }

        public static string NoteString(this List<int> data) {
            StringBuilder sb = new StringBuilder();
            foreach (var item in data) {
                var midiStr = (item != UndefinedMidi) ? item.ToString() : "?";
                sb.AppendFormat("{0} midi:{1} ({2:F1}Hz)\n", MidiToNote(item), midiStr, MidiToHerz(item));
            }
            return sb.ToString();
        }

    }
}