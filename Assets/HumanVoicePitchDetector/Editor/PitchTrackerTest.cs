using UnityEngine;
using UnityEditor;
using UnityEngine.TestTools;
using NUnit.Framework;
using System;
using System.Collections;
using System.IO;
using System.Linq;
using PitchDetector;

using System.Collections.Generic;

public class PitchTrackerTest {

	#if UNITY_EDITOR
	const AudioType preferredAudioType = AudioType.OGGVORBIS;
	const string audioFilePostfix = ".ogg";
	#else
	const AudioType preferredAudioType = AudioType.MPEG;
	const string audioFilePostfix = ".mp3";
	#endif

	[Test]
	public void PitchTrackerTestRemoveGlitch() {
		var testInput1 = new int[] { 1 };
		var testInput2 = new int[] { 1, 2, 1 };
		var testInput3 = new int[] { 1, 1, 1, 2, 3, 1, 1, 2, 1, 2, 2, 2 };
		var expected1 = new int[] { 1 };
		var expected2 = new int[] { 1, 1, 1 };
		var expected3 = new int[] { 1, 1, 1, 2, 3, 1, 1, 1, 1, 2, 2, 2 };
		var expected4 = new int[] { 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2 };
		var testList = testInput1.ToList ();
		testList.RemoveGlitches (1,1); // should not change anything
		CollectionAssert.AreEqual (expected1.ToList(), testList);
		testList = testInput2.ToList ();
		testList.RemoveGlitches (1,1); // should change index[1] to 1
		CollectionAssert.AreEqual (expected2.ToList(), testList);
		testList = testInput3.ToList ();
		testList.RemoveGlitches (1,1); // should change index[7] to 1
		CollectionAssert.AreEqual (expected3.ToList(), testList);
		testList = testInput3.ToList ();
		testList.RemoveGlitches (1,2); // should change index[3..4] and index[7] to 1
		CollectionAssert.AreEqual (expected4.ToList(), testList);
		testList = testInput3.ToList ();
		testList.RemoveGlitches (1,3); // same as previous: index[3..4] and index[7] from 2 to 1
		CollectionAssert.AreEqual (expected4.ToList(), testList);
	}

	[Test]
	public void TestListReverse() {
		var testInput1 = new float[] { 1f };
		var expected1 = new float[] { 1f };
		var testInput2 = new float[] { 1f, 2f };
		var expected2 = new float[] { 2f, 1f };
		var testInput3 = new float[] { 1f, 2f, 3f };
		var expected3 = new float[] { 3f, 2f, 1f };
		var testInput4 = new float[] { 1f, 2f, 3f, 4f };
		var expected4 = new float[] { 4f, 3f, 2f, 1f };
		var testInput5 = new float[] { 1f, 2f, 3f, 4f, 5f };
		var expected5 = new float[] { 5f, 4f, 3f, 2f, 1f };
		var expected6 = new float[] { 1f, 2f, 5f, 4f, 3f };
		var l = testInput1.ToList ();
		l.PartialReverse (0, 1);
		CollectionAssert.AreEqual (expected1.ToList(), l);
		l = testInput2.ToList ();
		l.PartialReverse (0, 2);
		CollectionAssert.AreEqual (expected2.ToList(), l);
		l = testInput3.ToList ();
		l.PartialReverse (0, 3);
		CollectionAssert.AreEqual (expected3.ToList(), l);
		l = testInput4.ToList ();
		l.PartialReverse (0, 4);
		CollectionAssert.AreEqual (expected4.ToList(), l);
		l = testInput5.ToList ();
		l.PartialReverse (0, 5);
		CollectionAssert.AreEqual (expected5.ToList(), l);
		l = testInput5.ToList ();
		l.PartialReverse (2, 3);
		CollectionAssert.AreEqual (expected6.ToList(), l);


	}

	[Test]
	public void PitchTrackerTestDominantValues() {
		var testInput = new int[] { 1, 1, 1, 2, 3, 1, 1, 2, 1, 2, 2, 2, 4, 4};
		var expected1 = new int[] { 1 };
		var expected2 = new int[] { 1, 2, 4 };
		var expected3 = new int[] { 1, 2, 4, 3 };
		var testList = testInput.ToList ().GetDominantValues (1);
		CollectionAssert.AreEqual (expected1.ToList(), testList);
		testList = testInput.ToList ().GetDominantValues (3);
		CollectionAssert.AreEqual (expected2.ToList(), testList);
		testList = testInput.ToList ().GetDominantValues (4);
		CollectionAssert.AreEqual (expected3.ToList(), testList);
		testList = testInput.ToList ().GetDominantValues (5);
		CollectionAssert.AreEqual (expected3.ToList(), testList);
	}


	[Test]
	public void PitchTrackerTestFrequencySweep() {
		// 0.001 and 0.01 clip lengths are too small
		float[] clipLengths = new float[] { 0.001f, 0.01f, 0.1f, 1f }; 
		float freqStep = 50f;
		float freq = 50f;
		float endFreq = 700f;
		float sampleRate = 22050;
		var pitchTracker = new RAPTPitchDetector (22050,40f,800f);
		float min, max, mean, median;
		float db = 0f;
		while (freq < endFreq) {
			foreach (var clipLen in clipLengths) {
				var inputSignal = GenerateWaveform (freq, clipLen, sampleRate);
				try {
					var pitchContour = pitchTracker.getPitch (inputSignal, ref db, sampleRate);
					if (pitchContour.Count > 0) {
						GetStats (pitchContour, out mean, out median, out min, out max);
						Assert.Greater(0.05, (max - min) / max); // verify that pitch values are consistent
						Debug.Log ("input:" + freq + "Hz measured:" + max + "Hz delta:" + 100 * Mathf.Abs(freq - max) / freq + "%");
						Assert.Greater(0.05, Mathf.Abs(freq - max) / freq); // verify measured pitch is same as input
					} else {
						Debug.Log ("Could not detect pitch. Clip length:" + clipLen + "s freq:" + freq + "Hz");
						Assert.Fail();
					}
				} catch (System.ArgumentException e) {
					// do not assert here since this is expected for small clip Lengths
					Debug.Log (e.Message);
				}
			}
			freq += freqStep;
		}
	}
		
	// Load audio clips from TestData directory and analyzes the midi notes of the clip
	[UnityTest]
	public IEnumerator PitchTrackerTestClips() {
		var pitchTracker = new RAPTPitchDetector ();
		// audioclip name and midi notes that are played in each clip
		Dictionary<string, int[]> testClips = new Dictionary<string, int[]>() { 
			{"piano", new int[] { 60, 62, 64, 65, 67, 69, 71, 72 }}, // C Major
			{"childvoice", new int[] { 60, 62, 64, 65, 67, 69, 71, 72 }}, // C Major
			{"violin", new int[]{ 55, 57, 59, 60, 62, 64, 66, 67, 69, 71, 72, 74, 76, 78}} // G Major (2 octaves)
		};
		foreach (var item in testClips) {
			var path = Path.Combine (Application.dataPath, "HumanVoicePitchDetector/TestData/" + item.Key + audioFilePostfix);
			using (WWW www = new WWW ("file://" + path)) {
				Debug.Log ("load clip:" + path);
				while (!www.isDone) {
					yield return null; 
				}
				if (string.IsNullOrEmpty (www.error)) {
					var clip = www.GetAudioClip (false, false, preferredAudioType);
					var midis = GetMidisFromClip (clip, pitchTracker);
					midis.RemoveGlitches (1, 2);
					var mostCommonNotes = midis.GetDominantValues (item.Value.Length);
					mostCommonNotes.Sort ();
					Debug.Log (mostCommonNotes.NoteString ());
                    CollectionAssert.AreEqual(item.Value.ToList(), mostCommonNotes);
                } else {
                    Debug.LogError("Could not load file:" + path);
                }

			}
		}

	}


	[UnityTest]
	public IEnumerator PitchTrackerTestStreaming() {
		var pitchTracker = new RAPTPitchDetector ();
		// audioclip name and midi notes that are played in each clip
		Dictionary<string, int[]> testClips = new Dictionary<string, int[]>() { 
			{"piano", new int[] { 60, 62, 64, 65, 67, 69, 71, 72 }}, // C Major
			{"childvoice", new int[] { 60, 62, 64, 65, 67, 69, 71, 72 }}, // C Major
			{"violin", new int[]{ 55, 57, 59, 60, 62, 64, 66, 67, 69, 71, 72, 74, 76, 78}} // G Major (2 octaves)
		};
		foreach (var item in testClips) {
            var path = Path.Combine (Application.dataPath, "HumanVoicePitchDetector/TestData/" + item.Key + audioFilePostfix);
			using (WWW www = new WWW ("file://" + path)) {
				Debug.Log ("load clip:" + path);
				while (!www.isDone) {
					yield return null; 
				}
				if (string.IsNullOrEmpty (www.error)) {
					var clip = www.GetAudioClip (false, false, preferredAudioType);
					var midis = GetMidisFromClipStreaming (clip, pitchTracker);
					midis.RemoveGlitches (1, 2);
					var mostCommonNotes = midis.GetDominantValues (item.Value.Length);
					mostCommonNotes.Sort ();
					Debug.Log (mostCommonNotes.NoteString ());
					CollectionAssert.AreEqual (item.Value.ToList(), mostCommonNotes);
                } else {
                    Debug.LogError("Could not load file:" + path);
                }

			}
		}

	}

	void GetStats(System.Collections.Generic.List<float> input, out float mean, out float median, out float min, out float max) {
		mean = input.Average ();
		min = input.Min ();
		max = input.Max ();
		median = input.OrderBy (sample => sample).ToList ().ElementAt (input.Count / 2);

	}

	List<int> GetMidisFromClip(AudioClip clip, RAPTPitchDetector tracker) {
		float[] audioSamples = new float[clip.samples]; 
		clip.GetData (audioSamples, 0);
		float db = 0f;
		var pitchContour = tracker.getPitch (audioSamples, ref db, (float)clip.frequency);
		return RAPTPitchDetectorExtensions.HerzToMidi (pitchContour);
	}

	List<int> GetMidisFromClipStreaming(AudioClip clip, RAPTPitchDetector tracker) {
		float[] audioSamples = new float[clip.samples];
		clip.GetData (audioSamples, 0);
		int frameSize = (int)(0.1f * (float)clip.frequency);
		var r = new List<float> ();
		int start = 0;
		int end = start + frameSize;
		float db = 0f;
		while (end < audioSamples.Length) {
			int prevEnd = end;
			r.AddRange(tracker.getPitch (audioSamples, start, ref end, ref db, (float)clip.frequency, true, false));
			if (start == end) {
				end = prevEnd + frameSize;
			} else {
				start = end;
				end += frameSize;
			}
		}

		if (audioSamples.Length - start > tracker.RequiredSampleCount ((float)clip.frequency)) {
			end = audioSamples.Length;
			r.AddRange (tracker.getPitch (audioSamples, start, ref end, ref db, (float)clip.frequency, true, true));
		}

		return RAPTPitchDetectorExtensions.HerzToMidi (r);
	}



	float[] GenerateWaveform(float freq, float lengthSec, float sampleRate) {
		int samples = (int)(lengthSec * sampleRate);
		float[] r = new float[samples];

		for (int i = 0; i < r.Length; i++) {
			r[i] = Mathf.Sin((2 * Mathf.PI * i * freq) / sampleRate);
		}
		return r;
	}
}
