/* MicrophonePitchDetector.cs
*
* Copyright 2017 Outloud Oy
* Written by: mikko.k.koivisto@outloud.fi
* 
* Analyzes pitch of Microphone input. After analysis, onPitchDetected event is invoked.
*/
using System;
using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using UnityEngine.Events;

namespace PitchDetector {


    [System.Serializable]
    public class PitchEvent : UnityEvent<List<float>, int, float> {
    }

    public class MicrophonePitchDetector : MonoBehaviour {

        public PitchEvent onPitchDetected;
        public int micSampleRate = 16000;

        private RAPTPitchDetector pitchDetector;
        private bool recording;

        public bool Record {
            get {
                return Microphone.IsRecording(null);
            }
            set {
                if (value && !Microphone.IsRecording(null)) {
                    StartMic();
                } else if (!value && Microphone.IsRecording(null)) {
                    StopMic();
                }
            }

        }

        private RAPTPitchDetector Detector {
            get {
                if (pitchDetector == null) {
                    pitchDetector = new RAPTPitchDetector((float)micSampleRate, 50f, 800f);
                }
                return pitchDetector;
            }
        }

        public void ToggleRecord() {
            Record = !Record;
        }

        private void StartMic() {
            StartCoroutine(RecordingCoroutine());
        }


        private void StopMic() {
            recording = false;
        }


        private IEnumerator RecordingCoroutine() {
            recording = true;
            AudioClip rec = Microphone.Start(null, true, 1, micSampleRate);
            float[] readBuffer = new float[rec.samples];
            int recPos = 0;
            int prevRecPos = 0;
            Func<bool> enoughSamples = () => { 
                int count = (readBuffer.Length + Microphone.GetPosition(null) - prevRecPos) % readBuffer.Length;
                return count > Detector.RequiredSampleCount((float)micSampleRate);
            };
            while (recording) {
                prevRecPos = recPos;
                yield return new WaitUntil(enoughSamples);
                recPos = Microphone.GetPosition(null);
                rec.GetData(readBuffer, 0);
                int sampleCount = (readBuffer.Length + recPos - prevRecPos) % readBuffer.Length;
                float db = 0f;
                List<float> pitchValues = Detector.getPitch(readBuffer, prevRecPos, ref recPos, ref db, (float)micSampleRate, true, !recording);
                sampleCount = (readBuffer.Length + recPos - prevRecPos) % readBuffer.Length;
                if (sampleCount > 0) {
                    onPitchDetected.Invoke(pitchValues, sampleCount, db);
                }
            }
            Microphone.End(null);
        }

    }
}