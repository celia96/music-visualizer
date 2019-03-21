using System;
using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using PitchDetector;


public class MicrophoneDemo : MonoBehaviour {

	public MicrophonePitchDetector detector;
	public int maxMidi = 100;
	public int minMidi = 0;
	public GameObject noteIndicatorPrefab;
//	public GameObject noteTitlePrefab;

//	private float scrollSpeed = 1f;
	private float appTime = 0f;
	private float analysisTime = 0f;
	private const float secsPerScreen = 10f;
	private float dbThres = -40;

	private Queue<PitchTime> drawQueue = new Queue<PitchTime> ();
	// Class that holds pitch value and its duration
	private class PitchTime {
		public PitchTime(float pitch, float dt) {
			this.pitch = pitch;
			this.dt = dt;
		}
		public float pitch;
		public float dt;
	};

	public float DbThreshold {
		get {
			return dbThres;
		}
		set {
			dbThres = value;
		}
	}

	// set event handlers 
	void Start () {
		if (detector == null) {
			Debug.LogWarning ("Pitch detector not set in MicrophoneDemo");
			return;
		}
		#if LOG_PITCH
		detector.onPitchDetected.AddListener (LogPitch);
		#endif
		Debug.Log("hello");
		if (noteIndicatorPrefab != null) {
			Debug.Log("draw");
			detector.onPitchDetected.AddListener (DrawPitch);
		} else {
			Debug.LogWarning ("Note indicator not set in MicrophoneDemo");
		}

//		scrollSpeed = 2 * Camera.main.orthographicSize * Camera.main.aspect / secsPerScreen;
//		AutoDestroy.Lifetime = 0.5f * secsPerScreen;
//		PopulateNoteTitles ();
	}



	// Print pitch values to console
	public void LogPitch (List<float> pitchList, int samples, float db) {
		var midis = RAPTPitchDetectorExtensions.HerzToMidi (pitchList);
//		Debug.Log ("detected " + pitchList.Count + " values from " + samples + " samples, db:" + db);
		Debug.Log (midis.NoteString ());
	}

	// Set pitch values to draw queue
	public void DrawPitch(List<float> pitchList, int samples, float db) {
		// instead of drawing, do something else...
		// map the pitch values in the pitch list into colors
		var duration = (float)samples / (float)detector.micSampleRate;
//		Debug.Log ("draw pitch" + duration);
//		Debug.Log ("pitch list" + pitchList);
		if (pitchList.Count > 0 && db > dbThres) {
			foreach (var pitchVal in pitchList) {
//				Debug.Log ("pitch value" + pitchVal);
				drawQueue.Enqueue (new PitchTime (pitchVal, duration / pitchList.Count));
			}
		} else {
			drawQueue.Enqueue(new PitchTime(0f, duration));
		}
	}

	// Get pitch values from queue and draw on screen
	void Update() {
		if (!detector.Record) {
			return;
		}
		appTime += Time.deltaTime;
		while (analysisTime < appTime && drawQueue.Count > 0) {
			var item = drawQueue.Dequeue ();
			Debug.Log ("pitch item" + item.pitch);
			var midi = RAPTPitchDetectorExtensions.HerzToFloatMidi (item.pitch);
//			Debug.Log ("HUE" + convertToHue (midi));

//			Debug.Log("Color" + Color.HSVToRGB(m_Hue, m_Saturation, m_Value))
//			Debug.Log ("hi"+ RAPTPitchDetectorExtensions.HerzToFloatMidi (4186) + "+" + RAPTPitchDetectorExtensions.HerzToFloatMidi (27.5f));
			// 21 - 108
			if (!float.IsInfinity(midi)) {
				var hue = convertToHue (midi);
//				float y = MidiToScreenY (midi);
//				Debug.Log ("y" + y);
//				m_Renderer.material.color = Color.HSVToRGB(m_Hue, m_Saturation, m_Value);
//				float x = analysisTime * scrollSpeed;
				// instantiate 3D object corresponding to the pitch
//				GameObject newNote = Instantiate<GameObject> (noteIndicatorPrefab);
				Debug.Log ("hue" + hue);
				MeshRenderer m_Renderer = noteIndicatorPrefab.GetComponent<MeshRenderer> ();
				Debug.Log ("color" + Color.HSVToRGB (hue, 0.5f, 0.5f));
				m_Renderer.material.color = Color.HSVToRGB (hue, 1f, 0.8f);
//				
//				newNote.transform.position = new Vector3 (x, y);
//				newNote.transform.SetParent (transform, false);
			}
			analysisTime += item.dt;
		}

//		MoveLeft (scrollSpeed * Time.deltaTime);
	}
		
//	// Populate note titles to left side of the display
//	private void PopulateNoteTitles() {
//		if (noteTitlePrefab == null) {
//			Debug.LogWarning ("Note Title Prefab not set in MicrophoneDemo");
//			return;
//		}
//		float x = -Camera.main.orthographicSize * Camera.main.aspect;
//		for (int i = minMidi; i <= maxMidi; ++i) {
//			float y = MidiToScreenY ((float)i);
//			GameObject newTitle = Instantiate<GameObject> (noteTitlePrefab);
//			newTitle.transform.position = new Vector3 (x, y);
//			newTitle.GetComponent<TextMesh> ().text = RAPTPitchDetectorExtensions.MidiToNote (i);
//		}
//	}
//
//	private float MidiToScreenY(float midiVal) {
//		if (float.IsInfinity (midiVal)) {
//			midiVal = RAPTPitchDetectorExtensions.HerzToMidi (1f);
//		}
//		float viewHeight = 2 * Camera.main.orthographicSize;
//		float dy = viewHeight / (float)(maxMidi - minMidi);
//		float middleMidi = 0.5f * (float)(maxMidi + minMidi);
//		return dy * (midiVal - middleMidi);
//
//	}

	private float convertToHue(float midiVal) {
		float hue = 0f;
//		OldRange = (OldMax - OldMin)  35
//		NewRange = (NewMax - NewMin)  1
//		NewValue = (((OldValue - OldMin) * NewRange) / OldRange) + NewMin
		if (midiVal >= 21 && midiVal <= 108) {
			if (midiVal < 35) {
				hue = 0f;
			} else if (midiVal > 70) {
				hue = 1f;
			} else {
				// 35 - 70
				hue = (((midiVal - 35f) * 1f) / 35f);
			}
		}
		return hue;
	}

	// scroll game object to left so that new pitch values are drawn on the same position of the screen
	private void MoveLeft(float amount) {
		var tmp = transform.position;
		tmp.x -= amount;
		transform.position = tmp;
	}
}
