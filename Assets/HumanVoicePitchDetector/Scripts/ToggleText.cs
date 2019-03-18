using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using UnityEngine.UI;

public class ToggleText : MonoBehaviour {

	public string[] texts = new string[] { "REC", "STOP" };

	private Text btnText;
	private int textIdx;
	// Use this for initialization
	void Start () {
		textIdx = 0;
		btnText = GetComponentInChildren<Text> ();
		SetText ();
	}

	public void Toggle() {
		if (texts.Length > 0) {
			textIdx += 1;
			textIdx = textIdx % texts.Length;
			SetText ();
		}
	}

	// Update is called once per frame
	private void SetText () {
		if (btnText != null && texts != null && texts.Length > 0) {
			btnText.text = texts[textIdx];
		}
	}
}
