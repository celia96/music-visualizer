using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class TextMeshSharpener : MonoBehaviour {

	// Use this for initialization
	void Start () {
		var textMesh = GetComponent<TextMesh> ();
		float pixelRatio = (Camera.main.orthographicSize * 2.0f) / Camera.main.pixelHeight;
		float fuzziness = 32.0f; // lower number is better

		int previousFontSize = textMesh.fontSize;
		textMesh.fontSize = Mathf.RoundToInt(textMesh.fontSize / (pixelRatio * fuzziness));

		float ratio = previousFontSize / ((float) textMesh.fontSize);
		transform.localScale *= ratio;
	}
	

}
