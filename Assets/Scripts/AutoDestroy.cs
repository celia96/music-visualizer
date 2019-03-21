using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class AutoDestroy : MonoBehaviour {

	public static float Lifetime = 1f;
	// Use this for initialization
	void Start () {
		Invoke ("DelayedDestroy", Lifetime);
	}
	
	// Update is called once per frame
	void DelayedDestroy () {
		Destroy (gameObject);
	}
}
