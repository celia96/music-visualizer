# Music Visualizer and Synesthesia Simulation in VR
A VR application that simulates the condition synesthesia and to visualize music.

## Introduction
Music visualization has become increasingly popular with the rise of VR, which has allowed for a more immersive and 
interactive music experience. Synesthesia is a rare neurological condition that affects 1 in 2000 people in which 
those affected perceive sensory data through multiple senses. Many Synesthetes possess a form in which the pitch of 
any sound they hear appears as specific bursts of color in their vision. Thus, VR allows for an exciting opportunity 
to simulate this condition for those that do not have it.

## Problem Statement
When audio and visual stimuli are coordinated, it creates a more immersive music experience in which the user not only 
hears music and sounds, but “sees” them. Listening to music is shown to increase levels of endorphins in the brain, 
lower cortisol levels and boost mood. Adding an ex- tra layer of stimulation by transforming the music listening experience 
into VR boosts these mind enhancing benefits even more. With this in mind, we created an application that stimulates the senses 
by synchronizing components of audio signals with visuals that dance across the environment in an aesthetic and captivating way.

## Methods
### Unity
We used Unity as a platform to build the 3D environment. Skybox and terrain components were used to set up the customizable backgrounds 
and particle systems, which simulate fluid entities by animating a number of 2D sprites, in order to create 4 different visual effects.

### Oculus Go
We used Oculus as our main platform to deploy our VR application. Unity provides built-in VR support for Oculus Go. 
We imported Oculus Unity Integration package that contains prefabs, scripts for camera behavior and API for controllers to supplement 
Unity’s built-in support. The application is built in Android since Oculus Go supports Android applications.

### Anaylzing Audio
We first converted the microphone input into the audio spectrum data using the 512 samples and BlackmanHarris FFT window. 
Once the sound is detected, the frequency elements are populated through the array of 512, where each element represents a 
frequency resolution that is approximately 40 Hz. We then distributed the frequency elements to 16 spectrum objects that form 
a circle of frequency bars, in which each object represents the certain range of the frequencies. To find the dominant frequency,
we also find the max amplitude element and multiply its index by the frequency resolution.

### Mapping Pitch to Color
We converted the sound waves picked up by the microphone into discrete MIDI notes that determine the pitch of a sound. 
We mapped MIDI notes from lowest pitch hearable by the human ear to highest pitch hearable to the HSL spectrum when deciding 
what color the visuals produced would be.

### Spawning Visuals
In order to prevent the visuals from instantiating right in front of the users
and overlapping with each other, we created 18 different spawn areas. Once the pitch is detected, it will choose 2 random 
unoccupied positions from the list of spawn areas and generate the visuals.
![1](https://user-images.githubusercontent.com/33583168/62816395-4853a400-baec-11e9-8221-426f209f6df4.png)

## Workflow
Scene opens up to a ring of bars that respond to changes in frequency and volume. 
They surround the user who is placed in the center of the ring and is able to observe visuals spawning all around them.
![2](https://user-images.githubusercontent.com/33583168/62816450-67066a80-baed-11e9-9ed4-b0ee9279e405.png)
![3](https://user-images.githubusercontent.com/33583168/62816451-67066a80-baed-11e9-84d3-b97bbf636dfa.png)
![4](https://user-images.githubusercontent.com/33583168/62816452-67066a80-baed-11e9-9427-d1d157a8cff7.png)
![5](https://user-images.githubusercontent.com/33583168/62816453-67066a80-baed-11e9-9f11-c03fb988ae4a.png)

## Demo
Demo of the application can be found here: https://youtu.be/b-mC2WIxsm0
