# CoherentFDM
A Coherent Frequency-Division Multiplex Voice Modem in Java

This was an exercise to see if I could write a DSP application in Java that operated in real-time.
The software was converted from the Codec 2 Digital Voice and Modem programs, with the goal of being native Java.

To test the modem:
Usage: java -jar FDMTest.jar hts.raw

The modem-out.raw is the test output that can be imported into Audacity (8000 sample rate, 1 channel)   
modem-demod.raw is the test output from decoding the modem output into voice. Use Audacity at 8000 also.   
These two files are created when running the FDMTest.jar program.

These programs were last compiled with Java SE 21.0.2 and Apache Netbeans IDE 21
#### Theory
The modem sends and receives a row of subcarriers 75 times a second. However, it takes six of these rows to make up a modem frame. First, two pilot reference-phase rows (28 bits), then two speech vocoder rows (28 bits), and finally two more rows for the second speech vocoder frame (28 bits). The process then repeats as long as the transmitter Push-To-Talk (PTT) is keyed.

Thus, a modem frame is 84 bits total. 56 bits are used for speech, and 28 bits are used for the reference-phase pilots. These pilots are what makes this a coherent modem. They are used to correct the received data bit phases. The data rate is 1050 bit/s (75 Baud × 14 bits). The effective data rate is 700 bit/s (75 Baud / 6 or 12.5 Baud × 56 bits). Each row of 14 bits is sent as seven QPSK carriers (2 bits per carrier).

The modem timings are also relevant, in that each speech vocoder frame outputs 28 bits every 40 ms. Since the modem has an 80 ms modem frame, it can transport two speech vocoder frames.

There are 100 complex IQ (In-Phase and Quadrature-Phase) audio samples for each row, at a 7500 Hz rate. 600 samples total for the modem frame. Thus, 100×6 * 12.5 equals the 7500 Hz sample rate. Using a rate conversion filter, the application is provided an 8 kHz interface, which is much more compatible with sound cards. There are 640 complex audio samples at the 8 kHz rate. This rate conversion would not be necessary in firmware.

The modem operates with a center frequency of 1500 Hz. The initial FDM subcarrier frequencies are set using a spreading function. This changes the spacing of each subcarrier a little bit more each subcarrier further to the left. About 105 Hz apart on the right, to about 109 Hz apart on the left. This design, along with spectrum clipping, improves the Peak to Average Power Ratio (PAPR). The measured Crest factor is about 8.3 dB with clipping, and about 10.3 dB without clipping.

The modem waveform consumes a different amount of bandwidth, depending on whether the diversity channel is enabled. About 750 Hz per group of seven subcarriers. Normally you would want to use diversity on shortwave, but optionally on VHF and above.
