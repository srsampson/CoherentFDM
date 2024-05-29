# CoherentFDM
A Coherent Frequency-Division Multiplex Voice Modem in Java

This was an exercise to see if I could write a DSP application in Java that operated in real-time.
The software was converted from the Codec 2 Digital Voice and Modem programs, with the goal of being native Java.

To test the modem:
Usage: java -jar FDMTest.jar hts.raw

The modem-out.raw is the test output that can be imported into Audacity (8000 sample rate, 1 channel)   
modem-demod.raw is the test output from decoding the modem output into voice. Use Audacity at 8000 also.   
These two files are created when running the FDMTest.jar program.

These programs were last compiled with JDK17 and Apache Netbeans IDE 20
