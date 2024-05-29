/*
 * Coherent FDM Modem Tester Application
 *
 * Steve Sampson, K5OKC
 */
package test;

import complex.Complex;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;

public final class Main {

    private static API api;
    private static FileInputStream fin;
    private static FileInputStream ain;
    private static FileOutputStream fout;
    private static String filename_in;
    private static String filename_out;

    public static void main(String[] args) {
        short mod;
        int msb, lsb, i, lnin, nout;

        if (args.length != 1) {
            System.err.println("Usage: java -jar FDMTest.jar filename.raw");
            System.exit(-1);
        }

        try {
            ain = new FileInputStream(args[0]);
        } catch (FileNotFoundException e0) {
            System.err.println("Fatal - Couldn't open audio file: " + e0.toString());
            System.exit(-1);
        }

        api = new API();
        api.setSNRSquelchThreshold(-100.0f);
        api.setSquelchBoolean(true);

        // create filenames in the current directory
        try {
            filename_in = "modem-demod.raw";
            filename_out = "modem-out.raw";
            fout = new FileOutputStream(filename_out);
        } catch (IOException e0) {
            System.err.println("Fatal - Couldn't open work file: " + e0.toString());
            System.exit(-1);
        }

        int n_speech_samples = api.getSpeechSamples();          // 640
        int n_nom_modem_samples = api.getNominalModemSamples();    // 640 (with 8 kHz conversion)
        short[] speech_in = new short[n_speech_samples];            // PCM 16-bit audio
        Complex[] mod_out = new Complex[n_nom_modem_samples];       // Complex Modulated Signal

        for (i = 0; i < n_nom_modem_samples; i++) {
            mod_out[i] = new Complex();
        }

        // OK, write a modem RAW file of the input PCM speech
        try {
            // Process audio in groups of two frames (320 * 2 = 640 samples * 2 bytes)
            while (ain.available() >= (n_speech_samples * 2)) {
                // Audio is little-endian PCM Signed 16-bit
                for (i = 0; i < n_speech_samples; i++) {
                    lsb = ain.read() & 0xFF;
                    msb = ain.read() & 0xFF;
                    speech_in[i] = (short) (msb << 8 | lsb);
                }

                // Send the Audio slice pair to the transmit API
                // We expect a complex signal back
                nout = api.send(mod_out, speech_in);

                // Take the real part of the signal, convert to signed 16-bits PCM
                // and write this to a file as little-endian, 8 kHz sample rate
                for (i = 0; i < nout; i++) {
                    mod = (short) mod_out[i].real();
                    lsb = mod & 0xFF;
                    msb = (mod >>> 8) & 0xFF;

                    fout.write(lsb);
                    fout.write(msb);
                }
            }

            fout.close();
        } catch (IOException e1) {
            System.err.println("Encoding modem file: " + e1.toString());
            System.exit(-1);
        }

        // OK, now see if we can decode the modem RAW file just written above

        try {
            fin = new FileInputStream(filename_out);
            fout = new FileOutputStream(filename_in);
        } catch (FileNotFoundException e2) {
            System.err.println("Fatal - Couldn't open work files: " + e2.toString());
            System.exit(-1);
        }

        Complex[] signal = new Complex[n_nom_modem_samples];
        short[] speech_out = new short[n_speech_samples];
        boolean[] sync = new boolean[1];    // pointer
        float[] snr = new float[1];       // pointer
        int frame = 0;

        try {
            lnin = api.getNIN();

            while (fin.available() >= (lnin * 2)) {  // bytes * 2 for 640 shorts
                frame++;

                for (i = 0; i < lnin; i++) {
                    lsb = fin.read();
                    msb = fin.read();
                    mod = (short) (msb << 8 | lsb);
                    signal[i] = new Complex((float) mod, 0.0f);
                }

                nout = api.receive(speech_out, signal);
                lnin = api.getNIN();

                for (i = 0; i < nout; i++) {
                    lsb = speech_out[i] & 0xFF;
                    msb = (speech_out[i] >>> 8) & 0xFF;

                    fout.write(lsb);
                    fout.write(msb);
                }

                api.getModemStats(sync, snr);

                System.out.printf("Frame: %d Nout: %d Sync: %s SNR: %.4f dB%n",
                        frame, nout, sync[0], snr[0]);
            }

            fin.close();
            fout.close();
        } catch (IOException e3) {
            System.err.println("Fatal - Decoding modem file: " + e3.toString());
            System.exit(-1);
        }
    }
}