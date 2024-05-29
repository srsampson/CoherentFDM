/*
 * Copyright (C) 2015 David Rowe
 *
 * All rights reserved.
 *
 * Java 17 Version by Steve Sampson
 */
package test;

import codec2.Codec2;
import cohfdm.CohfdmRX;
import cohfdm.CohfdmTX;
import cohfdm.Resample;
import cohfdm.IDefines;
import complex.Complex;
import complex.ComplexMath;

public final class API implements IDefines {

    private static final float NORM_PWR = 1.74f;     // power output fudge factor
    private static final float MODEM_SCALE = 1000.0f;
    //
    private final int mspeechSamples;
    private final int msamplesPerCodecFrame;
    private final int mnominalModemSamples; // size of tx and most rx modem sample buffers

    private final int mbitsPerCodecFrame;
    private final int mbitsPerModemFrame;
    private final int mbytesPerCodecFrame;

    private final byte[] mpackedCodecBits;
    private final boolean[] mcodecBits;

    private float msnrSquelchThreshold;
    private boolean msquelchEnable;
    private int mfnin;

    private final CohfdmTX cohfdmTX;
    private final CohfdmRX cohfdmRX;
    private final Codec2 codec2;
    private final Resample Filter7500;
    private final Resample Filter8000;

    public API() {
        cohfdmRX = new CohfdmRX();
        cohfdmTX = new CohfdmTX();

        codec2 = new Codec2();
        Filter7500 = new Resample();
        Filter8000 = new Resample();

        mpackedCodecBits = new byte[4];
        mcodecBits = new boolean[COHPSK_BITS_PER_FRAME];                  // 56
        mfnin = 0;
        mnominalModemSamples = COHPSK_NOM_TX_SAMPLES_PER_FRAME * 8000 / 7500;   // 640

        mbitsPerCodecFrame = codec2.codec2_getBitsPerFrame();           // 28
        mbitsPerModemFrame = COHPSK_BITS_PER_FRAME;            // 56
        mbytesPerCodecFrame = (mbitsPerCodecFrame + 7) / 8;    // 4

        msamplesPerCodecFrame = codec2.codec2_getSamplesPerFrame();     // 320
        mspeechSamples = 2 * msamplesPerCodecFrame;              // 320 * 2 = 640
    }

    public int send(Complex[] fdm, short[] speech_in) {
        Complex[] tx_fdm = new Complex[COHPSK_NOM_TX_SAMPLES_PER_FRAME];
        int bit, nbyte, index, i, j;

        for (i = 0; i < COHPSK_NOM_TX_SAMPLES_PER_FRAME; i++) {
            tx_fdm[i] = new Complex();
        }

        index = 0;

        // Send the audio frame pair to the Vocoder to convert to digital
        for (j = 0; j < mbitsPerModemFrame; j += mbitsPerCodecFrame) {
            codec2.codec2_encode(mpackedCodecBits, index, speech_in);

            // Get ready for the second frame of the pair
            index += msamplesPerCodecFrame;

            /* unpack bits, MSB first */
            bit = 7;
            nbyte = 0;

            for (i = 0; i < mbitsPerCodecFrame; i++) {
                mcodecBits[j + i] = ((mpackedCodecBits[nbyte] >> bit) & 0x1) != 0;
                bit--;
                if (bit < 0) {
                    bit = 7;
                    nbyte++;
                }
            }
        }

        // Modulate the two digital frames using the PSK modulator
        // Result will be a complex waveform in tx_fdm
        cohfdmTX.modulate(tx_fdm, mcodecBits);

        // Amplify the signal in the complex array
        for (i = 0; i < COHPSK_NOM_TX_SAMPLES_PER_FRAME; i++) {
            fdm[i] = ComplexMath.times(tx_fdm[i], MODEM_SCALE * NORM_PWR);
        }

        // Convert the sample rate from the modulators 7500 Hz to standard 8000 Hz
        // returns number of samples
        return Filter7500.interpDecim(fdm, COHPSK_NOM_TX_SAMPLES_PER_FRAME, 16, 15);
    }

    public void getModemStats(boolean[] sync, float[] snr) {
        sync[0] = cohfdmRX.getSync();
        snr[0] = cohfdmRX.getSNR();
    }

    public int getNIN() {
        return (16 * mfnin + Filter8000.getDecimationIndex()) / 15;
    }

    public int receive(short[] speech_out, Complex[] signal) {
        boolean[] rx_bits = new boolean[COHPSK_BITS_PER_FRAME];
        boolean[] sync = new boolean[1];
        int i, j, bit, nbyte, index, nout;

        // echo samples back out as default (say if sync not found)

        for (i = 0; i < getNIN(); i++) {
            speech_out[i] = (short) signal[i].real();
        }

        mfnin = Filter8000.interpDecim(signal, getNIN(), 15, 16);

        for (i = 0; i < mfnin; i++) {
            signal[i] = ComplexMath.divide(signal[i], MODEM_SCALE);      // 0..599 samples
        }

        mfnin = cohfdmRX.demodulate(rx_bits, sync, signal, mfnin);

        if (sync[0] == true) {
            index = 0;

            for (j = 0; j < COHPSK_BITS_PER_FRAME; j += mbitsPerCodecFrame) {
                bit = 7;
                nbyte = 0;

                for (i = 0; i < mbytesPerCodecFrame; i++) {
                    mpackedCodecBits[i] = 0;
                }

                for (i = 0; i < mbitsPerCodecFrame; i++) {
                    mpackedCodecBits[nbyte] |= ((rx_bits[j + i] ? 1 : 0) << bit);
                    bit--;
                    if (bit < 0) {
                        bit = 7;
                        nbyte++;
                    }
                }

                codec2.codec2_decode(speech_out, index, mpackedCodecBits);

                if ((msquelchEnable == true) && (Double.compare(cohfdmRX.getSNR(), msnrSquelchThreshold) < 0)) {
                    // get rid of both pairs
                    for (i = 0; i < mspeechSamples; i++) {
                        speech_out[i] = (short) 0;
                    }
                }

                index += msamplesPerCodecFrame;
            }

            nout = mspeechSamples;
        } else {
            // no sync

            nout = getNIN();

            if (msquelchEnable == true) {
                for (i = 0; i < nout; i++) {
                    speech_out[i] = (short) 0;
                }
            }
        }

        return nout;
    }

    public void setSquelchBoolean(boolean val) {
        msquelchEnable = val;
    }

    public void setSNRSquelchThreshold(float val) {
        msnrSquelchThreshold = val;
    }

    public int getSamplesperCodecFrame() {
        return msamplesPerCodecFrame;
    }

    public int getSpeechSamples() {
        return mspeechSamples;
    }

    public int getNominalModemSamples() {
        return mnominalModemSamples;
    }
}