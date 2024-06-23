/*
 * Copyright (C) 1993-2017 David Rowe
 *
 * All rights reserved
 *
 * Licensed under GNU LGPL V2.1
 * See LICENSE file for information
 *
 * Java Version by Steve Sampson
 */
package codec2;

import complex.Complex;

public final class Codec2 implements IDefines {

    private final Model m_encodeModel;
    private final Model[] m_decodeModels;
    //
    private final FFT m_fftencode;
    private final FFT m_fftdecode;
    private final FFT m_fftphase;
    private final Pack m_bitEncode;
    private final Pack m_bitDecode;
    private final Analyze m_analyze;
    private final Synthesize m_synthesize;
    private final Amp m_amp;

    private boolean m_equalizer;

    public Codec2() {
        m_fftdecode = new FFT(FFT_SIZE);     // for full duplex
        m_fftencode = new FFT(FFT_SIZE);     // we need separate FFT's
        m_fftphase = new FFT(AMP_PHASE_NFFT);

        m_bitEncode = new Pack();    // same for pack
        m_bitDecode = new Pack();    // independant for duplex

        m_synthesize = new Synthesize(m_fftdecode);
        m_analyze = new Analyze(m_fftencode);
        m_amp = new Amp(m_fftphase);

        m_equalizer = false;          // defaults to off

        /*
         * Pre-allocate the four voice models for decode
         */
        m_decodeModels = new Model[4];

        for (int i = 0; i < 4; i++) {
            m_decodeModels[i] = new Model();
        }

        m_encodeModel = new Model();
    }

    public int codec2_getAMPK() {
        return m_amp.getAMPK();
    }

    public boolean codec2_getEQBoolean() {
        return m_equalizer;
    }

    public void codec2_setEQBoolean(boolean val) {
        m_equalizer = val;
    }

    public float[] codec2_getEQValues() {
        return m_amp.getEQ();
    }

    public float codec2_getVariance() {
        long count = m_amp.getSquaredCount();
        
        if (count != 0L) {
            return m_amp.getSquaredError() / (float) count;
        } else {
            return 0.0f;
        }
    }

    public int codec2_getBitsPerFrame() {
        return 28;
    }

    public int codec2_getBytesPerFrame() {
        return (codec2_getBitsPerFrame() + 7) / 8;
    }

    public int codec2_getSamplesPerFrame() {
        return 320;
    }

    /*
     * The speech output as packed bytes from the 16-bit PCM buffer array
     *
     * The byte array is (codec2_getBitsPerFrame() + 7) / 8 bytes long
     */
    public void codec2_encode(byte[] bits, int index, short[] speech) {
        encode(bits, index, speech);
    }

    /*
     * The speech output as packed bytes from the 16-bit PCM buffer array
     * No offset version
     */
    public void codec2_encode(byte[] bits, short[] speech) {
        encode(bits, 0, speech);
    }

    /*
     * The output of a PCM 16-bit buffer array from an input of packed bytes
     */
    public void codec2_decode(short[] speech, int index, byte[] bits) {
        decode(speech, index, bits);
    }

    /*
     * The output of a PCM 16-bit buffer array from an input of packed bytes
     * No offset version.
     */
    public void codec2_decode(short[] speech, byte[] bits) {
        decode(speech, 0, bits);
    }

    /*
     * This provides the user a way to reset everything back to
     * defaults for the current mode.
     */
    public void codec2_reset() {
        m_bitEncode.reset();
        m_bitDecode.reset();
        m_analyze.reset();
        m_synthesize.reset();
        m_encodeModel.reset();

        m_amp.resetEQ();  // reset the VQ equalizer

        for (int i = 0; i < 4; i++) {
            m_decodeModels[i].reset();
        }
    }

    /*
     * Encodes 320 speech samples (40ms of speech) into 28 bits. Bits are shifted
     * left into the bytes.
     *
     * In the Coherent Modem, this is called with a 0, or a 320 sample
     * offset index. As the modem processes two codec frames for each modem
     * frame.
     */
    private void encode(byte[] bits, int offset, short[] speech) {
        int[] indexes = new int[4];

        for (int i = 0; i < codec2_getBytesPerFrame(); i++) {
            bits[i] = 0;
        }

        m_bitEncode.reset();

        m_analyze.analyze_one_frame(m_encodeModel, speech, offset);
        m_analyze.analyze_one_frame(m_encodeModel, speech, offset + N_SAMP);
        m_analyze.analyze_one_frame(m_encodeModel, speech, offset + N_SAMP * 2);
        m_analyze.analyze_one_frame(m_encodeModel, speech, offset + N_SAMP * 3);

        m_amp.amp_model_to_indexes(m_encodeModel, indexes, m_equalizer);

        m_bitEncode.pack(bits, indexes[0], 9);
        m_bitEncode.pack(bits, indexes[1], 9);
        m_bitEncode.pack(bits, indexes[2], 4);
        m_bitEncode.pack(bits, indexes[3], 6);
    }

    private void decode(short[] speech, int offset, byte[] bits) {
        Complex[][] HH = new Complex[4][MAX_AMP + 1];
        int[] indexes = new int[4];

        m_bitDecode.reset();

        indexes[0] = m_bitDecode.unpack(bits, 9);
        indexes[1] = m_bitDecode.unpack(bits, 9);
        indexes[2] = m_bitDecode.unpack(bits, 4);
        indexes[3] = m_bitDecode.unpack(bits, 6);

        for (int i = 0; i < 4; i++) {
            for (int j = 0; j <= MAX_AMP; j++) {
                HH[i][j] = new Complex();
            }
        }

        m_amp.amp_indexes_to_model(m_decodeModels, HH, indexes);

        for (int i = 0; i < 4; i++) {
            m_synthesize.synthesize_one_frame(m_decodeModels[i], speech, offset + N_SAMP * i, HH[i]);
        }
    }

    public float codec2_energy(byte[] bits) {
        int[] indexes = new int[4];

        m_bitDecode.reset();

        indexes[0] = m_bitDecode.unpack(bits, 9);
        indexes[1] = m_bitDecode.unpack(bits, 9);
        indexes[2] = m_bitDecode.unpack(bits, 4);
        indexes[3] = m_bitDecode.unpack(bits, 6);

        float mean = CODES0[indexes[2]] - 10.0f;

        if (indexes[3] == 0) {
            mean -= 10.0f;
        }

        return (float) Math.pow(10.0, mean / 10.0);
    }
}
