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
import complex.ComplexMath;

public final class Analyze implements IDefines {

    private static final float VOICING_THRESHOLD_DB = 6.0f;
    private static final float R = TAU / FFT_SIZE;
    private static final float WO_MAX = TAU / P_MIN;
    private static final float WO_MIN = TAU / P_MAX;
    private static final float CNLP = 0.3f;         // post processor constant
    private static final float COEFF = 0.95f;       // notch filter parameter
    //
    private static final int M_PITCH = 320;
    private static final int NW = 280;
    private static final int DEC = 5;               // decimation factor
    //
    private static final int MIN_BIN = FFT_SIZE * DEC / P_MAX;   // 16
    private static final float PREV_BIN = (4000.0f / (float) Math.PI) * (float) (FFT_SIZE * DEC) / FS;
    //
    private final FFT m_fft;
    //
    private final Complex[] m_fw;             // DFT of squared signal (input)
    //
    private final float[] m_mem_fir;         // decimation FIR filter memory
    private final float[] m_sq;              // squared speech samples
    private final float[] m_coswintab;       // optimization
    private final float[] m_Sn;
    //
    private float m_mem_x;
    private float m_mem_y;
    //
    private float m_previous_Wo;

    /*
     * 48 tap 600Hz low pass FIR filter coefficients
     */
    private static final float[] NLPCOEF = {
        -0.0010818124f,
        -0.0011008344f,
        -0.00092768838f,
        -0.00042289438f,
        0.0005503419f,
        0.0020029849f,
        0.0037058509f,
        0.0051449415f,
        0.0055924666f,
        0.0043036754f,
        0.00080284511f,
        -0.004820461f,
        -0.01170581f,
        -0.018199275f,
        -0.022065282f,
        -0.02092061f,
        -0.012808831f,
        0.0032204775f,
        0.026683811f,
        0.055520624f,
        0.086305944f,
        0.11480192f,
        0.13674206f,
        0.14867556f,
        0.14867556f,
        0.13674206f,
        0.11480192f,
        0.086305944f,
        0.055520624f,
        0.026683811f,
        0.0032204775f,
        -0.012808831f,
        -0.02092061f,
        -0.022065282f,
        -0.018199275f,
        -0.01170581f,
        -0.004820461f,
        0.00080284511f,
        0.0043036754f,
        0.0055924666f,
        0.0051449415f,
        0.0037058509f,
        0.0020029849f,
        0.0005503419f,
        -0.00042289438f,
        -0.00092768838f,
        -0.0011008344f,
        -0.0010818124f
    };

    private static final float[] HAMMING = {
        0.0f,
        0.0f,
        0.0f,
        0.0f,
        0.0f,
        0.0f,
        0.0f,
        0.0f,
        0.0f,
        0.0f,
        0.0f,
        0.0f,
        0.0f,
        0.0f,
        0.0f,
        0.0f,
        0.0f,
        0.0f,
        0.0f,
        0.0f,
        0.0f,
        0.000001f,
        0.000002f,
        0.000005f,
        0.000009f,
        0.000014f,
        0.000020f,
        0.000027f,
        0.000035f,
        0.000044f,
        0.000055f,
        0.000066f,
        0.000078f,
        0.000092f,
        0.000106f,
        0.000122f,
        0.000139f,
        0.000156f,
        0.000175f,
        0.000195f,
        0.000215f,
        0.000237f,
        0.000260f,
        0.000283f,
        0.000308f,
        0.000333f,
        0.000360f,
        0.000387f,
        0.000415f,
        0.000445f,
        0.000475f,
        0.000505f,
        0.000537f,
        0.000570f,
        0.000603f,
        0.000637f,
        0.000672f,
        0.000708f,
        0.000744f,
        0.000781f,
        0.000819f,
        0.000857f,
        0.000896f,
        0.000936f,
        0.000977f,
        0.001018f,
        0.001059f,
        0.001101f,
        0.001144f,
        0.001187f,
        0.001231f,
        0.001275f,
        0.001320f,
        0.001365f,
        0.001410f,
        0.001456f,
        0.001502f,
        0.001549f,
        0.001595f,
        0.001642f,
        0.001690f,
        0.001737f,
        0.001785f,
        0.001833f,
        0.001881f,
        0.001930f,
        0.001978f,
        0.002027f,
        0.002075f,
        0.002124f,
        0.002172f,
        0.002221f,
        0.002270f,
        0.002318f,
        0.002367f,
        0.002415f,
        0.002463f,
        0.002511f,
        0.002559f,
        0.002607f,
        0.002655f,
        0.002702f,
        0.002749f,
        0.002795f,
        0.002842f,
        0.002888f,
        0.002933f,
        0.002979f,
        0.003023f,
        0.003068f,
        0.003112f,
        0.003155f,
        0.003198f,
        0.003240f,
        0.003282f,
        0.003324f,
        0.003364f,
        0.003404f,
        0.003444f,
        0.003483f,
        0.003521f,
        0.003558f,
        0.003595f,
        0.003631f,
        0.003666f,
        0.003701f,
        0.003734f,
        0.003767f,
        0.003799f,
        0.003831f,
        0.003861f,
        0.003891f,
        0.003919f,
        0.003947f,
        0.003974f,
        0.004000f,
        0.004025f,
        0.004049f,
        0.004072f,
        0.004094f,
        0.004116f,
        0.004136f,
        0.004155f,
        0.004173f,
        0.004190f,
        0.004206f,
        0.004222f,
        0.004236f,
        0.004249f,
        0.004261f,
        0.004271f,
        0.004281f,
        0.004290f,
        0.004298f,
        0.004304f,
        0.004310f,
        0.004314f,
        0.004317f,
        0.004319f,
        0.004320f,
        0.004320f,
        0.004319f,
        0.004317f,
        0.004314f,
        0.004310f,
        0.004304f,
        0.004298f,
        0.004290f,
        0.004281f,
        0.004271f,
        0.004261f,
        0.004249f,
        0.004236f,
        0.004222f,
        0.004206f,
        0.004190f,
        0.004173f,
        0.004155f,
        0.004136f,
        0.004116f,
        0.004094f,
        0.004072f,
        0.004049f,
        0.004025f,
        0.004000f,
        0.003974f,
        0.003947f,
        0.003919f,
        0.003891f,
        0.003861f,
        0.003831f,
        0.003799f,
        0.003767f,
        0.003734f,
        0.003701f,
        0.003666f,
        0.003631f,
        0.003595f,
        0.003558f,
        0.003521f,
        0.003483f,
        0.003444f,
        0.003404f,
        0.003364f,
        0.003324f,
        0.003282f,
        0.003240f,
        0.003198f,
        0.003155f,
        0.003112f,
        0.003068f,
        0.003023f,
        0.002979f,
        0.002933f,
        0.002888f,
        0.002842f,
        0.002795f,
        0.002749f,
        0.002702f,
        0.002655f,
        0.002607f,
        0.002559f,
        0.002511f,
        0.002463f,
        0.002415f,
        0.002367f,
        0.002318f,
        0.002270f,
        0.002221f,
        0.002172f,
        0.002124f,
        0.002075f,
        0.002027f,
        0.001978f,
        0.001930f,
        0.001881f,
        0.001833f,
        0.001785f,
        0.001737f,
        0.001690f,
        0.001642f,
        0.001595f,
        0.001549f,
        0.001502f,
        0.001456f,
        0.001410f,
        0.001365f,
        0.001320f,
        0.001275f,
        0.001231f,
        0.001187f,
        0.001144f,
        0.001101f,
        0.001059f,
        0.001018f,
        0.000977f,
        0.000936f,
        0.000896f,
        0.000857f,
        0.000819f,
        0.000781f,
        0.000744f,
        0.000708f,
        0.000672f,
        0.000637f,
        0.000603f,
        0.000570f,
        0.000537f,
        0.000505f,
        0.000475f,
        0.000445f,
        0.000415f,
        0.000387f,
        0.000360f,
        0.000333f,
        0.000308f,
        0.000283f,
        0.000260f,
        0.000237f,
        0.000215f,
        0.000195f,
        0.000175f,
        0.000156f,
        0.000139f,
        0.000122f,
        0.000106f,
        0.000092f,
        0.000078f,
        0.000066f,
        0.000055f,
        0.000044f,
        0.000035f,
        0.000027f,
        0.000020f,
        0.000014f,
        0.000009f,
        0.000005f,
        0.000002f,
        0.000001f,
        0.0f,
        0.0f,
        0.0f,
        0.0f,
        0.0f,
        0.0f,
        0.0f,
        0.0f,
        0.0f,
        0.0f,
        0.0f,
        0.0f,
        0.0f,
        0.0f,
        0.0f,
        0.0f,
        0.0f,
        0.0f,
        0.0f,
        0.0f,
        0.0f
    };

    private static final float[] W = {
        0.0f,
        0.0f,
        0.0f,
        0.0f,
        0.0f,
        0.0f,
        0.0f,
        0.0f,
        0.0f,
        0.0f,
        0.0f,
        0.0f,
        0.0f,
        0.0f,
        0.0f,
        0.0f,
        0.000001f,
        0.0f,
        0.0f,
        0.0f,
        0.000001f,
        0.000001f,
        0.000001f,
        0.000001f,
        0.000001f,
        0.0f,
        0.0f,
        -0.000001f,
        -0.000001f,
        -0.000001f,
        0.000001f,
        0.000001f,
        0.0f,
        0.0f,
        0.0f,
        0.0f,
        0.0f,
        0.000001f,
        0.000002f,
        0.0f,
        -0.000001f,
        -0.000001f,
        0.0f,
        0.000001f,
        0.000001f,
        0.0f,
        -0.000001f,
        -0.000003f,
        -0.000002f,
        0.0f,
        0.000002f,
        0.000002f,
        0.000002f,
        0.0f,
        -0.000003f,
        -0.000004f,
        -0.000001f,
        0.000001f,
        -0.000001f,
        -0.000001f,
        0.000002f,
        0.000004f,
        0.000003f,
        -0.000001f,
        -0.000002f,
        0.0f,
        0.0f,
        -0.000002f,
        -0.000002f,
        0.000002f,
        0.000004f,
        0.000001f,
        -0.000002f,
        0.0f,
        0.000003f,
        0.000001f,
        -0.000003f,
        -0.000002f,
        0.000001f,
        -0.000001f,
        -0.000004f,
        -0.000002f,
        0.000002f,
        0.000001f,
        -0.000004f,
        -0.000006f,
        -0.000004f,
        -0.000003f,
        -0.000003f,
        0.000001f,
        0.000003f,
        0.000001f,
        -0.000001f,
        0.000001f,
        0.000002f,
        0.000001f,
        0.0f,
        0.000002f,
        0.000003f,
        -0.000001f,
        -0.000003f,
        0.0f,
        0.000003f,
        0.000001f,
        -0.000001f,
        0.0f,
        0.000004f,
        0.000006f,
        0.000003f,
        -0.000002f,
        -0.000003f,
        0.000001f,
        0.000002f,
        -0.000002f,
        -0.000004f,
        0.000001f,
        0.000006f,
        0.000005f,
        -0.000001f,
        -0.000001f,
        0.000005f,
        0.000008f,
        0.000001f,
        -0.000006f,
        -0.000006f,
        0.000001f,
        0.000005f,
        0.000003f,
        -0.000001f,
        -0.000002f,
        0.0f,
        0.000002f,
        0.000002f,
        0.000001f,
        -0.000001f,
        -0.000003f,
        -0.000003f,
        -0.000001f,
        0.0f,
        0.0f,
        0.0f,
        0.0f,
        0.0f,
        0.000002f,
        0.000003f,
        0.000002f,
        -0.000004f,
        -0.000008f,
        -0.000005f,
        0.0f,
        0.0f,
        -0.000001f,
        0.000003f,
        0.000007f,
        0.000005f,
        -0.000002f,
        -0.000005f,
        -0.000004f,
        -0.000004f,
        -0.000006f,
        -0.000004f,
        -0.000004f,
        -0.000004f,
        -0.000001f,
        0.000005f,
        0.000005f,
        -0.000002f,
        -0.000003f,
        0.000006f,
        0.000011f,
        0.000008f,
        0.000006f,
        0.000008f,
        0.000005f,
        -0.000001f,
        0.000001f,
        0.000005f,
        0.000002f,
        -0.000002f,
        0.000002f,
        0.000006f,
        0.0f,
        -0.000008f,
        -0.000007f,
        -0.000003f,
        0.0f,
        0.000004f,
        0.000003f,
        -0.000007f,
        -0.000015f,
        -0.000008f,
        -0.000001f,
        -0.000006f,
        -0.000009f,
        -0.000001f,
        0.000002f,
        -0.000005f,
        -0.000006f,
        -0.000003f,
        -0.000008f,
        -0.000010f,
        0.000003f,
        0.000011f,
        -0.000001f,
        -0.000007f,
        0.000005f,
        0.000007f,
        -0.000008f,
        -0.000007f,
        0.000007f,
        -0.000002f,
        -0.000018f,
        -0.000003f,
        0.000015f,
        -0.000003f,
        -0.000016f,
        0.000007f,
        0.000013f,
        -0.000015f,
        -0.000010f,
        0.000024f,
        0.000009f,
        -0.000031f,
        -0.000009f,
        0.000021f,
        -0.000019f,
        -0.000036f,
        0.000030f,
        0.000038f,
        -0.000049f,
        -0.000039f,
        0.000064f,
        0.000018f,
        -0.000097f,
        -0.000003f,
        0.000115f,
        -0.000051f,
        -0.000153f,
        0.000114f,
        0.000172f,
        -0.000229f,
        -0.000184f,
        0.000399f,
        0.000137f,
        -0.000693f,
        0.000027f,
        0.001214f,
        -0.000504f,
        -0.002207f,
        0.002046f,
        0.004537f,
        -0.008338f,
        -0.012548f,
        0.063965f,
        0.261109f,
        0.495783f,
        0.602722f,
        0.495783f,
        0.261109f,
        0.063965f,
        -0.012548f,
        -0.008338f,
        0.004537f,
        0.002046f,
        -0.002207f,
        -0.000504f,
        0.001214f,
        0.000027f,
        -0.000693f,
        0.000137f,
        0.000399f,
        -0.000184f,
        -0.000229f,
        0.000172f,
        0.000114f,
        -0.000153f,
        -0.000051f,
        0.000115f,
        -0.000003f,
        -0.000097f,
        0.000018f,
        0.000064f,
        -0.000039f,
        -0.000049f,
        0.000038f,
        0.000030f,
        -0.000036f,
        -0.000019f,
        0.000021f,
        -0.000009f,
        -0.000031f,
        0.000009f,
        0.000024f,
        -0.000010f,
        -0.000015f,
        0.000013f,
        0.000007f,
        -0.000016f,
        -0.000003f,
        0.000015f,
        -0.000003f,
        -0.000018f,
        -0.000002f,
        0.000007f,
        -0.000007f,
        -0.000008f,
        0.000007f,
        0.000005f,
        -0.000007f,
        -0.000001f,
        0.000011f,
        0.000003f,
        -0.000010f,
        -0.000008f,
        -0.000003f,
        -0.000006f,
        -0.000005f,
        0.000002f,
        -0.000001f,
        -0.000009f,
        -0.000006f,
        -0.000001f,
        -0.000008f,
        -0.000015f,
        -0.000007f,
        0.000003f,
        0.000004f,
        0.0f,
        -0.000003f,
        -0.000007f,
        -0.000008f,
        0.0f,
        0.000006f,
        0.000002f,
        -0.000002f,
        0.000002f,
        0.000005f,
        0.000001f,
        -0.000001f,
        0.000005f,
        0.000008f,
        0.000006f,
        0.000008f,
        0.000011f,
        0.000006f,
        -0.000003f,
        -0.000002f,
        0.000005f,
        0.000005f,
        -0.000001f,
        -0.000004f,
        -0.000004f,
        -0.000004f,
        -0.000006f,
        -0.000004f,
        -0.000004f,
        -0.000005f,
        -0.000002f,
        0.000005f,
        0.000007f,
        0.000003f,
        -0.000001f,
        0.0f,
        0.0f,
        -0.000005f,
        -0.000008f,
        -0.000004f,
        0.000002f,
        0.000003f,
        0.000002f,
        0.0f,
        0.0f,
        0.0f,
        0.0f,
        0.0f,
        -0.000001f,
        -0.000003f,
        -0.000003f,
        -0.000001f,
        0.000001f,
        0.000002f,
        0.000002f,
        0.0f,
        -0.000002f,
        -0.000001f,
        0.000003f,
        0.000005f,
        0.000001f,
        -0.000006f,
        -0.000006f,
        0.000001f,
        0.000008f,
        0.000005f,
        -0.000001f,
        -0.000001f,
        0.000005f,
        0.000006f,
        0.000001f,
        -0.000004f,
        -0.000002f,
        0.000002f,
        0.000001f,
        -0.000003f,
        -0.000002f,
        0.000003f,
        0.000006f,
        0.000004f,
        0.0f,
        -0.000001f,
        0.000001f,
        0.000003f,
        0.0f,
        -0.000003f,
        -0.000001f,
        0.000003f,
        0.000002f,
        0.0f,
        0.000001f,
        0.000002f,
        0.000001f,
        -0.000001f,
        0.000001f,
        0.000003f,
        0.000001f,
        -0.000003f,
        -0.000003f,
        -0.000004f,
        -0.000006f,
        -0.000004f,
        0.000001f,
        0.000002f,
        -0.000002f,
        -0.000004f,
        -0.000001f,
        0.000001f,
        -0.000002f,
        -0.000003f,
        0.000001f,
        0.000003f,
        0.0f,
        -0.000002f,
        0.000001f,
        0.000004f,
        0.000002f,
        -0.000002f,
        -0.000002f,
        0.0f,
        0.0f,
        -0.000002f,
        -0.000001f,
        0.000003f,
        0.000004f,
        0.000002f,
        -0.000001f,
        -0.000001f,
        0.000001f,
        -0.000001f,
        -0.000004f,
        -0.000003f,
        0.0f,
        0.000002f,
        0.000002f,
        0.000002f,
        0.0f,
        -0.000002f,
        -0.000003f,
        -0.000001f,
        0.0f,
        0.000001f,
        0.000001f,
        0.0f,
        -0.000001f,
        -0.000001f,
        0.0f,
        0.000002f,
        0.000001f,
        0.0f,
        0.0f,
        0.0f,
        0.0f,
        0.0f,
        0.000001f,
        0.000001f,
        -0.000001f,
        -0.000001f,
        -0.000001f,
        0.0f,
        0.0f,
        0.000001f,
        0.000001f,
        0.000001f,
        0.000001f,
        0.000001f,
        0.0f,
        0.0f,
        0.0f,
        0.000001f,
        0.0f,
        0.0f,
        0.0f,
        0.0f,
        0.0f,
        0.0f,
        0.0f,
        0.0f,
        0.0f,
        0.0f,
        0.0f,
        0.0f,
        0.0f,
        0.0f,
        0.0f
    };

    public Analyze(FFT fftencode) {
        m_fft = fftencode;
        m_mem_fir = new float[NLPCOEF.length];
        m_coswintab = new float[M_PITCH / DEC];
        m_sq = new float[M_PITCH];   // these are init to 0.0
        m_Sn = new float[M_PITCH];   // ditto
        m_previous_Wo = 0.0f;

        for (int i = 0; i < M_PITCH; i++) {
            m_Sn[i] = 1.0f;
        }

        m_mem_x = 0.0f;
        m_mem_y = 0.0f;

        m_fw = new Complex[FFT_SIZE];

        for (int i = 0; i < FFT_SIZE; i++) {
            m_fw[i] = new Complex();
        }

        /*
         * Pre-compute decimation window
         */
        for (int i = 0; i < (M_PITCH / DEC); i++) {
            m_coswintab[i] = 0.5f - 0.5f * (float) Math.cos(TAU * i / (float) (M_PITCH / DEC - 1));
        }
    }

    public void reset() {
        for (int i = 0; i < M_PITCH; i++) {
            m_Sn[i] = 1.0f;
        }
    }

    /*
     * Determines the pitch in samples using a Non Linear Pitch (NLP) algorithm.
     *
     * Returns the pitch in Hz.
     */
    private float nlp() {
        /*
         * Square, notch filter at DC, and LP filter vector
         */
        for (int i = (M_PITCH - N_SAMP); i < M_PITCH; i++) {
            m_sq[i] = (m_Sn[i] * m_Sn[i]);

            float notch = m_sq[i] - m_mem_x;
            notch += (COEFF * m_mem_y);

            m_mem_x = m_sq[i];
            m_mem_y = notch;
            m_sq[i] = notch + 1.0f;

            System.arraycopy(m_mem_fir, 1, m_mem_fir, 0, (NLPCOEF.length - 1));   // shift memory left
            m_mem_fir[NLPCOEF.length - 1] = m_sq[i];

            m_sq[i] = 0.0f;

            for (int j = 0; j < NLPCOEF.length; j++) {
                m_sq[i] += (m_mem_fir[j] * NLPCOEF[j]);
            }
        }

        // decimate and DFT
        for (int i = 0; i < (M_PITCH / DEC); i++) {
            m_fw[i] = new Complex(m_sq[i * DEC] * m_coswintab[i], 0.0f);
        }

        // fill the rest with 0.0
        for (int i = (M_PITCH / DEC); i < FFT_SIZE; i++) {
            m_fw[i] = new Complex();
        }

        m_fft.transform(m_fw);

        float[] temp = new float[FFT_SIZE];

        for (int i = 0; i < FFT_SIZE; i++) {
            temp[i] = ComplexMath.csqr(m_fw[i]);
        }

        /* find global peak */
        float gmax = 0.0f;
        int gmax_bin = FFT_SIZE * DEC / P_MAX;      // 16

        /*
         * 512 * 5 / 160(P_MAX) = 16  and 512 * 5 / 20(P_MIN) = 128
         */
        for (int i = FFT_SIZE * DEC / P_MAX; i <= FFT_SIZE * DEC / P_MIN; i++) {
            if (temp[i] > gmax) {
                gmax = temp[i];
                gmax_bin = i;
            }
        }

        /* Shift samples in buffer to make room for new samples */
        System.arraycopy(m_sq, N_SAMP, m_sq, 0, M_PITCH - N_SAMP);

        int prev_fo_bin = (int) (m_previous_Wo * PREV_BIN);

        /*
         * post process estimate by searching submultiples
         */
        int mult = 2;
        int cmax_bin = gmax_bin;
        int lmax_bin;
        float thresh, lmax;

        while ((gmax_bin / mult) >= MIN_BIN) {
            int b = (gmax_bin / mult);			// determine search interval
            int bmin = (int) (0.8f * b);
            int bmax = (int) (1.2f * b);

            if (bmin < MIN_BIN) {
                bmin = MIN_BIN;
            }

            /*
             * lower threshold to favour previous frames pitch estimate,
             * this is a form of pitch tracking
             */
            if ((prev_fo_bin > bmin) && (prev_fo_bin < bmax)) {
                thresh = CNLP * 0.5f * gmax;
            } else {
                thresh = CNLP * gmax;
            }

            lmax = 0.0f;
            lmax_bin = bmin;

            for (b = bmin; b <= bmax; b++) {        // look for maximum in interval
                if (temp[b] > lmax) {
                    lmax = temp[b];
                    lmax_bin = b;
                }
            }

            if (lmax > thresh) {
                if ((lmax > temp[lmax_bin - 1]) && (lmax > temp[lmax_bin + 1])) {
                    cmax_bin = lmax_bin;
                }
            }

            mult++;
        }

        // pitch = sample rate / best Fo
        return FS / ((float) cmax_bin * FS / (float) (FFT_SIZE * DEC));
    }

    /*
     * Extract sinusoidal model parameters from input speech samples. 10
     * milliseconds of new speech samples are added during each call.
     */
    public void analyze_one_frame(Model model, short[] speech, int index) {
        Complex[] Swe = new Complex[FFT_SIZE];   // complex representation

        /*
         * Initialize the working buffer with complex zero's
         */
        for (int i = 0; i < FFT_SIZE; i++) {
            Swe[i] = new Complex();
        }

        /*
         * Prepare for new 80 samples.
         *
         * The Sn array is initialized to all 1.0 values in reset()
         */
        System.arraycopy(m_Sn, N_SAMP, m_Sn, 0, M_PITCH - N_SAMP);    // M = 320, N = 80, M-N = 240

        /*
         * Now add the new samples to the end
         */
        for (int i = 0; i < N_SAMP; i++) {                       // N = 80
            m_Sn[i + (M_PITCH - N_SAMP)] = (float) speech[i + index];
        }

        /*
         * Center analysis window on time axis. We need to arrange input
         * time domain to make FFT phases correct
         */
        // move 2nd half to start of FFT input vector (160,299 --> 0,139)
        for (int i = 0; i < NW / 2; i++) {  // NW/2 = 140
            Swe[i] = new Complex(m_Sn[i + M_PITCH / 2] * HAMMING[i + M_PITCH / 2], 0.0f);
        }

        // move 1st half to end of FFT input vector (20,159 --> 372,511)
        for (int i = 0; i < NW / 2; i++) {
            Swe[(FFT_SIZE - NW / 2 + i)] = new Complex(m_Sn[i + M_PITCH / 2 - NW / 2]
                    * HAMMING[i + M_PITCH / 2 - NW / 2], 0.0f);
        }

        m_fft.transform(Swe);       // now in frequency domain

        float Wo = TAU / nlp();

        model.setWo(Wo);
        model.setL((int) ((float) Math.PI / Wo));

        // calculate model parameters
        two_stage_pitch_refinement(model, Swe);   // refine the pitch
        estimate_amplitudes(model, Swe);
        est_voicing_mbe(model, Swe);

        m_previous_Wo = model.getWo();         // used in pitch on next pass
    }

    /*
     * Refines the current pitch estimate using the harmonic sum pitch
     * estimation technique.
     */
    private void two_stage_pitch_refinement(Model model, Complex[] Swe) {
        // Coarse refinement
        float pmax = TAU / model.getWo() + 5.0f;
        float pmin = TAU / model.getWo() - 5.0f;
        hs_pitch_refinement(model, Swe, pmin, pmax, 1.0f);    // step 1.0 course

        // Fine refinement
        pmax = TAU / model.getWo() + 1.0f;
        pmin = TAU / model.getWo() - 1.0f;
        hs_pitch_refinement(model, Swe, pmin, pmax, 0.25f);   // step 0.25 fine

        float Wo = model.getWo();

        // Limit range
        if (Wo < WO_MIN) {
            Wo = WO_MIN;
            model.setWo(Wo);
        } else if (Wo > WO_MAX) {    // changed to if-else
            Wo = WO_MAX;
            model.setWo(Wo);
        }

        model.setL((int) Math.floor(Math.PI / Wo));   // was floor()
    }

    /*
     * Harmonic sum pitch refinement function.
     *
     * pmin pitch search range minimum pmax pitch search range maximum step
     * pitch search step size model current pitch estimate in model.Wo
     *
     * model refined pitch estimate in model.Wo
     */
    private void hs_pitch_refinement(Model model, Complex[] Swe, float pmin, float pmax, float pstep) {
        float p;	// current pitch

        float Wom = model.getWo();	// Wo that maximizes E
        model.setL((int) ((float) Math.PI / Wom));	// use initial pitch est. for L

        float Em = 0.0f;	// maximum energy

        // Determine harmonic sum for a range of Wo values
        for (p = pmin; p <= pmax; p += pstep) {
            float E = 0.0f;	// energy for current pitch
            float Wo = TAU / p;	// current "test" fundamental freq.

            // Sum harmonic magnitudes
            for (int m = 1; m <= model.getL(); m++) {
                int b = (int) (m * Wo / R + 0.5f);
                E += ComplexMath.csqr(Swe[b]);
            }

            // Compare to see if this is a maximum
            if (E > Em) {
                Em = E;
                Wom = Wo;
            }
        }

        model.setWo(Wom);
    }

    /*
     * Estimates the complex amplitudes of the harmonics.
     */
    private void estimate_amplitudes(Model model, Complex[] Swe) {
        float tmp = model.getWo() / R;
        int l = model.getL();

        for (int m = 1; m <= l; m++) {
            int am = (int) ((m - 0.5f) * tmp + 0.5f);     // lower
            int bm = (int) ((m + 0.5f) * tmp + 0.5f);     // upper

            // Estimate ampltude of harmonic
            float den = 0.0f;

            for (int i = am; i < bm; i++) {
                den += ComplexMath.csqr(Swe[i]);
            }

            model.setA(m, (float) Math.sqrt(den));
        }
    }

    /*
     * Returns the error of the MBE cost function for a given F0.
     */
    private void est_voicing_mbe(Model model, Complex[] Swe) {
        float signal = .0001f;
        int l = model.getL();

        for (int i = 1; i <= (l / 4); i++) {
            float aval = model.getA(i);
            signal += (aval * aval);
        }

        float Wo = model.getWo();
        float error = .0001f;

        /* Just test across the harmonics in the first 1000 Hz (L/4) */
        for (int i = 1; i <= (l / 4); i++) {
            Complex Am = new Complex();

            int al = (int) Math.ceil((i - 0.5) * Wo * FFT_SIZE / TAU);
            int bl = (int) Math.ceil((i + 0.5) * Wo * FFT_SIZE / TAU);

            int offset = (int) (FFT_SIZE / 2.0f - i * Wo * FFT_SIZE / TAU + 0.5f);

            for (int m = al; m < bl; m++) {
                Am = ComplexMath.add(Am, ComplexMath.times(Swe[m], W[offset + m]));
            }

            /*
             * Determine error between estimated harmonic and original
             */
            for (int m = al; m < bl; m++) {
                Complex Ew = ComplexMath.minus(Swe[m], ComplexMath.times(Am, W[offset + m]));
                error += ComplexMath.csqr(Ew);
            }
        }

        if (10.0f * (float) Math.log10(signal / error) > VOICING_THRESHOLD_DB) {
            model.setVoiced(true);
        } else {
            model.setVoiced(false);
        }

        /*
         * post processing, helps clean up some voicing errors
         *
         * Determine the ratio of low freqency to high frequency energy,
         * voiced speech tends to be dominated by low frequency energy,
         * unvoiced by high frequency. This measure can be used to
         * determine if we have made any gross errors.
         */
        float elow = .0001f;
        float ehigh = .0001f;

        for (int i = 1; i <= (l / 2); i++) {
            float aval = model.getA(i);
            elow += (aval * aval);
        }

        for (int i = (l / 2); i <= l; i++) {
            float aval = model.getA(i);
            ehigh += (aval * aval);
        }

        float eratio_dB = 10.0f * (float) Math.log10(elow / ehigh);

        if (model.getVoiced() == false) {
            /*
             * Look for Type 1 errors, strongly Voiced speech that has been
             * accidentally declared Un-Voiced
             */

            if (eratio_dB > 10.0f) {
                model.setVoiced(true);
            }
        } else if (model.getVoiced() == true) {
            /*
             * Look for Type 2 errors, strongly Un-Voiced speech that has been
             * accidentally declared Voiced
             */

            if (eratio_dB < -10.0f) {
                model.setVoiced(false);
            }

            /*
             * A common source of Type 2 errors is the pitch estimator
             * gives a low (50Hz) estimate for Un-Voiced speech, which gives a
             * good match with noise due to the close harmoonic spacing.
             * These errors are much more common than people with 50Hz
             * pitch, so we have just a small eratio_dB threshold.
             */
            if ((eratio_dB < -4.0f) && (model.getWo() <= ((60.0f * TAU) / FS))) {
                model.setVoiced(false);
            }
        }
    }
}
