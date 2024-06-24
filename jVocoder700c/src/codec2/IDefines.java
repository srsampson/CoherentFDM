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

public interface IDefines {
    static final float TAU = (float) (Math.PI * 2.0);
    static final float FS = 8000.0f;
    //
    static final int P_MIN = 20;
    static final int P_MAX = 160;
    //
    static final int FFT_SIZE = 512;
    static final int AMP_PHASE_NFFT = 128;
    static final int N_SAMP = 80;
    static final int MAX_AMP = 80;
    //
    static final int AMP_K = 20;
    static final int AMP_M = 512;
    static final int AMP_LOG2 = 9;
    //
    static final int[] FFT_FACTORS = {4, 2, 3, 5};
    //
    static final byte[] IDEAL = {
        8, 10, 12, 14, 14, 14, 14,
        14, 14, 14, 14, 14, 14, 14,
        14, 14, 14, 14, 14, -20
    };
    //
    static final float[] CODES0 = {
        10.0f, 12.5f, 15.0f, 17.5f, 20.0f,
        22.5f, 25.0f, 27.5f, 30.0f, 32.5f,
        35.0f, 37.5f, 40.0f, 42.5f, 45.0f,
        47.5f
    };
    //
    static final float[] PARZEN = {
        0.0f,
        0.0125f,
        0.025f,
        0.0375f,
        0.05f,
        0.0625f,
        0.075f,
        0.0875f,
        0.1f,
        0.1125f,
        0.125f,
        0.1375f,
        0.15f,
        0.1625f,
        0.175f,
        0.1875f,
        0.2f,
        0.2125f,
        0.225f,
        0.2375f,
        0.25f,
        0.2625f,
        0.275f,
        0.2875f,
        0.3f,
        0.3125f,
        0.325f,
        0.3375f,
        0.35f,
        0.3625f,
        0.375f,
        0.3875f,
        0.4f,
        0.4125f,
        0.425f,
        0.4375f,
        0.45f,
        0.4625f,
        0.475f,
        0.4875f,
        0.5f,
        0.5125f,
        0.525f,
        0.5375f,
        0.55f,
        0.5625f,
        0.575f,
        0.5875f,
        0.6f,
        0.6125f,
        0.625f,
        0.6375f,
        0.65f,
        0.6625f,
        0.675f,
        0.6875f,
        0.7f,
        0.7125f,
        0.725f,
        0.7375f,
        0.75f,
        0.7625f,
        0.775f,
        0.7875f,
        0.8f,
        0.8125f,
        0.825f,
        0.8375f,
        0.85f,
        0.8625f,
        0.875f,
        0.8875f,
        0.9f,
        0.9125f,
        0.925f,
        0.9375f,
        0.95f,
        0.9625f,
        0.975f,
        0.9875f,
        1.0f,
        0.9875f,
        0.975f,
        0.9625f,
        0.95f,
        0.9375f,
        0.925f,
        0.9125f,
        0.9f,
        0.8875f,
        0.875f,
        0.8625f,
        0.85f,
        0.8375f,
        0.825f,
        0.8125f,
        0.8f,
        0.7875f,
        0.775f,
        0.7625f,
        0.75f,
        0.7375f,
        0.725f,
        0.7125f,
        0.7f,
        0.6875f,
        0.675f,
        0.6625f,
        0.65f,
        0.6375f,
        0.625f,
        0.6125f,
        0.6f,
        0.5875f,
        0.575f,
        0.5625f,
        0.55f,
        0.5375f,
        0.525f,
        0.5125f,
        0.5f,
        0.4875f,
        0.475f,
        0.4625f,
        0.45f,
        0.4375f,
        0.425f,
        0.4125f,
        0.4f,
        0.3875f,
        0.375f,
        0.3625f,
        0.35f,
        0.3375f,
        0.325f,
        0.3125f,
        0.3f,
        0.2875f,
        0.275f,
        0.2625f,
        0.25f,
        0.2375f,
        0.225f,
        0.2125f,
        0.2f,
        0.1875f,
        0.175f,
        0.1625f,
        0.15f,
        0.1375f,
        0.125f,
        0.1125f,
        0.1f,
        0.0875f,
        0.075f,
        0.0625f,
        0.05f,
        0.0375f,
        0.025f,
        0.0125f
    };
    /*
     * 48 tap 600Hz low pass FIR filter coefficients
     */
    static final float[] NLPCOEF = {
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
    //
    static final float[] HAMMING = {
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
    //
    static final float[] W = {
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
}
