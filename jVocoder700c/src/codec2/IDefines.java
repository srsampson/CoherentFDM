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
}
