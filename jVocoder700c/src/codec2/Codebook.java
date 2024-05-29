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

public final class Codebook {

    private final int m_k;         // dimension of the vector
    private final int m_log2m;     // log base 2 of size m
    private final int m_m;         // number of vector elements
    private final float[] m_cb;    // the codebook array reference pointer

    public Codebook(int kval, int log2mval, int mval, float[] cbval) {
        m_k = kval;
        m_log2m = log2mval;
        m_m = mval;
        m_cb = cbval;
    }

    public int getDimension() {
        return m_k;
    }

    public int getLogBase2Elements() {
        return m_log2m;
    }

    public int getNumberOfElements() {
        return m_m;
    }

    public float[] getCodeBookArray() {
        return m_cb.clone();
    }

    public float getCodeBookArray(int index) {
        return m_cb[index];
    }
}
