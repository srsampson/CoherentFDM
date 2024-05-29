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
import java.util.Random;

public final class Model implements IDefines {

    private final static float BG_THRESH_DB = 40.0f;    // only consider low levels signals for bg_est
    private final static float BG_BETA = 0.1f;          // averaging filter constant
    private final static float BG_MARGIN_DB = 6.0f;
    private final static int MAX_HARMONIC = 80;         // maximum number of harmonics
    //
    private final float[] m_A;                           // amplitude of each harmonic
    private final float[] m_phi;                         // phase of each harmonic
    private float m_Wo;                                  // fundamental frequency estimate in rad/s
    private int m_L;                                      // number of harmonics
    private boolean m_voiced;                             // if this frame is voiced
    //
    private final Random m_rand;

    /*
     * Class to hold model parameters for one 10 ms frame
     */
    public Model() {
        m_rand = new Random(System.currentTimeMillis());     // seed the noise generator
        m_phi = new float[MAX_HARMONIC+1];  // 0..80
        m_A = new float[MAX_HARMONIC+1];
        m_Wo = 0.0f;
        m_L = 3;          // (int) Math.PI
        m_voiced = false;
    }

    public void reset() {
        for (int i = 0; i <= MAX_HARMONIC; i++) {
            m_A[i] = 0.0f;
            m_phi[i] = 0.0f;
        }
    }

    public boolean getVoiced() {
        return m_voiced;
    }

    public void setVoiced(boolean val) {
        m_voiced = val;
    }

    public float getA(int index) {
        return m_A[index];
    }

    public void setA(int index, float val) {
        m_A[index] = val;
    }

    public float getPhi(int index) {
        return m_phi[index];
    }

    public void setPhi(int index, float val) {
        m_phi[index] = val;
    }

    public int getL() {
        return m_L;
    }

    public void setL(int val) {
        m_L = val;
    }

    public float getWo() {
        return m_Wo;
    }

    public void setWo(float val) {
        m_Wo = val;
    }

    /*
     * Postfilter to improve sound quality for speech with high levels of
     * background noise. Unlike mixed-excitation models requires no bits to be
     * transmitted to handle background noise.
     */
    public float postfilter(float bg_est_val) {
        float bg_est = bg_est_val;
        float e = 1E-12f;

        for (int m = 1; m <= m_L; m++) {
            e += (m_A[m] * m_A[m]);
        }

        e = 10.0f * (float) Math.log10(e / m_L);

        if ((e < BG_THRESH_DB) && !m_voiced) {
            bg_est *= (1.0f - BG_BETA) + (e * BG_BETA); // IIR
        }

        float thresh = (float) Math.pow(10.0f, (bg_est + BG_MARGIN_DB) / 20.0f);

        if (m_voiced == true) {
            for (int m = 1; m <= m_L; m++) {
                if (m_A[m] < thresh) {
                    m_phi[m] = TAU * m_rand.nextFloat(); // random value 0.0..TAU
                }
            }
        }

        return bg_est;
    }

    /*
     * Synthesizes phases based on SNR and a rule based approach. No phase
     * parameters are required apart from the SNR (which can be reduced to a 1
     * bit V/UV decision per frame).
     */
    public float phase_synth_zero_order(float ex_phase_val, Complex[] Aw) {
        float ex_phase = ex_phase_val;
        Complex Ex;

        /*
         * Update excitation fundamental phase track, this sets the position
         * of each pitch pulse during voiced speech.
         */
        ex_phase += (m_Wo * N_SAMP);
        ex_phase -= TAU * (float) Math.floor(ex_phase / TAU + 0.5f);

        for (int m = 1; m <= m_L; m++) {

            /*
             * generate excitation
             */
            if (m_voiced == true) {
                Ex = ComplexMath.cexp(new Complex(0.0f, ex_phase * (float) m));
            } else {
                /*
                 * When a few samples were tested I found that LPC filter
                 * phase is not needed in the un-voiced case, but no harm in
                 * keeping it.
                 */

                Ex = ComplexMath.cexp(new Complex(0.0f, TAU * m_rand.nextFloat()));
            }

            /*
             * filter using LPC filter
             * modify sinusoidal phase
             */
            m_phi[m] = ComplexMath.carg(ComplexMath.times(Aw[m], Ex));
        }

        return ex_phase;
    }
}
