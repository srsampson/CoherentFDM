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

public final class Synthesize implements IDefines {

    private final FFT m_fft;
    //
    private final Complex[] m_Swd;        // the decode speech spectrum
    //
    private final float[] m_Sn_;
    //
    private static final float[] PARZEN = {
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

    private float m_bg_est;      // background noise estimate for post filter
    private float m_ex_phase;    // phase tracking

    public Synthesize(FFT fftdecode) {
        m_fft = fftdecode;
        m_Swd = new Complex[FFT_SIZE];
        m_Sn_ = new float[N_SAMP * 2];
    }

    public void reset() {
        m_bg_est = 0.0f;
        m_ex_phase = 0.0f;
    }

    /*
     * Synthesize 80 speech samples (10ms) from model parameters. Limits output
     * level to protect ears when there are bit errors or the input is over
     * driven.
     *
     * This doesn't correct or mask bit errors, just reduces the effects.
     */
    public void synthesize_one_frame(Model model, short[] speech, int index, Complex[] Aw) {
        m_ex_phase = model.phase_synth_zero_order(m_ex_phase, Aw);
        m_bg_est = model.postfilter(m_bg_est);
        synthesize(model);                  // get updated Sn_

        /*
         * find maximum sample in frame for ear protection
         */
        float max_sample = 0.0f;

        for (int i = 0; i < N_SAMP; i++) {
            if (m_Sn_[i] > max_sample) {
                max_sample = m_Sn_[i];
            }
        }

        /*
         * determine how far above set point
         */
        float over = max_sample / 30000.0f;

        /*
         * If we are x dB over set point we reduce level by 2x dB, this
         * attenuates major excursions in amplitude (likely to be caused
         * by bit errors) more than smaller ones
         */
        if (over > 1.0f) {
            float gain = (over * over);

            for (int i = 0; i < N_SAMP; i++) {
                m_Sn_[i] /= gain;
            }
        }

        for (int i = 0; i < N_SAMP; i++) {           // 80
            if (m_Sn_[i] > 32767.0f) {
                speech[i + index] = 32760;
            } else if (m_Sn_[i] < -32767.0f) {
                speech[i + index] = -32760;
            } else {
                speech[i + index] = (short) m_Sn_[i];
            }
        }
    }

    /*
     * Synthesize a speech signal in the frequency domain from the sinusoidal
     * model parameters. Uses overlap-add with a trapezoidal window to smoothly
     * interpolate between frames.
     */
    private void synthesize(Model model) {
        /*
         * Initialize the working buffer to complex zero
         */
        for (int i = 0; i < FFT_SIZE; i++) {
            m_Swd[i] = new Complex();
        }

        // Update window by shifting 80 samples left (10 ms)
        System.arraycopy(m_Sn_, N_SAMP, m_Sn_, 0, N_SAMP - 1);
        m_Sn_[N_SAMP - 1] = 0.0f;

        /*
         * Now set up frequency domain synthesized speech
         */
        float tmp = model.getWo() * FFT_SIZE / TAU;

        for (int i = 1; i <= model.getL(); i++) {
            int b = (int) (i * tmp + 0.5);

            if (b > ((FFT_SIZE / 2) - 1)) {
                b = (FFT_SIZE / 2) - 1;
            }

            m_Swd[b] = ComplexMath.times(ComplexMath.cexp(new Complex(0.0f, model.getPhi(i))), model.getA(i));
            m_Swd[FFT_SIZE - b] = ComplexMath.conjugate(m_Swd[b]);
        }

        // Perform inverse DFT
        m_fft.itransform(m_Swd);

        // Overlap add to previous samples
        for (int i = 0; i < (N_SAMP - 1); i++) {
            m_Sn_[i] += (m_Swd[FFT_SIZE - N_SAMP + 1 + i].real() * PARZEN[i]);
        }

        // put the new data on the end of the window
        for (int i = N_SAMP - 1, j = 0; i < (N_SAMP * 2); i++, j++) {
            m_Sn_[i] = (m_Swd[j].real() * PARZEN[i]);
        }
    }
}
