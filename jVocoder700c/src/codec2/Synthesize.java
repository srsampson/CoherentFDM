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
    private float m_bg_est;      // background noise estimate for post filter
    private float m_ex_phase;    // phase tracking

    protected Synthesize(FFT fftdecode) {
        m_fft = fftdecode;
        m_Swd = new Complex[FFT_SIZE];
        m_Sn_ = new float[N_SAMP * 2];
    }

    protected void reset() {
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
    protected void synthesize_one_frame(Model model, short[] speech, int index, Complex[] Aw) {
        m_ex_phase = model.phase_synth_zero_order(m_ex_phase, Aw);
        m_bg_est = model.postfilter(m_bg_est);
        synthesize(model);                  // get updated Sn_

        /*
         * find maximum sample in frame for ear protection
         */
        float max_sample = m_Sn_[0];

        for (int i = 1; i < N_SAMP; i++) {
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

        // Perform inverse DFT no scale
        m_fft.itransform(m_Swd, false);

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
