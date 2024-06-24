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
import java.util.ArrayList;

public final class Amp implements IDefines {

    private static final int NS = (AMP_PHASE_NFFT / 2 + 1);
    private static final int MBEST_ENTRIES = 5;
    private static final int MBEST_STAGES = 4;
    //
    private static final float WO_MIN = (float) Math.log10((TAU / P_MAX));
    private static final float WO_MAX = (float) Math.log10((TAU / P_MIN));
    //
    private static final float MEL200 = (float) Math.floor(2595.0 * Math.log10(1.0 + 200 / 700.0) + 0.5);
    private static final float MEL3700 = (float) Math.floor(2595.0 * Math.log10(1.0 + 3700 / 700.0) + 0.5);
    private static final float MELSTEP = (MEL3700 - MEL200);
    private static final float PHASE_SCALE = (float) (20.0 / Math.log(10.0));
    //
    private final CodebookVQ m_codebookVQ;
    private final FFT fft;
    //
    private final float[] m_rate_K_sample_freqs_kHz;
    private final float[] m_rate_K_vec;
    private final float[] m_prev_rate_K_vec_;
    //
    private final float[][] m_interpolated_surface_;
    //
    private float[] m_eq;   // equalization
    private float m_se;     // squared error
    private long m_nse;      // counter
    //
    private float m_Wo_left;
    private boolean m_voicing_left;

    protected class MBest {

        private final int[] mm_index;
        private float mm_error;

        protected MBest() {
            mm_index = new int[MBEST_STAGES];
            mm_error = 1E32f;
        }

        protected void reset() {
            for (int i = 0; i < MBEST_STAGES; i++) {
                mm_index[i] = 0;
            }

            mm_error = 1E32f;
        }

        protected void resetIndex() {
            for (int i = 0; i < MBEST_STAGES; i++) {
                mm_index[i] = i;
            }
        }

        protected void setIndex(int ind, int val) {
            mm_index[ind] = val;
        }

        protected int getIndex(int ind) {
            return mm_index[ind];
        }

        protected void setError(float val) {
            mm_error = val;
        }

        protected float getError() {
            return mm_error;
        }
    }

    protected Amp(FFT fftphase) {
        fft = fftphase;
        m_interpolated_surface_ = new float[4][AMP_K];
        m_rate_K_sample_freqs_kHz = new float[AMP_K];
        m_rate_K_vec = new float[AMP_K];
        m_prev_rate_K_vec_ = new float[AMP_K];
        m_eq = new float[AMP_K];
        m_se = 0.0f;
        m_nse = 0L;

        float step = MELSTEP / (float) (AMP_K - 1);
        float mel = MEL200;

        m_codebookVQ = new CodebookVQ();

        for (int k = 0; k < AMP_K; k++) {
            m_rate_K_sample_freqs_kHz[k] = 0.7f * (float) (Math.pow(10.0, (mel / 2595.0f)) - 1.0f);
            mel += step;
        }
    }

    protected int getAMPK() {
        return AMP_K;
    }

    /**
     * Method to reset the Equalizer values to zero.
     */
    protected void resetEQ() {
        m_eq = new float[AMP_K]; // init to 0.0
        m_se = 0.0f;
        m_nse = 0L;
    }

    /**
     * Method to return a copy of the EQ array at this instant.
     * Used for debugging purposes or graphic display.
     *
     * @return float array containing the Equalizer values.
     */
    protected float[] getEQ() {
        return m_eq.clone();
    }

    protected float getSquaredError() {
        return m_se;
    }

    protected long getSquaredCount() {
        return m_nse;
    }

    protected void amp_indexes_to_model(Model[] model_, Complex[][] HH, int[] indexes) {
        float[] rate_K_vec_ = new float[AMP_K];
        float[] left_vec;
        float[] right_vec;
        float[] aWo_ = new float[4];
        boolean[] avoicing_ = new boolean[4];
        int[] aL_ = new int[4];
        float Wo_right;
        boolean voicing_right;

        // extract latest rate K vector
        amp_indexes_to_rate_K_vec(rate_K_vec_, indexes);

        // decode latest Wo and voicing
        if (indexes[3] != 0) {
            Wo_right = decode_log_Wo(indexes[3], 6);
            voicing_right = true;
        } else {
            Wo_right = TAU / 100.0f;
            voicing_right = false;
        }

        left_vec = m_prev_rate_K_vec_;
        right_vec = rate_K_vec_;

        // interpolate 25Hz rate K vec back to 100Hz
        amp_interpolate(left_vec, right_vec);

        /* interpolate 25Hz v and Wo back to 100Hz */
        interp_Wo_v(aWo_, aL_, avoicing_, m_Wo_left, Wo_right, m_voicing_left, voicing_right);

        // back to rate L amplitudes, synthesis phase for each frame
        for (int i = 0; i < 4; i++) {
            model_[i].setWo(aWo_[i]);
            model_[i].setL(aL_[i]);
            model_[i].setVoiced(avoicing_[i]);

            resample_rate_L(model_[i], m_interpolated_surface_, i);
            determine_phase(HH, i, model_[i]);
        }

        // update memories for next time
        System.arraycopy(rate_K_vec_, 0, m_prev_rate_K_vec_, 0, AMP_K);

        m_Wo_left = Wo_right;
        m_voicing_left = voicing_right;
    }

    private void amp_eq(float[] rate_K_vec_no_mean, boolean eq_en) {
        float gain = 0.02f;

        for (int k = 0; k < AMP_K; k++) {
            float update = rate_K_vec_no_mean[k] - (float) IDEAL[k];

            // Always update eq array independant of eq_en
            m_eq[k] = (1.0f - gain) * m_eq[k] + gain * update;

            if (m_eq[k] < 0.0f) {
                m_eq[k] = 0.0f;
            }

            // Only update no_mean array if eq_en is true
            if (eq_en == true) {
                rate_K_vec_no_mean[k] -= m_eq[k];
            }
        }
    }

    protected void amp_model_to_indexes(Model model, int[] indexes, boolean eq_en) {
        float[] rate_K_vec_no_mean = new float[AMP_K];
        float[] rate_K_vec_no_mean_ = new float[AMP_K];

        // convert variable rate L to fixed rate K
        resample_const_rate_f(model, m_rate_K_vec);

        // remove mean and two stage VQ
        float sum = 0.0f;

        for (int i = 0; i < AMP_K; i++) {
            sum += m_rate_K_vec[i];
        }

        float mean = sum / AMP_K;

        for (int i = 0; i < AMP_K; i++) {
            rate_K_vec_no_mean[i] = m_rate_K_vec[i] - mean;
        }

        // update eq and optionally run equaliser on no_mean array
        amp_eq(rate_K_vec_no_mean, eq_en);

        rate_K_mbest_encode(indexes, rate_K_vec_no_mean, rate_K_vec_no_mean_);

        // running sum of squared error for variance calculation
        for (int i = 0; i < AMP_K; i++) {
            float tmp = (rate_K_vec_no_mean[i] - rate_K_vec_no_mean_[i]);
            m_se += (tmp * tmp);
        }

        m_nse += AMP_K;

        // scalar quantise the mean (effectively the frame energy)
        indexes[2] = quantizeAmplitude(mean);

        // scalar quantise Wo.  We steal the smallest Wo index
        // to signal an unvoiced frame
        if (model.getVoiced()) {
            int index = encode_log_Wo(model.getWo(), 6);

            if (index == 0) {
                index = 1;
            }

            indexes[3] = index;
        } else {
            indexes[3] = 0;
        }
    }

    /*
     * Insert the results in a vector for codebook entry comparison. The
     * list is ordered by error, so those entries with the smallest error
     * will be first on the list.
     */
    private void mbest_insert(MBest[] mbest, int[] index, float error) {
        for (int i = 0; i < MBEST_ENTRIES; i++) {
            if (error < mbest[i].getError()) {
                for (int j = MBEST_ENTRIES - 1; j > i; j--) {
                    mbest[j] = mbest[j - 1];
                }

                for (int j = 0; j < MBEST_STAGES; j++) {
                    mbest[i].setIndex(j, index[j]);
                }

                mbest[i].setError(error);
                break;  // break out of for loop, we're done
            }
        }
    }

    /*
     * Searches vec[] to a codebook of vectors, and maintains a list of the mbest
     * closest matches.
     */
    private void mbest_search(ArrayList<Float> cb, float[] vec, MBest[] mbest, int[] index) {
        for (int j = 0; j < AMP_M; j++) {
            float e = 0.0f;

            for (int i = 0; (i < AMP_K) && (e < mbest[MBEST_ENTRIES - 1].getError()); i++) {
                float diff = cb.get(j * AMP_K + i) - vec[i];
                e += (diff * diff);
            }

            index[0] = j;
            mbest_insert(mbest, index, e);
        }
    }

    private void mag_to_phase(float[] phase, float[] Gdbfk) {
        Complex[] Sdb = new Complex[AMP_PHASE_NFFT];
        Complex[] cf = new Complex[AMP_PHASE_NFFT];

        // initialize values to zero
        for (int i = 0; i < AMP_PHASE_NFFT; i++) {
            Sdb[i] = new Complex();
            cf[i] = new Complex();
        }

        // install negative frequency components
        Sdb[0] = new Complex(Gdbfk[0], 0.0f);

        for (int i = 1; i < NS; i++) {
            Sdb[i] = new Complex(Gdbfk[i], 0.0f);
            Sdb[AMP_PHASE_NFFT - i] = new Complex(Gdbfk[i], 0.0f);
        }

        // compute real cepstrum from log magnitude spectrum
        fft.itransform(Sdb, true);

        // Fold cepstrum to reflect non-min-phase zeros inside unit circle
        cf[0] = Sdb[0];

        for (int i = 1; i < NS - 1; i++) {
            cf[i] = ComplexMath.add(Sdb[i], Sdb[AMP_PHASE_NFFT - i]);
        }

        cf[NS - 1] = Sdb[NS - 1];

        // cf = dB_magnitude + j * minimum_phase
        fft.transform(cf);

        /*
         * The maths says we are meant to be using log(x), not 20*log10(x), so
         * we need to scale the phase to account for this: log(x) =
         * 20*log10(x)/scale
         */
        for (int i = 0; i < NS; i++) {
            phase[i] = cf[i].imag() / PHASE_SCALE;
        }
    }

    /*
     * We add some variables here to translate the C pointers into indexes.
     */
    private void interp_para(float[] y, int yindex, float[] xp, int xpindex,
            float[] yp, int ypindex, int np, float[] x, int xindex, int n) {
        float xi, x1, y1, x2, y2, x3, y3, a, b;
        int k = 0;

        for (int i = 0; i < n; i++) {
            xi = x[i + xindex];

            /* k is index into xp of where we start 3 points used to form parabola */
            while ((xp[k + xpindex + 1] < xi) && (k < (np - 3))) {
                k++;
            }

            x1 = xp[k + xpindex];
            y1 = yp[k + ypindex];

            x2 = xp[k + xpindex + 1];
            y2 = yp[k + ypindex + 1];

            x3 = xp[k + xpindex + 2];
            y3 = yp[k + ypindex + 2];

            a = ((y3 - y2) / (x3 - x2) - (y2 - y1) / (x2 - x1)) / (x3 - x1);
            b = ((y3 - y2) / (x3 - x2) * (x2 - x1) + (y2 - y1) / (x2 - x1) * (x3 - x2)) / (x3 - x1);

            y[i + yindex] = a * (xi - x2) * (xi - x2) + b * (xi - x2) + y2;
        }
    }

    private void resample_const_rate_f(Model model, float[] rate_K_vec) {
        float[] AmdB = new float[MAX_AMP + 1];
        float[] rate_L_sample_freqs_kHz = new float[MAX_AMP + 1];
        float AmdB_peak;

        /* convert rate L=pi/Wo amplitude samples to fixed rate K */
        AmdB_peak = -100.0f;

        for (int m = 1; m <= model.getL(); m++) {
            AmdB[m] = 20.0f * (float) Math.log10(model.getA(m) + 1E-16);

            if (AmdB[m] > AmdB_peak) {
                AmdB_peak = AmdB[m];
            }

            rate_L_sample_freqs_kHz[m] = m * model.getWo() * 4.0f / (float) Math.PI;
        }

        /* clip between peak and peak -50dB, to reduce dynamic range */
        for (int m = 1; m <= model.getL(); m++) {
            if (AmdB[m] < (AmdB_peak - 50.0f)) {
                AmdB[m] = AmdB_peak - 50.0f;
            }
        }

        interp_para(rate_K_vec, 0, rate_L_sample_freqs_kHz, 1, AmdB, 1, model.getL(),
                m_rate_K_sample_freqs_kHz, 0, AMP_K);
    }

    private void rate_K_mbest_encode(int[] indexes, float[] x, float[] xq) {
        ArrayList<Float> codebook1 = m_codebookVQ.getCodebook(0).CodeBookArray();
        ArrayList<Float> codebook2 = m_codebookVQ.getCodebook(1).CodeBookArray();

        MBest[] mbest_stage1 = new MBest[MBEST_ENTRIES];
        MBest[] mbest_stage2 = new MBest[MBEST_ENTRIES];

        for (int i = 0; i < MBEST_ENTRIES; i++) {
            mbest_stage1[i] = new MBest();
            mbest_stage2[i] = new MBest();
        }
        
        int[] index = new int[MBEST_STAGES];
        
        /* Stage 1 */
        mbest_search(codebook1, x, mbest_stage1, index);

        /* Stage 2 */
        for (int j = 0; j < MBEST_ENTRIES; j++) {
            float[] target = new float[AMP_K];
            
            index[1] = mbest_stage1[j].getIndex(0);

            for (int i = 0; i < AMP_K; i++) {
                target[i] = x[i] - codebook1.get(AMP_K * index[1] + i);
            }

            mbest_search(codebook2, target, mbest_stage2, index);
        }

        int n1 = mbest_stage2[0].getIndex(1);
        int n2 = mbest_stage2[0].getIndex(0);

        for (int i = 0; i < AMP_K; i++) {
            xq[i] = (codebook1.get(AMP_K * n1 + i) + codebook2.get(AMP_K * n2 + i));
        }

        indexes[0] = n1;
        indexes[1] = n2;
    }

    private void post_filter_amp(float[] vec, float pf_gain) {
        float[] pre = new float[AMP_K];

        /*
         * vec is rate K vector describing spectrum of current frame lets pre-emp
         * before applying PF. 20dB/dec over 300Hz. Postfilter affects energy of
         * frame so we measure energy before and after and normalise. Plenty of room
         * for experiment here as well.
         */
        float e_before = 0.0f;
        float e_after = 0.0f;

        for (int i = 0; i < AMP_K; i++) {
            pre[i] = 20.0f * (float) Math.log10(m_rate_K_sample_freqs_kHz[i] / 0.3);
            vec[i] += pre[i];
            e_before += (float) Math.pow(10.0, 2.0 * vec[i] / 20.0);
            vec[i] *= pf_gain;
            e_after += (float) Math.pow(10.0, 2.0 * vec[i] / 20.0);
        }

        float gaindB = 10.0f * (float) Math.log10(e_after / e_before);

        for (int i = 0; i < AMP_K; i++) {
            vec[i] -= gaindB;
            vec[i] -= pre[i];
        }
    }

    private void interp_Wo_v(float[] Wo_, int[] L_, boolean[] voicing_, float Wo1,
            float Wo2, boolean voicing1, boolean voicing2) {
        float omega = TAU / 100.0f;

        // interpolation rate
        for (int i = 0; i < 4; i++) {
            voicing_[i] = false;
        }

        if (!voicing1 && !voicing2) {
            for (int i = 0; i < 4; i++) {
                Wo_[i] = omega;
            }
        }

        if (voicing1 && !voicing2) {
            Wo_[0] = Wo1;
            Wo_[1] = Wo1;
            Wo_[2] = omega;
            Wo_[3] = omega;
            voicing_[0] = true;
            voicing_[1] = true;
        }

        if (!voicing1 && voicing2) {
            Wo_[0] = omega;
            Wo_[1] = omega;
            Wo_[2] = Wo2;
            Wo_[3] = Wo2;
            voicing_[2] = true;
            voicing_[3] = true;
        }

        if (voicing1 && voicing2) {
            float c;
            int j;

            for (j = 0, c = 1.0f; j < 4; j++, c -= .250f) {
                Wo_[j] = Wo1 * c + Wo2 * (1.0f - c);
                voicing_[j] = true;
            }
        }

        for (int i = 0; i < 4; i++) {
            L_[i] = (int) Math.floor(Math.PI / Wo_[i]);
        }
    }

    private void resample_rate_L(Model model, float[][] rate_K_vec, int index) {
        float[] rate_K_vec_term = new float[AMP_K + 2];
        float[] rate_K_sample_freqs_kHz_term = new float[AMP_K + 2];
        float[] AmdB = new float[MAX_AMP + 1];
        float[] rate_L_sample_freqs_kHz = new float[MAX_AMP + 1];

        /* terminate either end of the rate K vecs with 0dB points */
        rate_K_vec_term[0] = 0.0f;
        rate_K_vec_term[AMP_K + 1] = 0.0f;
        rate_K_sample_freqs_kHz_term[0] = 0.0f;
        rate_K_sample_freqs_kHz_term[AMP_K + 1] = 4.0f;

        for (int i = 0; i < AMP_K; i++) {
            rate_K_vec_term[i + 1] = rate_K_vec[index][i];
            rate_K_sample_freqs_kHz_term[i + 1] = m_rate_K_sample_freqs_kHz[i];
        }

        for (int i = 1; i <= model.getL(); i++) {
            rate_L_sample_freqs_kHz[i] = i * model.getWo() * 4.0f / (float) Math.PI;
        }

        interp_para(AmdB, 1, rate_K_sample_freqs_kHz_term, 0, rate_K_vec_term, 0,
                AMP_K + 2, rate_L_sample_freqs_kHz, 1, model.getL());

        for (int i = 1; i <= model.getL(); i++) {
            model.setA(i, (float) Math.pow(10.0, AmdB[i] / 20.0));
        }
    }

    private void determine_phase(Complex[][] H, int index, Model model) {
        float[] Gdbfk = new float[NS];
        float[] sample_freqs_kHz = new float[NS];
        float[] phase = new float[NS];
        float[] AmdB = new float[MAX_AMP + 1];
        float[] rate_L_sample_freqs_kHz = new float[MAX_AMP + 1];

        for (int i = 1; i <= model.getL(); i++) {
            AmdB[i] = 20.0f * (float) Math.log10(model.getA(i));
            rate_L_sample_freqs_kHz[i] = (float) i * model.getWo() * 4.0f / (float) Math.PI;
        }

        for (int i = 0; i < NS; i++) {
            sample_freqs_kHz[i] = (FS / 1000.0f) * (float) i / (float) AMP_PHASE_NFFT;
        }

        interp_para(Gdbfk, 0, rate_L_sample_freqs_kHz, 1, AmdB, 1, model.getL(), sample_freqs_kHz, 0, NS);
        mag_to_phase(phase, Gdbfk);

        for (int i = 1; i <= model.getL(); i++) {
            int temp = (int) Math.floor(0.5 + (float) i * model.getWo() * (float) AMP_PHASE_NFFT / TAU);
            H[index][i] = ComplexMath.cexp(new Complex(0.0f, phase[temp]));
        }
    }

    private void amp_interpolate(float[] left_vec, float[] right_vec) {
        float c;
        int i;

        // (linearly) interpolate 25Hz amplitude vectors back to 100Hz
        for (i = 0, c = 1.0f; i < 4; i++, c -= .250f) {
            for (int k = 0; k < AMP_K; k++) {
                m_interpolated_surface_[i][k] = left_vec[k] * c + right_vec[k] * (1.0f - c);
            }
        }
    }

    private void amp_indexes_to_rate_K_vec(float[] rate_K_vec_, int[] indexes) {
        float[] rate_K_vec_no_mean_ = new float[AMP_K];
        ArrayList<Float> codebook1 = m_codebookVQ.getCodebook(0).CodeBookArray();
        ArrayList<Float> codebook2 = m_codebookVQ.getCodebook(1).CodeBookArray();

        for (int k = 0; k < AMP_K; k++) {
            rate_K_vec_no_mean_[k] = (codebook1.get(AMP_K * indexes[0] + k)
                    + codebook2.get(AMP_K * indexes[1] + k));
        }

        post_filter_amp(rate_K_vec_no_mean_, 1.5f);

        float mean_ = CODES0[indexes[2]];

        for (int k = 0; k < AMP_K; k++) {
            rate_K_vec_[k] = rate_K_vec_no_mean_[k] + mean_;
        }
    }
    
    private int quantizeAmplitude(float mean) {
        int besti = 0;
        float beste = 1E32f;

        for (int j = 0; j < CODES0.length; j++) {
            float diff = CODES0[j] - mean;
            float e = (diff * diff);

            if (e < beste) {
                beste = e;
                besti = j;
            }
        }

        return besti;
    }

    private int encode_log_Wo(float Wo, int bits) {
        int Wo_levels = 1 << bits;

        float norm = (float) (Math.log10(Wo) - WO_MIN) / (WO_MAX - WO_MIN);
        int index = (int) Math.floor(Wo_levels * norm + 0.5f);

        if (index < 0) {
            index = 0;
        }

        if (index > (Wo_levels - 1)) {
            index = Wo_levels - 1;
        }

        return index;
    }

    private float decode_log_Wo(int index, int bits) {
        float step = (WO_MAX - WO_MIN) / (1 << bits);

        return (float) Math.pow(10.0, WO_MIN + step * index);
    }
}
