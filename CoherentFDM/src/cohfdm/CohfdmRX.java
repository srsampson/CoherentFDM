/*
 * Copyright (C) 2015 David Rowe
 *
 * All rights reserved
 *
 * Java Version by Steve Sampson
 */
package cohfdm;

import complex.Complex;
import complex.ComplexMath;

/**
 * Coherent FDM Modem Receiver Class
 */
public final class CohfdmRX implements IDefines {

    private final Complex[][] mrxSymbols;
    private final Complex[][] mcourseTimingSymbolBuffer;
    private final Complex[][] mcourseTimingFrameBuffer;
    private final Complex[][] mrxFilterMemory;
    private final Complex[][] mrxFilterMemoryTiming;
    //
    private final Complex[] mchFrameBuffer;
    private final Complex[] mprevRxSymbols;
    private final Complex[] mrxPhase;
    private final Complex[] mrxFrequency;
    //
    private Complex mbasebandRxPhase;
    //
    private final float[][] mpilots;
    //
    private float mfreqEstimate;
    private float mratio;
    private float msignalRMS;
    private float mnoiseRMS;
    private float msnrEstimate;
    private float mfrequencyOffsetFiltered;
    private float mrxTiming;
    private float mfrequencyFineEstimate;
    private float mfrequencySeparation;
    //
    private boolean msync;
    //
    private int msampleCenter;
    private int mnin;
    private int msyncTimer;
    private int mdiversityFactor;

    private class Point {   // Basically just a struct
        Complex slope;
        Complex intercept;
    }

    /**
     * Instantiate with Diversity True
     */
    public CohfdmRX() {
        this(true);
    }

    /**
     * Instantiate with Diversity option true = Diversity, false = No Diversity
     *
     * @param div boolean to enable or disable diversity
     */
    public CohfdmRX(boolean div) {
        int r, c, p;

        mdiversityFactor = (div) ? 2 : 1;

        mfrequencySeparation = 75.0f * 1.5f;    // 75 Hz + .5 Excess Bandwidth

        mchFrameBuffer = new Complex[NSW * NSYMROWPILOT * COHPSK_M];
        mrxSymbols = new Complex[NSYMROW][COHPSK_NC * mdiversityFactor];
        mcourseTimingSymbolBuffer = new Complex[NCT_SYMB_BUF][COHPSK_NC * mdiversityFactor];
        mcourseTimingFrameBuffer = new Complex[NSYMROWPILOT + 2][COHPSK_NC * mdiversityFactor];
        mrxFilterMemory = new Complex[COHPSK_NC * mdiversityFactor][COHPSK_NFILTER];
        mrxFilterMemoryTiming = new Complex[COHPSK_NC * mdiversityFactor][NT * P];
        mpilots = new float[NPILOTSFRAME * 2][COHPSK_NC];

        mprevRxSymbols = new Complex[COHPSK_NC * mdiversityFactor];
        mrxPhase = new Complex[COHPSK_NC * mdiversityFactor];
        mrxFrequency = new Complex[COHPSK_NC * mdiversityFactor];

        // We have to instantiate all of our complex arrays
        // initialization is a bit long-winded.
        for (c = 0; c < (NSW * NSYMROWPILOT * COHPSK_M); c++) {
            mchFrameBuffer[c] = new Complex();
        }

        for (c = 0; c < NSYMROW; c++) {
            for (r = 0; r < (COHPSK_NC * mdiversityFactor); r++) {
                mrxSymbols[c][r] = new Complex();
            }
        }

        for (c = 0; c < NCT_SYMB_BUF; c++) {
            for (r = 0; r < (COHPSK_NC * mdiversityFactor); r++) {
                mcourseTimingSymbolBuffer[c][r] = new Complex();
            }
        }

        for (c = 0; c < (NSYMROWPILOT + 2); c++) {
            for (r = 0; r < (COHPSK_NC * mdiversityFactor); r++) {
                mcourseTimingFrameBuffer[c][r] = new Complex();
            }
        }

        for (c = 0; c < (COHPSK_NC * mdiversityFactor); c++) {
            for (r = 0; r < COHPSK_NFILTER; r++) {
                mrxFilterMemory[c][r] = new Complex();
            }
        }

        for (c = 0; c < (COHPSK_NC * mdiversityFactor); c++) {
            for (r = 0; r < (NT * P); r++) {
                mrxFilterMemoryTiming[c][r] = new Complex();
            }
        }

        for (r = 0; r < (2 * NPILOTSFRAME);) {
            for (p = 0; p < NPILOTSFRAME; r++, p++) {
                for (c = 0; c < COHPSK_NC; c++) {
                    mpilots[r][c] = PILOTS[p][c];
                }
            }
        }

        Complex val = ComplexMath.cexp(new Complex());
        mbasebandRxPhase = val;
        mnin = COHPSK_M;

        for (c = 0; c < (COHPSK_NC * mdiversityFactor); c++) {
            mprevRxSymbols[c] = val;
            mrxPhase[c] = val;

            /*
             * Spread the frequencies out by different spreading values
             */
            float freqHz = mfrequencySeparation * (-(COHPSK_NC * mdiversityFactor) / 2.0f - 0.5f + (float) Math.pow(c + 1, 0.98));

            mrxFrequency[c] = ComplexMath.cexp(new Complex(0.0f, TAU * freqHz / COHPSK_FS));
        }
    }

    /**
     * Getter for modem SNR Estimate
     *
     * @return float SNR in dB of the received signal
     */
    public float getSNR() {
        float new_snr_est = 20.0f * (float) Math.log10((msignalRMS + 1E-6)
                / (mnoiseRMS + 1E-6f)) - 10.0f * (float) Math.log10(3000.0 / 700.0);

        msnrEstimate = 0.9f * msnrEstimate + 0.1f * new_snr_est;

        return msnrEstimate;
    }

    /**
     * Getter for Sync boolean
     *
     * True = Modem Sync, false = Modem Not Sync
     *
     * @return boolean for sync status
     */
    public boolean getSync() {
        return msync;
    }

    /**
     * Getter for Diversity Status
     *
     * True = Diversity Enabled, False = Diversity disabled
     *
     * @return boolean for diversity enabled or disabled
     */
    public boolean getDiversity() {
        return (mdiversityFactor == 2);
    }

    /**
     * Getter for number of symbols returned
     *
     * @return int for number of decoded symbols
     */
    public int getNIN() {
        return mnin;
    }

    /**
     * Getter for current Frequency Separation value
     *
     * @return float for the frequency separation value
     */
    public float getFreqSeparation() {
        return mfrequencySeparation;
    }

    /**
     * Method to demodulate a signal to baseband, and show sync state
     *
     * @param rx_bits a boolean array of the demodulated bits
     * @param sync_good a pointer to a boolean showing sync state
     * @param signal a complex array of the modulated signal
     * @param nin_frame an integer containing the frame size
     * @return a possibly updated frame size
     */
    public int demodulate(boolean[] rx_bits, boolean[] sync_good, Complex[] signal, int nin_frame) {
        Complex[][] ch_symb = new Complex[NSW * NSYMROWPILOT][COHPSK_NC * mdiversityFactor];
        boolean[] nextSync = new boolean[1];    // java pointer
        boolean[] anextSync = new boolean[1];   // java pointer
        int i, j, r, c;

        boolean lsync = msync;
        nextSync[0] = msync;

        for (i = 0; i < NSW * NSYMROWPILOT * COHPSK_M - nin_frame; i++) {
            mchFrameBuffer[i] = mchFrameBuffer[i + nin_frame];
        }

        for (j = 0; i < NSW * NSYMROWPILOT * COHPSK_M; i++, j++) {
            mchFrameBuffer[i] = signal[j];      // 0..599
        }

        if (lsync == false) {
            float maxRatio = 0.0f;
            float estimate = 0.0f;

            for (mfreqEstimate = COHPSK_CENTER - 40.0f; mfreqEstimate <= COHPSK_CENTER + 40.0f; mfreqEstimate += 40.0f) {
                receiveProcessor(ch_symb, mchFrameBuffer, NSW * NSYMROWPILOT, COHPSK_M, false);

                for (i = 0; i < NSW - 1; i++) {
                    updateCoarseTimingSymbolBuffer(ch_symb, (i * NSYMROWPILOT));
                }

                frameSyncFineFreqEstimate(ch_symb, ((NSW - 1) * NSYMROWPILOT), lsync, anextSync);

                if (anextSync[0] == true) {
                    if (mratio > maxRatio) {
                        maxRatio = mratio;
                        estimate = (mfreqEstimate - mfrequencyFineEstimate);
                        nextSync[0] = anextSync[0];
                    }
                }
            }

            if (nextSync[0] == true) {
                mfreqEstimate = estimate;
                receiveProcessor(ch_symb, mchFrameBuffer, NSW * NSYMROWPILOT, COHPSK_M, false);

                for (i = 0; i < NSW - 1; i++) {
                    updateCoarseTimingSymbolBuffer(ch_symb, (i * NSYMROWPILOT));
                }

                frameSyncFineFreqEstimate(ch_symb, ((NSW - 1) * NSYMROWPILOT), lsync, nextSync);

                if (Math.abs(mfrequencyFineEstimate) > 2.0f) {
                    nextSync[0] = false;
                }
            }

            if (nextSync[0] == true) {
                for (r = 0; r < NSYMROWPILOT + 2; r++) {
                    for (c = 0; c < COHPSK_NC * mdiversityFactor; c++) {
                        mcourseTimingFrameBuffer[r][c] = mcourseTimingSymbolBuffer[msampleCenter + r][c];
                    }
                }
            }
        } else if (lsync == true) {
            receiveProcessor(ch_symb, signal, NSYMROWPILOT, mnin, true);
            frameSyncFineFreqEstimate(ch_symb, 0, lsync, nextSync);

            for (r = 0; r < 2; r++) {
                for (c = 0; c < COHPSK_NC * mdiversityFactor; c++) {
                    mcourseTimingFrameBuffer[r][c] = mcourseTimingFrameBuffer[r + NSYMROWPILOT][c];
                }
            }

            for (r = 2; r < NSYMROWPILOT + 2; r++) {
                for (c = 0; c < COHPSK_NC * mdiversityFactor; c++) {
                    mcourseTimingFrameBuffer[r][c] = mcourseTimingSymbolBuffer[msampleCenter + r][c];
                }
            }
        }

        sync_good[0] = false;

        if ((nextSync[0] == true) || (lsync == true)) {
            symbolsToBits(rx_bits);
            sync_good[0] = true;
        }

        syncStateMachine(lsync, nextSync);
        msync = nextSync[0];
        lsync = nextSync[0];

        int lnin = COHPSK_M;

        // COHPSK_M / P = (100 / 4) = 25
        // The chances of g_rx_timing to be greater than this is probably nil
        if (lsync == true) {
            if (mrxTiming > COHPSK_M / P) {
                lnin = COHPSK_M + COHPSK_M / P;
            } else if (mrxTiming < -COHPSK_M / P) {
                lnin = COHPSK_M - COHPSK_M / P;
            }
        }

        mnin = lnin;

        return (NSYMROWPILOT - 1) * COHPSK_M + lnin;     // 600 normally
    }

    private void symbolsToBits(boolean[] bitPairs) {
        Complex[] y = new Complex[NPILOTSFRAME + 2];
        Complex[] rxSymbolLinear = new Complex[NSYMROW * COHPSK_NC * mdiversityFactor];
        float[][] pskPhase = new float[NSYMROW][COHPSK_NC * mdiversityFactor];

        for (int column = 0; column < (COHPSK_NC * mdiversityFactor); column++) {
            for (int pilot = 0; pilot < NPILOTSFRAME + 2; pilot++) {
                y[pilot] = ComplexMath.times(mcourseTimingFrameBuffer[SAMPLINGPOINTS[pilot]][column], mpilots[pilot][column % COHPSK_NC]);
            }

            final Point point = linearRegression(y);

            for (int row = 0; row < NSYMROW; row++) {     // 4 Data Rows
                Complex yfit = ComplexMath.add(point.intercept, ComplexMath.times(point.slope, row + NPILOTSFRAME)); // m*x+b
                pskPhase[row][column] = ComplexMath.carg(yfit);
            }
        }

        for (int column = 0; column < (COHPSK_NC * mdiversityFactor); column++) {
            for (int row = 0; row < NSYMROW; row++) {
                Complex phi_rect = ComplexMath.conjugate(ComplexMath.cexp(new Complex(0.0f, pskPhase[row][column])));
                mrxSymbols[row][column] = ComplexMath.times(mcourseTimingFrameBuffer[NPILOTSFRAME + row][column], phi_rect);

                rxSymbolLinear[column * NSYMROW + row] = mrxSymbols[row][column];
            }
        }

        // load the 56 bits detected from the received symbols
        // these are the 4 rows of 7 bit-pairs in the modem frame
        for (int column = 0; column < COHPSK_NC; column++) {
            for (int row = 0; row < NSYMROW; row++) {
                Complex symbol = mrxSymbols[row][column];

                // add in the diversity channel symbol
                for (int d = 1; d < mdiversityFactor; d++) {
                    symbol = ComplexMath.add(symbol, mrxSymbols[row][column + COHPSK_NC * d]);
                }

                Complex rot = ComplexMath.times(symbol, PI4); // rotate by 45 degrees
                int i = column * NSYMROW + row;

                bitPairs[2 * i + 1] = rot.real() < 0.0f;
                bitPairs[2 * i] = rot.imag() < 0.0f;
            }
        }

        float mag = 0.0f;

        for (int i = 0; i < (NSYMROW * COHPSK_NC * mdiversityFactor); i++) {
            mag += ComplexMath.cabs(rxSymbolLinear[i]);
        }

        msignalRMS = mag / (NSYMROW * COHPSK_NC * mdiversityFactor);

        float sum_x = 0.0f;
        float sum_xx = 0.0f;
        int n = 0;

        for (int i = 0; i < (NSYMROW * COHPSK_NC * mdiversityFactor); i++) {
            Complex s = rxSymbolLinear[i];

            if (Math.abs(s.real()) > msignalRMS) {
                sum_x += s.imag();
                sum_xx += (s.imag() * s.imag());
                n++;
            }
        }

        if (n > 1) {
            mnoiseRMS = (float) Math.sqrt((n * sum_xx - sum_x * sum_x) / (n * (n - 1)));
        } else {
            mnoiseRMS = 0.0f;
        }
    }

    private void downconvert(Complex[][] filtered, Complex[] signal, int lnin) {
        Complex[][] baseband = new Complex[COHPSK_NC * mdiversityFactor][COHPSK_M + COHPSK_M / P];

        /*
         * Mix the incoming signal and get the baseband sum+difference
         * It will need filtering (below) to get rid of the harmonic sum
         */
        for (int column = 0; column < (COHPSK_NC * mdiversityFactor); column++) { // 100 samples
            for (int i = 0; i < lnin; i++) {
                mrxPhase[column] = ComplexMath.times(mrxPhase[column], mrxFrequency[column]);
                baseband[column][i] = ComplexMath.timesConjugate(signal[i], mrxPhase[column]);
            }
        }

        // Normalize carriers to combat drift
        for (int column = 0; column < (COHPSK_NC * mdiversityFactor); column++) {
            mrxPhase[column] = ComplexMath.normalize(mrxPhase[column]);
        }

        // Now Raised Root Cosine Filter (RRC) the signal

        int n = COHPSK_M / P;

        for (int i = 0, j = 0; i < lnin; i += n, j++) {
            for (int column = 0; column < (COHPSK_NC * mdiversityFactor); column++) {
                for (int k = COHPSK_NFILTER - n, l = i; k < COHPSK_NFILTER; k++, l++) {
                    mrxFilterMemory[column][k] = baseband[column][l];
                }
            }

            for (int column = 0; column < (COHPSK_NC * mdiversityFactor); column++) {
                filtered[column][j] = new Complex();

                for (int k = 0; k < COHPSK_NFILTER; k++) {
                    filtered[column][j] = ComplexMath.add(filtered[column][j],
                            ComplexMath.times(mrxFilterMemory[column][k], GTALPHA5ROOT[k]));
                }
            }

            for (int column = 0; column < (COHPSK_NC * mdiversityFactor); column++) {
                for (int k = 0, l = n; k < COHPSK_NFILTER - n; k++, l++) {
                    mrxFilterMemory[column][k] = mrxFilterMemory[column][l];
                }
            }
        }
    }

    private void frequencyShift(Complex[] offsetSignal, Complex[] signal, int offset, int lnin) {
        Complex rxPhase = ComplexMath.cexp(new Complex(0.0f, TAU * -mfreqEstimate / COHPSK_FS));

        for (int i = 0; i < lnin; i++) {
            mbasebandRxPhase = ComplexMath.times(mbasebandRxPhase, rxPhase);
            offsetSignal[i] = ComplexMath.times(signal[offset + i], mbasebandRxPhase);
        }

        mbasebandRxPhase = ComplexMath.normalize(mbasebandRxPhase);
    }

    private void receiveProcessor(Complex[][] symbols, Complex[] signal,
            int nsymb, int ninval, boolean freqTrack) {
        Complex[] offsetSignal = new Complex[COHPSK_M + COHPSK_M / P];
        Complex[][] rxFiltered = new Complex[COHPSK_NC * mdiversityFactor][P + 1];
        Complex[] rxOneSymbol = new Complex[COHPSK_NC * mdiversityFactor];

        int lnin = ninval;
        int index = 0;
        float adjustedRxTiming = 0.0f;

        for (int row = 0; row < nsymb; row++) {      // nsymb = 6 when in sync, else 24
            frequencyShift(offsetSignal, signal, index, lnin);
            index += lnin;                       // 0, 100, 200, 300, 400, 500  (0..599) in sync
                                                 // else (0..2399) when un-sync (100 * 24)
            downconvert(rxFiltered, offsetSignal, lnin);

            adjustedRxTiming = rxEstimatedTiming(rxOneSymbol, rxFiltered, lnin);

            System.arraycopy(rxOneSymbol, 0, symbols[row], 0, COHPSK_NC * mdiversityFactor);

            if (freqTrack == true) {
                /*
                 * Derive carrier frequency of each column using 4th power
                 */
                Complex modStrip = new Complex();

                for (int column = 0; column < (COHPSK_NC * mdiversityFactor); column++) {
                    Complex adiff = ComplexMath.timesConjugate(rxOneSymbol[column], mprevRxSymbols[column]);
                    mprevRxSymbols[column] = rxOneSymbol[column];

                    /*
                     * 4th power strips QPSK modulation, by multiplying phase by 4
                     * Using the abs value of the real coord was found to help
                     * non-linear issues when noise power was large.
                     */

                    Complex amodStrip = ComplexMath.times(adiff, adiff);
                    amodStrip = ComplexMath.times(amodStrip, amodStrip);  // adiff^4
                    amodStrip = new Complex(Math.abs(amodStrip.real()), amodStrip.imag());

                    modStrip = ComplexMath.add(modStrip, amodStrip);
                }

                mfrequencyOffsetFiltered = (1.0f - 0.005f) * mfrequencyOffsetFiltered + 0.005f
                        * ComplexMath.carg(modStrip);
                mfreqEstimate += (0.2f * mfrequencyOffsetFiltered);
            }

            if (lnin != COHPSK_M) {
                lnin = COHPSK_M;
            }
        }

        mrxTiming = adjustedRxTiming;
    }

    private float pilotCorrelation(float[] corr_out, int t, float f_fine) {
        Complex[] freqFinePhase = new Complex[NPILOTSFRAME + 2];

        // Move trig out of main loop
        for (int pilot = 0; pilot < NPILOTSFRAME + 2; pilot++) {
            freqFinePhase[pilot] = ComplexMath.cexp(new Complex(0.0f, f_fine * TAU * (SAMPLINGPOINTS[pilot] + 1.0f) / COHPSK_RS));
        }

        float corr = 0.0f;
        float mag = 1E-12f; // avoid NAN

        for (int column = 0; column < (COHPSK_NC * mdiversityFactor); column++) {
            Complex acorr = new Complex();
            int pc = column % COHPSK_NC;

            for (int pilot = 0; pilot < NPILOTSFRAME + 2; pilot++) {
                Complex freqCorr = ComplexMath.times(mcourseTimingSymbolBuffer[t + SAMPLINGPOINTS[pilot]][column], freqFinePhase[pilot]);
                acorr = ComplexMath.add(acorr, ComplexMath.times(freqCorr, mpilots[pilot][pc]));
                mag += ComplexMath.cabs(freqCorr);
            }

            corr += ComplexMath.cabs(acorr);
        }

        corr_out[0] = corr; // java pointer

        return mag;
    }

    private void frameSyncFineFreqEstimate(Complex[][] ch_symb, int offset, boolean sync, boolean[] nextSync) {
        updateCoarseTimingSymbolBuffer(ch_symb, offset);

        if (sync == false) {
            float[] corr = new float[1];
            float max_corr = 0.0f;
            float max_mag = 1E-12f; // avoid NAN

            for (int j = -80; j <= 80; j++) {
                float f_fine = (float) j * .25f;

                for (int i = 0; i < NSYMROWPILOT; i++) {
                    float mag = pilotCorrelation(corr, i, f_fine);

                    if (corr[0] >= max_corr) {
                        max_corr = corr[0];
                        max_mag = mag;
                        msampleCenter = i;
                        mfrequencyFineEstimate = f_fine;
                    }
                }
            }

            if (max_corr / max_mag > 0.9f) {
                msyncTimer = 0;
                nextSync[0] = true;
            } else {
                nextSync[0] = false;
            }

            mratio = max_corr / max_mag;
        }
    }

    private void updateCoarseTimingSymbolBuffer(Complex[][] symbol, int offset) {
        for (int row = 0; row < NCT_SYMB_BUF - NSYMROWPILOT; row++) {
            System.arraycopy(mcourseTimingSymbolBuffer[row + NSYMROWPILOT], 0, mcourseTimingSymbolBuffer[row], 0, COHPSK_NC * mdiversityFactor);
        }

        for (int row = NCT_SYMB_BUF - NSYMROWPILOT, i = 0; row < NCT_SYMB_BUF; row++, i++) {
            System.arraycopy(symbol[offset + i], 0, mcourseTimingSymbolBuffer[row], 0, COHPSK_NC * mdiversityFactor);
        }
    }

    private void syncStateMachine(boolean lsync, boolean[] nextSync) {
        if (lsync == true) {
            float[] corr = new float[1];    // java pointer

            float mag = pilotCorrelation(corr, msampleCenter, mfrequencyFineEstimate);
            mratio = Math.abs(corr[0]) / mag;

            if (mratio < 0.8f) {
                msyncTimer++;
            } else {
                msyncTimer = 0;  // validate sync
            }

            // Give up and set sync to false
            if (msyncTimer == 10) {
                nextSync[0] = false;
            }
        }
    }

    private float rxEstimatedTiming(Complex[] symbols, Complex[][] rxFiltered, int lnin) {
        int adjust = P - lnin * P / COHPSK_M;

        // update buffer of NT rate P filtered symbols
        for (int column = 0; column < (COHPSK_NC * mdiversityFactor); column++) {
            for (int i = 0, j = P - adjust; i < (NT - 1) * P + adjust; i++, j++) {
                mrxFilterMemoryTiming[column][i] = mrxFilterMemoryTiming[column][j];
            }
        }

        for (int column = 0; column < (COHPSK_NC * mdiversityFactor); column++) {
            for (int i = (NT - 1) * P + adjust, j = 0; i < NT * P; i++, j++) {
                mrxFilterMemoryTiming[column][i] = rxFiltered[column][j];
            }
        }

        // sum envelopes of all carriers

        float[] env = new float[NT * P];

        for (int i = 0; i < (NT * P); i++) {
            env[i] = 0.0f;

            for (int column = 0; column < (COHPSK_NC * mdiversityFactor); column++) {
                env[i] += ComplexMath.cabs(mrxFilterMemoryTiming[column][i]);
            }
        }

        /*
         * The envelope has a frequency component at the symbol rate. The phase
         * of this frequency component indicates the timing. So work out single
         * DFT at frequency 2*pi/P
         */
        Complex x = new Complex();
        Complex phase = ComplexMath.cexp(new Complex());

        for (int i = 0; i < (NT * P); i++) {
            x = ComplexMath.add(x, ComplexMath.times(phase, env[i]));
            phase = ComplexMath.times(phase, PPHASE);
        }

        /*
         * Map phase to estimated optimum timing instant at rate P. The P/4 part
         * was adjusted by experiment, I know not why....
         */
        float adjustedRxTiming = ComplexMath.carg(x) / TAU;
        float rx_timing = (adjustedRxTiming * P) + (P / 4);

        if (rx_timing > P) {
            rx_timing -= P;
        } else if (rx_timing < -P) {
            rx_timing += P;
        }

        /*
         * rx_filter_mem_timing contains Nt*P samples (Nt symbols at rate P),
         * where Nt is odd. Lets use linear interpolation to resample in the
         * center of the timing estimation window
         */
        rx_timing += Math.floor(NT / 2.0) * P;
        int low_sample = (int) Math.floor(rx_timing);

        float fract = rx_timing - low_sample;
        int high_sample = (int) Math.ceil(rx_timing);

        for (int column = 0; column < (COHPSK_NC * mdiversityFactor); column++) {
            symbols[column] = ComplexMath.add(ComplexMath.times(mrxFilterMemoryTiming[column][low_sample - 1], 1.0f - fract),
                    ComplexMath.times(mrxFilterMemoryTiming[column][high_sample - 1], fract)
            );
        }

        /*
         * This value will be +/- half a symbol
         * so will wrap around at +/- M/2 or
         * +/- 80 samples with M=160
         */
        return adjustedRxTiming * COHPSK_M;
    }

    private Point linearRegression(Complex[] y) {
        Complex sumxy = new Complex();
        Complex sumy = new Complex();
        float sumx = 0.0f;
        float sumx2 = 0.0f;

        for (int i = 0; i < SAMPLINGPOINTS.length; i++) {
            float x = (float) SAMPLINGPOINTS[i];
            
            sumx += x;
            sumx2 += (x * x);
            sumxy = ComplexMath.add(sumxy, ComplexMath.times(y[i], x));
            sumy = ComplexMath.add(sumy, y[i]);
        }

        float denom = (SAMPLINGPOINTS.length * sumx2 - sumx * sumx);

        Point point = new Point();

        // fits y = mx + b to the (x,y) data
        // x is the sampling points

        if (denom == 0.0f) { // never happen with the data used in this modem
            // no solution
            point.slope = new Complex();
            point.intercept = new Complex();
        } else {
            point.slope = ComplexMath.divide(ComplexMath.add(ComplexMath.times(sumxy, SAMPLINGPOINTS.length), ComplexMath.cneg(ComplexMath.times(sumy, sumx))), denom);
            point.intercept = ComplexMath.divide(ComplexMath.add(ComplexMath.times(sumy, sumx2), ComplexMath.cneg(ComplexMath.times(sumxy, sumx))), denom);
        }

        return point;
    }
}
